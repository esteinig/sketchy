// Sketchy Nextflow managing switchable tasks that require mo4e compute power (like bootstrapping)

log.info """
.
.
sketchy.nf v0.4a
.
fastq          =  $params.fastq
conda          =  $params.conda
resources      =  $params.resources
.
.
"""

// Bootstrap pipeline:

//
// reads = Channel.fromPath(params.fastq)
//      .splitFastq( by: 10, file: true )
//      .println()



import groovy.json.JsonSlurper

if (params.sketch_files){

  if (params.uuid){
    uuid = '--uuid'
  } else {
    uuid = ''
  }

  sketchies = Channel
      .fromPath(params.sketch_files)
      .map { file -> tuple(file.baseName, file) }


  process LinkSketchy {

      label "sketchy"

      input:
      set id, file(sketch_file) from sketchies

      output:
      set id, file("${id}.out") into sketch

      """
      sketchy link --data $sketch_file --outdir ${id}.out $uuid
      """
  }

  process SketchSketchy {

      label "sketchy"

      publishDir "${params.outdir}/", mode: "copy"

      input:
      set id, file(fasta_dir) from sketch
      each k from params.kmers
      each hash_size from params.hash_sizes

      output:
      file("*.msh")

      """
      sketchy sketch -f $fasta_dir --kmer $k --size $hash_size --prefix ${id}
      """
  }

}



if (params.bootstrap) {

fastq = Channel
        .fromPath(params.fastq)
        .map { file -> tuple(file.baseName, file.baseName.split('\\.')[0], file) }


process Bootstrap {

    // Create replicate shuffled data (--sample 1.0) with replacement (-r)
    // Predict on bootstrap top 10 (default)
    // Determine breakpoints at first 300 continous reads (--stable)

    // NB: plots secondary, primary summary are breakpoints

    label "sketchy"
    publishDir "${params.outdir}/bootstrap", mode: "copy"

    input:
    set id, sketch, file(fq) from fastq
    each bs from Channel.from(1..params.nboot)

    output:
    file("${id}.${bs}.boot.tsv")
    file("${id}.${bs}.bp.tsv")
    file("${id}.${bs}.png")

    """
    sketchy fq-sample -f $fq -o ${id}.bs.fq -r -s $params.reads
    sketchy predict -f ${id}.bs.fq -s $sketch -t $task.cpus --top $params.top \
    -o ${id}.${bs}.boot -r 0
    sketchy plot -d ${id}.${bs}.boot.tsv -b --stable $params.stable $params.boot_plot -t 5 -p ${id}.${bs} -f png
    """
}
}


if (params.mixture) {

  fastq = Channel
      .fromPath(params.fastq)
      .map { file -> tuple(file.baseName, file) }

  // Classify and extract S. aureus reads from Zymo Mocks

  process Kraken2 {

      label "kraken2"

      publishDir "${params.outdir}/kraken", mode: "copy"

      input:
      set id, file(fastq) from fastq

      output:
      set id, file("$fastq"), file("${id}.out") into kraken_filter

      """
      kraken2 --db ${params.resources}/${params.minikraken} --threads $task.cpus --output "${id}.out" \
      --gzip-compressed --report ${id}.report --use-names ${fastq}
      """
    }

    process KrakenFilter {

        label "sketchy"

        publishDir "${params.outdir}/species", mode: "copy"

        input:
        set id, file(fastq), file(kraken) from kraken_filter

        output:
        set file("${id}.species.reads.out"), file("${id}.species.fq"), file("$kraken"), file("reads.fq")

        """
        zcat $fastq > reads.fq
        cat $kraken | grep "$params.species" > ${id}.species.reads.out
        sketchy select-fastq -f reads.fq -i ${id}.species.reads.out -o ${id}.species.fq
        """
      }

}
