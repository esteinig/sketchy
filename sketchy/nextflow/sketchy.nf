// Sketchy Nextflow managing switchable tasks that require mo4e compute power (like bootstrapping)

log.info """
--------------------------------------------------------------------------------
fastq          =  $params.fastq
meta           =  $params.meta
conda          =  $params.conda
resources      =  $params.resources

prefix         =  $params.prefix
outdir         =  $params.outdir

preprint       =  $params.preprint
mixture        =  $params.mixture
reads          =  $params.reads
evaluation     =  $params.evaluation
ranks          =  $params.ranks
--------------------------------------------------------------------------------
"""

// Bootstrap pipeline:

//
// reads = Channel.fromPath(params.fastq)
//      .splitFastq( by: 10, file: true )
//      .println()



import groovy.json.JsonSlurper

if (params.preprint) {

  jsonSlurper = new JsonSlurper();

  meta = new File(params.meta)
  meta_json = meta.text
  meta_data = jsonSlurper.parseText(meta_json)

  fastq = Channel
      .fromPath(params.fastq)
      .map { file -> tuple(file.baseName, file) }

  process LineageCaller {

      label "sketchy"

      tag { "Sketchy predictions on $id and $read reads" }

      publishDir "${params.outdir}/${id}", mode: "copy"

      input:
      set id, file(fastq) from fastq

      output:
      set id, file("${id}_${params.reads}") into evaluate
      val id into concat2

      """
      sketchy predict -f $fastq -s ${meta_data[id].sketch} \
      -r $params.reads --mode $params.mode --tmp ${id}_${params.reads} --keep \
      --cores $task.cpus
      """
  }

  process Evaluation {

      label "sketchy"

      publishDir "${params.outdir}/${id}", mode: "copy"

      input:
      set id, file(indir) from evaluate
      each eval from params.evaluation

      output:
      file("${id}_${eval}_evaluation")
      file("${id}.${eval}.pdf") into concat

      """
      sketchy evaluate --indir $indir \
        --limit $eval --ranks $params.ranks \
        --primary ${meta_data[id].color} \
        --secondary ${meta_data[id].color} \
        --lineage ${meta_data[id].lineage} \
        --genotype ${meta_data[id].genotype} \
        --resistance ${meta_data[id].resistance} \
        --outdir ${id}_${eval}_evaluation

      mv ${id}_${eval}_evaluation/evaluation.pdf ${id}.${eval}.pdf
      """
  }

  process ConcatEvaluation {

      label "sketchy"

      publishDir "${params.outdir}/", mode: "copy"

      input:
      file('*.pdf') from concat.collect()
      val id from concat2

      output:
      file("${id}.evaluation.pdf")

      """
      sketchy concat --output ${id}.evaluations.pdf
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
