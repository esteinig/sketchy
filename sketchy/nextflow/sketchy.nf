// Sketchy Nextflow managing switchable tasks that require mo4e compute power (like bootstrapping)

log.info """
=======================================================
                      Sketchy
                  Development v0.1
=======================================================

fastq                   =           $params.fastq

sketchy                 =           $params.sketchy
mixture                 =           $params.mixture

prefix                  =           $params.prefix
outdir                  =           $params.outdir

reads                   =           $params.reads
sketch                  =           $params.sketch
index                   =           $params.index
lineage                 =           $params.true_lineage

conda                   =           $params.conda
resources               =           $params.resources
=======================================================
=======================================================
"""

// Bootstrap pipeline:

fastq = Channel
    .fromPath(params.fastq)
    .map { file -> tuple(file.baseName, file) }
//
// reads = Channel.fromPath(params.fastq)
//      .splitFastq( by: 10, file: true )
//      .println()

if (params.sketchy) {
  process Sketchy {

      label "sketchy"

      tag { "Sketchy predictions on $id and $read reads" }

      publishDir "${params.outdir}/sketchy", mode: "copy"

      input:
      set id, file(fastq) from fastq
      each read from params.reads

      output:
      file("${id}_${read}")
      file("${id}_${read}_evaluation")

      """
      sketchy predict -s $params.sketch -d $params.index -f $fastq \
        -r $read --keep --tmp ${id}_${read} \
        --output ${id}_${read}.out \
        --mode $params.mode --cores $task.cpus

      sketchy evaluate --indir ${id}_${read} \
        --limit $read \
        --color green \
        --outdir ${id}_${read}_evaluation \
        --lineage $params.true_lineage
      """
  }
}

if (params.mixture) {

  // Classify and extract S. aureus reads from Zymo Mocks

  process Kraken2 {

      label "kraken2"

      tag { "Kraken2 for Zymo mock community: $id" }

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

        tag { "Sketchy filter for Zymo mock community: $id" }

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
