// Sketchy Nextflow managing switchable tasks that require mo4e compute power (like bootstrapping)

log.info """
=======================================================
                      Sketchy
                  Development v0.1
=======================================================

fastq                   =           $params.fastq
sketch                  =           $params.sketch
sketch_data             =           $params.sketch_data
bootstraps              =           $params.bootstraps
sample_reads            =           $params.sample_reads
predict_reads           =           $params.predict_reads
prefix                  =           $params.prefix
outdir                  =           $params.outdir

=======================================================
=======================================================
"""

fastq = Channel
		.fromPath(params.fastq)
		.map { file -> tuple(file.baseName, file) }


process createBootstraps {

    label "sketchy"
    label "bootstrap"

    tag { "Sampling from shuffled reads with replacement." }

    input:
    set id, file(fastq) from fastq

    output:
    file("bootstraps/${id}/boot*.fastq") into predict_bootstrap mode flatten
    val id into predict_id

    """
    sketchy boot --fastq $fastq -b $params.bootstraps -r $params.sample_reads -o bootstraps
    """
}

process predictBootstrap {

        label "sketchy"
        label "predict"

        tag { "Predicting sketchy scores on bootstrap replicate." }
        publishDir "params.outdir/${id}/", mode: "copy"

        input:
        file(replicate) from predict_bootstrap
        val id from predict_id

        output:
        file("${replicate}.tab")

        """
        sketchy pboot -f $replicate -s $params.sketch -d $params.sketch_data -r $params.predict_reads \
        -c $task.cpus -o . -p $replicate
        """
    }


