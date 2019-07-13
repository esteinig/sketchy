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

sketch = file(params.sketch)
sketch_data = file(params.sketch_data)

// Bootstrap pipeline:

fastq = Channel
    .fromPath(params.fastq)
    .map { file -> tuple(file.baseName, file) }


process Bootstrapping {

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


process Bootstrap {

    label "sketchy"
    label "predict"

    publishDir "${params.outdir}/${id}", mode: "copy"

    input:
    set val(id), file(replica) from predict_id.combine(predict_bootstrap)

    output:
    file("${replica.baseName}.tab")

    """
    sketchy pboot -f $replica -s $sketch -d $sketch_data -r $params.predict_reads \
    -c $task.cpus -o . -p ${replica.baseName}
    """

}




