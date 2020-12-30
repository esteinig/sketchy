process SketchyStream {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    file(db)

    output:
    file("${id}.tsv")

    """
    sketchy-rs stream --fastx $fx --db $db --ranks $params.ranks --stability $params.stability --threads $task.cpus | sketchy-rs predict --db $db --limit $params.limit > ${id}.tsv
    """

}

process SketchyScreen {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    file(db)

    output:
    file("${id}.tsv")

    """
    sketchy-rs screen --fastx $fx --db $db --limit $params.limit --threads $task.cpu > ${id}.tsv
    """

}

process SketchyDist {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    file(db)

    output:
    file("${id}.tsv")

    """
    sketchy-rs dist --fastx $fx --db $db --limit $params.limit --threads $task.cpu > ${id}.tsv
    """

}