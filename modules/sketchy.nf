process SketchyStream {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream/${db.baseName}", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    each file(db)

    output:
    file("${id}.tsv")

    """
    sketchy stream --fastx $fx --db $db --reads $params.reads --ranks $params.ranks --stability $params.stability --threads $task.cpus | sketchy-rs predict --db $db --limit $params.limit > ${id}.tsv
    """

}

process SketchyScreen {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/screen/${db.baseName}", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    each file(db)

    output:
    file("${id}.tsv")

    """
    sketchy screen --fastx $fx --db $db --limit $params.limit --threads $task.cpu > ${id}.tsv
    """

}

process SketchyDist {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/dist/${db.baseName}", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    each file(db)

    output:
    file("${id}.tsv")

    """
    sketchy dist --fastx $fx --db $db --limit $params.limit --threads $task.cpu > ${id}.tsv
    """

}