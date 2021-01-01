process SketchyStream {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.tsv"
    publishDir "${params.outdir}/stream/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.sssh.tsv"

    input:
    tuple val(id), file(fx)
    each file(db)
    each read_limit

    output:
    file("${id}.tsv")

    """
    sketchy stream --fastx $fx --db $db --reads $read_limit --ranks $params.ranks --stability $params.stability --threads $task.cpus > ${id}.sssh.tsv
    cat  ${id}.sssh.tsv | sketchy predict --db $db --limit $params.limit > predict.tsv
    tail -$limit predict.tsv > ${id}.tsv
    """

}

process SketchyScreen {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/screen/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    each file(db)
    each read_limit

    output:
    file("${id}.tsv")

    script:

    _read_limit = 4*read_limit

    """
    head -$_read_limit $fx | sketchy screen --fastx - --db $db --limit $params.limit --threads $task.cpus > ${id}.tsv
    """

}

process SketchyDist {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/dist/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), file(fx)
    each file(db)
    each read_limit

    output:
    file("${id}.tsv")

    script:

    _read_limit = 4*read_limit

    """
    head -$_read_limit $fx | sketchy dist --fastx - --db $db --limit $params.limit --threads $task.cpus > ${id}.tsv
    """

}