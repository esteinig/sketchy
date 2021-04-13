process SketchyStream {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.tsv"
    publishDir "${params.outdir}/stream/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.sssh"
    publishDir "${params.outdir}/stream/${db.baseName}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), file(fx)
    each file(db)
    each read_limit

    output:
    file("${id}.tsv")
    file("${id}.sssh")
    file("header.txt")
    
    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy stream --db $db --ranks $params.ranks --stability $params.stability --threads $task.cpus > ${id}.sssh
    cat  ${id}.sssh | sketchy predict --db $db --limit $params.limit > predict.tsv
    tail -$params.limit predict.tsv > ${id}.tsv
    sketchy head --db $db > header.txt
    """

}

process SketchyScreen {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/screen/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.tsv"
    publishDir "${params.outdir}/screen/${db.baseName}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), file(fx)
    each file(db)
    each read_limit

    output:
    file("${id}.tsv")
    file("header.txt")
    
    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy screen --fastx - --db $db --limit $params.limit --threads $task.cpus > ${id}.tsv
    sketchy head --db $db > header.txt
    """

}

process SketchyDist {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/dist/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}.tsv"
    publishDir "${params.outdir}/dist/${db.baseName}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), file(fx)
    each file(db)
    each read_limit

    output:
    file("${id}.tsv")
    file("header.txt")

    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy dist --fastx - --db $db --limit $params.limit --threads $task.cpus > ${id}.tsv
    sketchy head --db $db > header.txt
    """

}