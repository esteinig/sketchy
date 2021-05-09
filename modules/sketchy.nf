process SketchyStream {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}_${rep}.tsv"
    publishDir "${params.outdir}/stream/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}_${rep}.sssh"
    publishDir "${params.outdir}/stream/${db.baseName}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), val(rep), val(read_limit), file(fx)
    each file(db)

    output:
    file("${id}_${rep}.tsv")
    file("${id}_${rep}.sssh")
    file("header.txt")
    
    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy stream --db $db --ranks $params.ranks --stability $params.stability --threads $task.cpus > ${id}_${rep}.sssh
    cat  ${id}_${rep}.sssh | sketchy predict --db $db --limit $params.limit --genotype > predict.tsv
    tail -n $params.limit predict.tsv > ${id}_${rep}.tsv
    sketchy head --db $db > header.txt
    """

}

process SketchyScreen {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/screen/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}_${rep}.tsv"
    publishDir "${params.outdir}/screen/${db.baseName}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), val(rep), val(read_limit), file(fx)
    each file(db)

    output:
    file("${id}_${rep}.tsv")
    file("header.txt")
    
    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy screen --fastx - --db $db --limit $params.limit --threads $task.cpus > ${id}_${rep}.tsv
    sketchy head --db $db > header.txt
    """

}


process SketchyScreenWinner {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/screen_winner/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}_${rep}.tsv"
    publishDir "${params.outdir}/screen_winner/${db.baseName}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), val(rep), val(read_limit), file(fx)
    each file(db)

    output:
    file("${id}_${rep}.tsv")
    file("header.txt")
    
    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy screen --fastx - --db $db --limit $params.limit --threads $task.cpus --winner > ${id}_${rep}.tsv
    sketchy head --db $db > header.txt
    """

}

process SketchyDist {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/dist/${db.baseName}/${read_limit}", mode: "copy", pattern: "${id}_${rep}.tsv"
    publishDir "${params.outdir}/dist/${db.baseName}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), val(rep), val(read_limit), file(fx)
    each file(db)

    output:
    file("${id}_${rep}.tsv")
    file("header.txt")

    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy dist --fastx - --db $db --limit $params.limit --threads $task.cpus > ${id}_${rep}.tsv
    sketchy head --db $db > header.txt
    """

}

process Bootstrap {

    label "sketchy"
    tag { id }

    input:
    tuple val(id), file(fx)
    each rep
    each read_limit

    output:
    tuple val(id), val(rep), val(read_limit), file("${id}_${read_limit}_${rep}.fastq")

    script:

    """
    sketchy-utils fastx sample --fastx $fx --output ${id}_${read_limit}_${rep}.fastq --sample $read_limit --replacement
    """

}


process ReadLimit {

    label "sketchy"
    tag { id }

    input:
    tuple val(id), file(fx)
    each read_limit

    output:
    tuple val(id), val(0), val(read_limit), file(fx)

    """
    echo "PLACEHOLDER"
    """

}
