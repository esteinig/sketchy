process BootstrapBuild {


    tag { "Bootstrap : $sample : $replicate" }

    label "bbuild"

    publishDir "${params.outdir}/bootstrap"

    input:
    file(fasta_directory)
    each sample
    each replicate

    output:
    tuple val(sample), val(replicate), file("s${sample}_r${replicate}")

    """
    sketchy-utils database bootstrap --fasta_directory $fasta_directory --bootstrap_sample $sample --outdir bootstrap_${replicate} --genotypes bootstrap_${replicate}
    mash sketch -k 15 -s 1000 -o bootstrap_$replicate replicate_${replicate}/*.fasta
    sketchy-utils database create --sketch bootstrap_${replicate}.msh --genotypes bootstrap_${replicate}.tsv $params.create_options --outdir s${sample}_r${replicate}
    """

}

process SketchyStream {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/${sample}/${read_limit}", mode: "copy", pattern: "${id}_${replicate}.tsv"
    publishDir "${params.outdir}/${sample}/${read_limit}", mode: "copy", pattern: "${id}_${replicate}.sssh"
    publishDir "${params.outdir}/${sample}/${read_limit}", mode: "copy", pattern: "header.txt"

    input:
    tuple val(id), val(read_limit), file(fx)
    tuple val(sample), val(replicate), file(db)

    output:
    file("${id}_${replicate}.tsv")
    file("${id}_${replicate}.sssh")
    file("header.txt")
    
    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy stream --db $db --ranks $params.ranks --stability $params.stability --threads $task.cpus > ${id}_${replicate}.sssh
    cat  ${id}_${replicate}.sssh | sketchy predict --db $db --limit $params.limit --genotype > predict.tsv
    tail -n $params.limit predict.tsv > ${id}_${replicate}.tsv
    sketchy head --db $db > header.txt
    """

}