process BootstrapBuild {


    tag { "Bootstrap : $sample : $replicate" }

    label "build"

    publishDir "${params.outdir}/bootstrap"

    input:
    file(fasta_directory)
    file(reference_database)
    each sample
    each replicate

    output:
    tuple val(sample), val(replicate), file("s${sample}_r${replicate}")

    """
    sketchy-utils database bootstrap --fasta_directory \$PWD/$fasta_directory --bootstrap_samples $sample --reference_database $reference_database --outdir \$PWD/bootstrap_${replicate} --genotypes bootstrap_${replicate}
    mash sketch -k 15 -s 1000 -o bootstrap_$replicate bootstrap_${replicate}/*.fasta
    sketchy-utils database create --sketch bootstrap_${replicate}.msh --genotypes bootstrap_${replicate}.tsv $params.create_options --outdir s${sample}_r${replicate}
    """

}

process SketchyStream {

    label "sketchy"
    tag { id }

    publishDir "${params.outdir}/stream/${sample}/${read_limit}", mode: "copy", pattern: "${id}_${replicate}.tsv"
    publishDir "${params.outdir}/stream/${sample}/${read_limit}", mode: "copy", pattern: "${id}_${replicate}.sssh"
    
    input:
    tuple val(id), val(read_limit), file(fx), val(sample), val(replicate), file(db)

    output:
    file("${id}_${replicate}.tsv")
    file("${id}_${replicate}.sssh")
    
    script:

    _read_limit = 4*read_limit

    """
    SKETCHY_PATH=\$PWD
    head -$_read_limit $fx | sketchy stream --db $db --ranks $params.ranks --stability $params.stability --threads $task.cpus > ${id}_${replicate}.sssh
    cat  ${id}_${replicate}.sssh | sketchy predict --db $db --limit $params.limit --genotype > predict.tsv
    tail -n $params.limit predict.tsv > ${id}_${replicate}.tsv
    """

}

process PublishHeader {

    label "sketchy"
    tag { sample }

    input:
    file(reference_database)
    each sample

    publishDir "${params.outdir}/stream/${sample}/", mode: "copy", pattern: "header.txt"

    """
    SKETCHY_PATH=\$PWD
    sketchy head --db $reference_database > header.txt
    """
}