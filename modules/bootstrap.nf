process BootstrapBuild {


    tag { "$sample - $replicate" }

    label "bbuild"

    publishDir "${params.outdir}/bootstrap"

    input:
    set file(fasta_index), file(fasta_directory)
    each sample from params.samples
    each replicate from params.bootstraps

    output:
    file(replicate)
    set sample, replicate, file("sketch_${sample}_${replicate}") into sketchy_bootstrap

    """
    sketchy-utils database bootstrap --fasta_index $fasta_index --fasta_directory $fasta_directory --sample $sample --output replicate${replicate}
    mash sketch -k 15 -s 1000 replicate_${replicate}/*.fasta -o $replicate
    sketchy-utils database create 
    """

}