process Sketch {

    tag { "Reference sketch ($prefix): k=$kmer_size s=$sketch_size" }
    label "sketch"

    memory { params.sketch_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/sketches", mode: "symlink", pattern: "*.msh"

    input:
    tuple val(prefix), val(fasta_glob), file(genotype_file), file(fasta_files)  // collected reference genomes, with glob to list into stdin
    each kmer_size
    each sketch_size

    output:
    tuple val(prefix), file("${prefix}_k${kmer_size}_s${sketch_size}.msh")

    script:

    """
    ls $fasta_glob | sketchy sketch -k $kmer_size -s $sketch_size -o ${prefix}_k${kmer_size}_s${sketch_size}.msh


    """ 

}