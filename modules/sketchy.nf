process Sketch {

    tag { "prefix=$prefix k=$kmer_size s=$sketch_size" }
    label "sketch"

    memory { params.sketch_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/sketches", mode: "symlink", pattern: "*.msh"

    input:
    val(prefix)
    val(fasta_glob) 
    file(genotype_file)
    file(fasta_files)  // collected reference genomes, with glob to list into stdin
    each kmer_size
    each sketch_size

    output:
    tuple val("${prefix}_k${kmer_size}_s${sketch_size}"), file("${prefix}_k${kmer_size}_s${sketch_size}.msh")

    script:

    sketchy = params.exec ?: "sketchy"

    """
    find . -name "$fasta_glob" | $sketchy sketch -k $kmer_size -s $sketch_size -o ${prefix}_k${kmer_size}_s${sketch_size}.msh
    $sketchy check -g $genotype_file -r ${prefix}_k${kmer_size}_s${sketch_size}.msh
    """ 

}

process PredictBatch {

    tag { "reads=$reads sketch=$sketch_name" }
    label "predict"

    memory { params.predict_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/batch_predict", mode: "symlink", pattern: "*.txt"

    input:
    val(read_limit)
    file(genotype_file)
    file(read_files)
    each sketch

    output:
    file("${sketch_name}.txt")

    script:

    sketchy = params.exec ?: "sketchy"
    consensus = params.batch_consensus ?: ""
    sketch_name = sketch.baseName

    """
    for file in $read_files; do
        name=\$(basename "\$file" | cut -d. -f1)
        prediction=\$($sketchy predict -g $genotype_file -i \$file -l $read_limit -r $sketch $consensus)
        echo -e "\${name}\t\${prediction}" >> ${sketch_name}.txt
    done
    """ 

}