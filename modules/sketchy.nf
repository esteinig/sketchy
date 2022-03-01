process Sketch {

    tag { "prefix=$prefix k=$kmer_size s=$sketch_size" }
    label "sketch"

    memory { params.sketch_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/sketches", mode: "copy", pattern: "*.msh"

    input:
    val(prefix)
    val(fasta_glob) 
    file(genotype_file)
    file(fasta_dir)  // fasta directory symlink
    each kmer_size
    each sketch_size

    output:
    tuple val("${prefix}_k${kmer_size}_s${sketch_size}"), file("${prefix}_k${kmer_size}_s${sketch_size}.msh")

    script:

    sketchy = params.exec ?: "sketchy"

    """
    find . -name "${fasta_dir.baseName}/$fasta_glob" | $sketchy sketch -k $kmer_size -s $sketch_size -o ${prefix}_k${kmer_size}_s${sketch_size}.msh
    """ 

}

process PredictBatch {

    tag { "reads=$read_limit sketch=$sketch_name" }
    label "predict"

    memory { params.predict_mem * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/batch_predict", mode: "copy", pattern: "*.txt"

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
    first_file = read_files[0]

    """
    header=\$($sketchy predict -g $genotype_file -i $first_file -l 1 -r $sketch -t 0 -H) 
    echo -e "name\t\${header}" > ${sketch_name}.txt

    for file in $read_files; do
        name=\$(basename "\$file" | cut -d. -f1)
        prediction=\$($sketchy predict -g $genotype_file -i \$file -l $read_limit -r $sketch $consensus)
        echo -e "\${name}\t\${prediction}" >> ${sketch_name}.txt
    done
    """ 

}