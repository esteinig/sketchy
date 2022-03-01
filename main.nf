nextflow.enable.dsl=2

params.outdir = "sketchy_pipelines"
params.exec = null

// Sketch building

params.prefix = "test"
params.kmer_min = 16
params.kmer_max = 31
params.sketch_sizes = [1000]
params.sketch_genomes_dir = "test/"
params.sketch_genomes_glob = "*.fasta"  // for large numbers of genomes, uses find pipe

// Batch predicing

params.batch_consensus = "-t 5 -c" // null to disable
params.batch_read_limit = 1000
params.batch_sketch_files = "*.msh"
params.batch_genotypes = "genotypes.tsv"
params.batch_read_files = "*.fastq"

include { Sketch } from './modules/sketchy'
include { PredictBatch } from './modules/sketchy'

workflow sketch {
    Sketch(
        params.prefix, 
        params.sketch_genomes_glob, 
        file(params.sketch_genomes_dir), 
        params.kmer_min..params.kmer_max, 
        params.sketch_sizes
    )
}

workflow batch_predict {
    PredictBatch(
        params.batch_read_limit, 
        file(params.batch_genotypes), 
        channel.fromPath(params.batch_read_files) | collect,
        channel.fromPath(params.batch_sketch_files) | collect
    )
}
