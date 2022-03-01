nextflow.enable.dsl=2

// General 

params.help = false
params.outdir = "sketchy_pipelines"
params.exec = null

// Sketch construction

params.prefix = "test"
params.kmer_min = 16
params.kmer_max = 31
params.sketch_sizes = "1000,10000"
params.sketch_genomes_dir = "test/"
params.sketch_genomes_glob = "*.fasta"  // for large numbers of genomes, uses find pipe

// Batch prediction

params.batch_consensus = null
params.batch_read_limit = 1000
params.batch_sketch_files = "*.msh"
params.batch_read_files = "*.fastq"
params.batch_genotype_file = "genotypes.tsv"

// Module imports

include { helpMessage } from './modules/help'
include { Sketch } from './modules/sketchy'
include { PredictBatch } from './modules/sketchy'

// Help message

if (params.help) {
    helpMessage()
}

sketch_sizes = param.sketch_sizes.tokenize(",")

System.exit(0)

// Workflow entry points

workflow sketch {
    Sketch(
        params.prefix, 
        params.sketch_genomes_glob, 
        file(params.sketch_genomes_dir), 
        params.kmer_min..params.kmer_max, 
        sketch_sizes
    )
}

workflow batch_predict {
    PredictBatch(
        params.batch_read_limit, 
        file(params.batch_genotype_file), 
        channel.fromPath(params.batch_read_files) | collect,
        channel.fromPath(params.batch_sketch_files) | collect
    )
}
