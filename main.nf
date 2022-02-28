nextflow.enable.dsl=2

params.outdir = "sketchy_pipelines"
params.exec = null

// Sketch building

params.prefix = "test"
params.kmer_min = 16
params.kmer_max = 21
params.sketch_sizes = [1000, 10000]
params.sketch_genomes = "*.fasta"
params.sketch_genomes_glob = "*.fasta"  // for large numbers of genomes, uses ls pipe
params.sketch_genotypes = "genotypes.tsv"

include { Sketch } from './modules/sketchy'

workflow sketch {
    fasta_files = channel.fromPath(params.sketch_genomes) | collect | map { file(it) }
    Sketch(tuple(params.prefix, params.sketch_genomes_glob), file(params.sketch_genotypes), fasta_files, params.kmer_min..params.kmer_max, params.sketch_sizes)
}
