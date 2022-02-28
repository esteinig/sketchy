nextflow.enable.dsl=2

params.outdir = "sketchy_pipelines"

// Sketch building

params.prefix = "test"
params.kmer_range = 16..31
params.sketch_sizes = [1000, 10000]
params.sketch_genomes = "*.fasta"
params.sketch_genotypes = "genotypes.tsv"

// GNT predictions

params.test_reads = "*.fq"

include { Sketch } from './modules/sketchy'

workflow sketch {
    fasta_files = channel.fromPath(params.sketch_genomes).collect()
    sketch_inputs = tuple(params.prefix, params.sketch_genomes, fasta_files)
    Sketch(sketch_inputs, params.kmer_range, params.sketch_sizes)
}
