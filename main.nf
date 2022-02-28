nextflow.enable.dsl=2

params.outdir = "sketchy_pipelines"

// Sketch building

params.prefix = "test"
params.kmer_min = 16
params.kmer_max = 21
params.sketch_sizes = [1000, 10000]
params.sketch_genomes = "*.fasta"
params.sketch_genotypes = "genotypes.tsv"

include { Sketch } from './modules/sketchy'



workflow sketch {
    fasta_files = channel.fromPath(params.sketch_genomes).collect()
    sketch_inputs = tuple(params.prefix, params.sketch_genomes, params.sketch_genotypes, fasta_files)
    Sketch(sketch_inputs, params.kmer_min..params.kmer_max, params.sketch_sizes)
}
