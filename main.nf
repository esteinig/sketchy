nextflow.enable.dsl=2

params.outdir = "sketchy_pipelines"
params.exec = null

// Sketch building

params.prefix = "test"
params.kmer_min = 16
params.kmer_max = 31
params.sketch_sizes = "1000,10000"
params.sketch_genomes = "*.fasta"
params.sketch_genomes_glob = "*.fasta"  // for large numbers of genomes, uses find pipe
params.sketch_genotypes = "genotypes.tsv"

def get_sketch_size_array(sketch_sizes){
    sketch_size_array = sketch_sizes.split(',');
    if (sketch_size_array){
        return sketch_size_array
    } else {
        return [1000]
    }
}

include { Sketch } from './modules/sketchy'

workflow sketch {
    sketch_sizes = get_sketch_size_array(params.sketch_sizes)
    Sketch(
        params.prefix, 
        params.sketch_genomes_glob, 
        file(params.sketch_genotypes), 
        channel.fromPath(params.sketch_genomes) | collect, 
        params.kmer_min..params.kmer_max, 
        sketch_sizes
    )
}
