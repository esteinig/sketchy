nextflow.enable.dsl=2

params.fasta = "ref_sketch/*.fa"
params.outdir = "sketch_builds"

include { Sketchy } from './modules/sketchy'


workflow build_sketch {

    fasta = channel.fromPath(params.fasta) | map { file -> tuple(file.baseName, file) }

}

workflow build_sketch {

    fasta = channel.fromPath(params.fasta) | map { file -> tuple(file.baseName, file) }

}

workflow cross_validation {
    // For each reference sketch, sample ten assemblies from the reference assemblies at random
    // then for each of those assemblies simulate nanopore reads, and build reference sketches
    // for each subsample with the samples removed from the reference sketch
    fasta = channel.fromPath(params.fasta) | map { file -> tuple(file.baseName, file) }


    // In the next step, predict selected traits across the simulated assembly samples, with and 
    // without them included in the reference sketch for evaluation of database representation

    
}
