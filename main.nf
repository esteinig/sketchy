#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
==============================================================================================================================
                                        S K E T C H Y   P I P E L I N E
==============================================================================================================================

 Genomic neighbor typing using MinHash

 Documentation: https://github.com/esteinig/sketchy

 Original development by Australian Institute of Tropical Health and Medicince, 
 Queensland Genomics, The Peter Doherty Intitute for Infection and Immunity

Developers:

    Eike Steinig  @esteinig  < @EikeSteinig >
    Michael Hall  @mbhall88  < @mbhall88    >

----------------------------------------------------------------------------------------
*/

import java.nio.file.Paths
import groovy.io.*

nextflow.enable.dsl=2




// Helper functions



version = '0.5.0'

def helpMessage() {

    log.info"""
    =========================================
          S K E T C H Y   v${version}
    =========================================

    Usage:

    Workflows:


    Help:

        nextflow run np-core/np-variants --workflow sketchy --help true

    =========================================

    """.stripIndent()
}


params.outdir = "nxf-sketchy"
params.workflow = "genotype"
params.fastq = "*.fastq"
params.db = ""                               // must be actual path
params.reads = "20,50,100,200,500,1000,2000,5000,10000"

params.limit = 1                             // output best prediction
params.ranks = 10
params.stability = 100
params.replicates = 0
params.fasta_directory = ""
params.reference_database = ""
params.create_options = ""

reps = 1..params.replicates
samples = [200, 500, 2000, 8000, 16000, 32000]

bootstrap_read_limit = 1000

if (params.db) {
    dbs = params.db.split(",").collect { file(it) }
}

if (params.reads) {
    read_limits = params.reads.split(",").collect { it.toInteger() }
}

if (params.fasta_directory){
    fasta_directory = file(params.fasta_directory)
}

if (params.reference_database){
    reference_database = file(params.reference_database)
}




def startMessage() {

    log.info"""
    =========================================
          S K E T C H Y   v${version}
    =========================================

    Fastq:                  ${params.fastq}
    Outdir:                 ${params.outdir}
    Databases:              ${dbs}
    Read limits:            ${read_limits}
    Bootstraps replicates:  ${reps}

    Prediction limit:       ${params.limit} [1 = best prediction]
    Prediction ranks:       ${params.ranks}
    Prediction stability:   ${params.stability}
    
    Fasta directory:        ${fasta_directory}
    Bootstrap samples:      ${samples}
    Bootstrap read limit:   ${bootstrap_read_limit}

    =========================================

    """.stripIndent()
}

include { SketchyStream } from './modules/sketchy'
include { SketchyScreen } from './modules/sketchy'
include { SketchyScreenWinner } from './modules/sketchy'
include { SketchyDist } from './modules/sketchy'
include { Bootstrap } from './modules/sketchy'
include { ReadLimit } from './modules/sketchy'

include { PublishHeader } from './modules/bootstrap'
include { BootstrapBuild } from './modules/bootstrap'
include { SketchyStream as BootstrapStream } from './modules/bootstrap'

workflow {

    if (params.workflow == "genotype"){
        ont = channel.fromPath("${params.fastq}", type: 'file').map { tuple(it.simpleName, it) }
        
        if (params.replicates > 0) {
            fq = Bootstrap(ont, reps, read_limits)
        } else {
            fq = ReadLimit(ont, read_limits)
        }

        SketchyStream(fq, dbs)
        SketchyScreen(fq, dbs)
        SketchyDist(fq, dbs)
        SketchyScreenWinner(fq, dbs)
    } 

    if (params.workflow == 'bootstrap_database') {
        ont = channel.fromPath("${params.fastq}", type: 'file').map { tuple(it.simpleName, bootstrap_read_limit, it) }
        
        bootstrap_dbs = BootstrapBuild(fasta_directory, reference_database, samples, reps)
        bootstrap_replicates = ont.combine(bootstrap_dbs)
        BootstrapStream(bootstrap_replicates)

        PublishHeader(reference_database, samples)


    }
   

}
