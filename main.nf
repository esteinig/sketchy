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

params.fastq = "*.fq"
params.db = ""                               // must be actual path
params.reads = "10,50,100,200,500,1000"

params.limit = 1                             // output best prediction
params.ranks = 10
params.stability = 100

if (params.db) {
    dbs = params.db.split(",").collect { file(it) }
} else {
    println("You need to specify one or multiple databases with argument: --db")
    System.exit(1)
}

if (params.reads) {
    reads = params.reads.split(",").collect { file(it) }
} else {
    println("You need to specify one or multiple prediction end points with argument: --reads")
    System.exit(1)
}


include { SketchyStream } from './modules/sketchy'
include { SketchyScreen } from './modules/sketchy'
include { SketchyDist } from './modules/sketchy'

workflow {

    ont = channel.fromPath("${params.fastq}", type: 'file').map { tuple(it.simpleName, it) }

    SketchyStream(ont, reads, dbs)
    SketchyScreen(ont, reads, dbs)
    SketchyDist(ont, reads, dbs)

}
