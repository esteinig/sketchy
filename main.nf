
log.info """

SKETCHY - NF v0.4.1
=======================

+++ UNDER DEVELOPMENT +++

PARALLEL SKETCHY COMPUTE
==========================

Run Sketchy on multiple samples in 
parallel.

sketchy       : $params.sketchy

METAGENOME : EXTRACT SPECIES
=============================

Extract a species by name from reads
classified with Kraken2

metagenome    : $params.metagenome
outdir        : $params.outdir
fastx         : $params.fastx
species       : $params.species
taxdb         : $params.taxdb
prefix        : $params.prefix
threads       : $params.threads

BUILD MASH SKETCH
=====================

Construct a Mash sketch from multiple 
Fasta files and merge.

build         : $params.build
outdir        : $params.outdir
fasta         : $params.fasta
prefix        : $params.prefix
kmer_size     : $params.kmer_size
sketch_size   : $params.sketch_size
threads       : $params.threads

"""



if (params.build) {

    // DISTRIBUTED MASH SKETCH

    genomes = Channel.fromPath(params.fasta, followLinks: true, glob: true)
                     .map { file -> tuple(file.baseName, file) }

    process GenomeSketch {

        // - sketch a single genome with Mash for merge in subsequent process

        label "build"
        
        input:
        set id, file(fasta) from genomes

        output:
        file("${id}.msh") into mash_files

        """
        echo $fasta $id
        mash sketch -k $params.kmer_size -s $params.sketch_size -o $id \
        -p $params.threads $fasta
        """
    }

    mash_list = mash_files.collect()

    process SketchAssembly {

        // - combine all sketches into a single sketch file

        label "build"
        publishDir "${params.outdir}", mode: "copy"

        input:
        file(sketches) from mash_list

        output:
        file("${params.prefix}.msh")

        """
        echo $sketches
        mash paste $params.prefix $sketches
        """

    }

}

if (params.metagenome) {

  fastx = Channel
      .fromPath(params.fastx)
      .map { file -> tuple(file.baseName, file) }

  // Classify and extract species reads from Zymo Mocks

  process Kraken2 {

        label "kraken2"

        publishDir "${params.outdir}/metagenome", mode: "copy", pattern: "*.report"

        input:
        set id, file(fastx) from fastx

        output:
        set id, file("$fastx"), file("${id}.out"), val("$params.taxdb") into kraken_filter
        set id, file("${id}.report") into kraken_plot

        """
        kraken2 --db $params.taxdb --threads $task.cpus --output "${id}.out" --report ${id}.report --use-names $fastx
        """
    }

    process KrakenFilter {

        label "kraken2"

        publishDir "${params.outdir}/metagenome", mode: "copy"

        input:
        set id, file(fastx), file(kraken), val(taxdb) from kraken_filter

        output:
        set file("${id}.species.reads.out"), file("${id}.species.fq")

        """
        cat $kraken | grep "$params.species" > ${id}.species.reads.out
        sketchy utils fx-filter -f $fastx -i ${id}.species.reads.out -o ${id}.species.fq
        """
    }

    process KrakenPlot {

        label "sketchy"

        publishDir "${params.outdir}/metagenome", mode: "copy"

        input:
        set id, file(report) from kraken_plot

        output:
        file("${id}.png")

        """
        sketchy utils plot-kraken --report $report --prefix $id
        """
    }


}
