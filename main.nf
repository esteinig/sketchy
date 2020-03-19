
log.info """

SKETCHY - NF v0.4.4
=======================

+++ UNDER DEVELOPMENT +++

PARALLEL SKETCHY COMPUTE
==========================

Run Sketchy on multiple samples and
sketch confirgurations in parallel

sketchy         : $params.sketchy
sketches        : $params.sketches
ranks           : $params.ranks
limits          : $params.limits

==========================

Additonally bootstrap each feature and 
compute 95% CI of prediction

bootstrap       : $params.bootstrap
samples         : $params.samples
limit           : $params.limit
seed            : $params.seed

METAGENOME : EXTRACT SPECIES
=============================

Extract a species by name from reads
classified with Kraken2

metagenome    : $params.metagenome
outdir        : $params.outdir
fastx         : $params.fastq
species       : $params.species
taxdb         : $params.taxdb
prefix        : $params.prefix

BUILD MASH SKETCH
=====================

Construct a Mash sketch from multiple 
Fasta files and merge

build         : $params.build
outdir        : $params.outdir
fasta         : $params.fasta
prefix        : $params.prefix
kmer_size     : $params.kmer_size
sketch_size   : $params.sketch_size

"""

if (params.sketchy){

    // DISTRIBUTED SKETCHY PREDICTION

    
    Channel
        .fromPath(params.fastq)
        .map { file -> tuple(file.baseName, file) }
        .into { fastq_nanopore; prestats_fastq }


    process StatsPrefilter {

        tag { id }
        label "ont"

        publishDir "$params.outdir/fastq", mode: "copy", pattern: "*.txt"

        input:
        set id, file(fq) from prestats_fastq

        output:
        file("${id}.prefiltered.stats.txt")

        shell:

        """
        if [[ $fq == *.gz ]]
        then
            zcat $fq | nanoq 2> ${id}.prefiltered.stats.txt
        else
            nanoq -f $fq 2> ${id}.prefiltered.stats.txt
        fi
        """

    }

    process Nanoq {
        
        tag { id }
        label "ont"

        publishDir "$params.outdir/fastq", mode: "copy", pattern: "*.txt"

        input:
        set id, file(fq) from fastq_nanopore

        output:
        set id, file("${id}.filtered.fq") into (poststats_fastq, rasusa_fastq, sketchy_fastq, sketchy_baseline, sketchy_bootstrap)

        """
        if [[ $fq == *.gz ]]
        then
            zcat $fq | nanoq -l $params.length -q $params.quality > ${id}.filtered.fq
        else
            nanoq -f $fq -l $params.length -q $params.quality > ${id}.filtered.fq
        fi
        """

    }

    process StatsPostfilter {

        tag { id }
        label "ont"

        publishDir "$params.outdir/fastq", mode: "copy", pattern: "*.txt"

        input:
        set id, file(filtered) from poststats_fastq

        output:
        file("${id}.filtered.stats.txt")

        """
        nanoq -f $filtered 2> ${id}.filtered.stats.txt
        """

    }


    process Sketchy {
        
        tag { "$id - $sketch - $rank - $limit" }
        label "sketchy"

        publishDir "$params.outdir/sketchy/${sketch}", mode: "copy"

        input:
        set id, file(fastq) from sketchy_fastq
        each sketch from params.sketches
        each rank from params.ranks
        each limit from params.limits
        
        output:
        file("${id}.${rank}.${limit}.*")
        
        when:
        params.sketchy

        script:

        """
        SKETCHY_PATH=$params.home

        sketchy run --fastq $fastq --sketch $sketch --ranks $rank --limit $limit --outdir sketchy --prefix ${id}.${rank}.${limit} --threads $task.cpus 
        mv sketchy/* .

        if [[ $params.time ]]
        then
            sketchy utils fx-time --fastq $fastq --evaluation ${id}.${rank}.${limit}.data.tsv --prefix ${id}.${rank}.${limit}
            rm *.fxi
        fi
     
        """

    }


}


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
        -p $task.cpus $fasta
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

    fastq = Channel
      .fromPath(params.fastq)
      .map { file -> tuple(file.baseName, file) }

    // Classify and extract species from classified reads

    process Kraken2 {

        label "kraken2"

        publishDir "${params.outdir}/metagenome", mode: "copy", pattern: "*.report"

        input:
        set id, file(fastq) from fastq

        output:
        set id, file("$fastq"), file("${id}.out"), val("$params.taxdb") into kraken_filter
        set id, file("${id}.report") into kraken_plot

        """
        kraken2 --db $params.taxdb --threads $task.cpus --output "${id}.out" --report ${id}.report --use-names $fastq
        """
    }

    process KrakenFilter {

        label "kraken2"

        publishDir "${params.outdir}/metagenome", mode: "copy"

        input:
        set id, file(fastq), file(kraken), val(taxdb) from kraken_filter

        output:
        set file("${id}.species.reads.out"), file("${id}.species.fq")

        """
        cat $kraken | grep "$params.species" > ${id}.species.reads.out
        sketchy utils fx-filter -f $fastq -i ${id}.species.reads.out -o ${id}.species.fq
        """
    }

}
