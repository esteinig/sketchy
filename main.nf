
log.info """

===============================================
             SKETCHY - NF v0.4.4
===============================================

outdir        : $params.outdir

===============================================
           PARALLEL SKETCHY COMPUTE          
===============================================

Run Sketchy on multiple samples and sketch 
configurations in parallel

sketchy         : $params.sketchy
fastq reads     : $params.fastq

min quality     : $params.quality
min length      : $params.length
sort by time    : $params.sort

sketches        : $params.sketches
ranks           : $params.ranks
limits          : $params.limits
time            : $params.time
delta           : $params.delta

Additonally bootstrap each feature and compute 
95% bootrap interval for each prediction

! Should be used on clusters only !

bootstrap       : $params.bootstrap
samples         : $params.samples
seed            : $params.seed

===============================================
         METAGENOME : EXTRACT SPECIES          
===============================================

Extract a species by scientific name from reads 
classified with Kraken2

metagenome    : $params.metagenome
fastq         : $params.fastq
species       : $params.species
taxdb         : $params.taxdb

===============================================
             BUILD MASH SKETCH            
===============================================

Construct Mash sketch in parallel

build         : $params.build
fasta         : $params.fasta
prefix        : $params.prefix
kmer_size     : $params.kmer_size
sketch_size   : $params.sketch_size

===============================================
===============================================
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
        set id, file("${id}.filtered.fq") into (poststats_fastq, rasusa_fastq, coverm_fastq, sketchy_fastq, sketchy_baseline, sketchy_bootstrap)


        """
        if [[ $params.length > 0 && $params.quality > 0 ]]
        then
            if [[ $fq == *.gz ]]
            then
                zcat $fq | nanoq -l $params.length -q $params.quality > ${id}.filtered.fq
            else
                nanoq -f $fq -l $params.length -q $params.quality > ${id}.filtered.fq
            fi
        else
            ln -s $fq ${id}.filtered.fq
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
        if [[ $params.length > 0 && $params.quality > 0 ]]
        then
            nanoq -f $filtered 2> ${id}.filtered.stats.txt
        else
            touch ${id}.filtered.stats.txt
        fi
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
        SKETCHY_PATH=$params.sketchy_path
        
        if [[ $params.sort = true ]]
        then
            sketchy utils fx-sort --fastx $fastq --output sketchy.fq
        else
            ln -s $fastq sketchy.fq
        fi

        sketchy run --fastq sketchy.fq --sketch $sketch --ranks $rank --limit $limit --stable $params.stable --outdir sketchy --prefix ${id}.${rank}.${limit} --threads $task.cpus 
        mv sketchy/* .

        if [[ $params.time = true ]]
        then
            sketchy utils fx-time --fastq sketchy.fq --evaluation ${id}.${rank}.${limit}.data.tsv --prefix ${id}.${rank}.${limit}
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

if (params.bbuild) {

    // Build a bootstrapped sketch from Pathfinder data by randomly sampling
    // <sample> genomes from the input files <bootstrap> times with replacement

    process BootstrapSketch {

        tag { "$sample - $replicate" }

        label "bbuild"

        publishDir "${params.outdir}/bbuild/sketches/${sample}"

        input:
        file(iid) from Channel.fromPath(params.iid)
        file(survey) from Channel.fromPath(params.survey)
        file(features) from Channel.fromPath(params.features)
        each sample from params.samples
        each replicate from params.bootstraps

        output:
        file(replicate)
        set sample, replicate, file("${replicate}.msh"), file("${replicate}.json"), file("${replicate}.tsv") into sketchy_bootstrap

        """
        sketchy survey link -i $iid -d $survey -b $sample -s -o $replicate
        mash sketch -k 15 -s 1000 $replicate/*.fasta -o $replicate
        sketchy feature merge -s ${replicate}.msh  -f $features -i key -p merged
        sketchy feature prepare -i merged.tsv -d key -p ${replicate}
        """

    }

    process BootstrapSketchy {

        tag { "$sample - $replicate" }
        label "sketchy"

        publishDir "$params.outdir/bbuild/sketchy/${sample}/${fastq}", mode: "copy"

        input:
        set sample, replicate, file(sketch), file(key), file(index) from sketchy_bootstrap
        each file(fastq) from Channel.fromPath(params.fastq).collect()

        output:
        set file("${replicate}.data.tsv"), file("${replicate}.ssh.tsv"), file("${replicate}.sssh.tsv")

        """
        sketchy run --fastq $fastq --sketch $replicate --ranks 10 --limit 10 --stable 10 --outdir sketchy --prefix $replicate --threads $task.cpus
        mv sketchy/* .

        """

    }



}
