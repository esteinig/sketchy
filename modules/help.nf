
def helpMessage(){
    log.info """
===================================================================================================
                                    S K E T C H Y - N F X                   
===================================================================================================                             
                                                                                    
                                                                                            
                    ██████        ██████████  ██████████        ██████                  
                    ████        ██████████████████████████████        ████                
                    ████    ██████████████████████████████████████    ████                                     
                    ██████████████████████████████████████████████████████                
                    ████████████████████████  ████████████████████████                  
                        ██████████████████          ██████████████████                    
                            ████████                      ████████                        
                                                                                                        
                                                                            

===================================================================================================
                                    General Command Line Arguments                   
===================================================================================================

--outdir    output directory for workflow results             [$params.outdir]
--exec      sketchy executable path, default env "sketchy"    [$params.exec]
--help      show this help message and exit                   [$params.help]


===================================================================================================
                                    Reference sketch construction                  
===================================================================================================

-entry sketch

Construct reference sketches from arbitrarily large genome
assembly collections across a range of k-mer sizes (k).

Command line options:

--prefix                 prefix for output files                    [$params.prefix]
--kmer_min               minimum k-mer size for range               [$params.kmer_min]
--kmer_max               maximum k-mer size for range               [$params.kmer_min]
--sketch_size            sketch size to construct sketch            [$params.sketch_size]
--sketch_genomes_dir     path to directory containing assemblies    [$params.sketch_genomes_dir]
--sketch_genomes_glob    glob to grab assemblies from directory     [$params.sketch_genomes_dir]

Example:

    nextflow run esteinig/sketchy -profile local -entry sketch --prefix test \\
        --kmer_min 16 \\
        --kmer_max 31 \\
        --sketch_sizes 1000,10000 \\
        --sketch_genomes_dir reference/assemblies \\
        --sketch-genomes_glob "*.fasta"

Notes:

    * glob with wildcard must be string quoted e.g. "*.fasta"

Resources:

--sketch_mem     memory for sketching process          [$params.sketch_mem]
--sketch_cpus    threads for sketching process         [$params.sketch_cpus]
--sketch_time    maximum time for sketching process    [$params.sketch_time]

===================================================================================================
                                        Batch prediction                   
===================================================================================================

-entry batch_predict

Conduct predictions over multiple input files and sketches; adds
a name column into the output table for each input file base name.

Command line options:

--batch_consensus        use consensus calling options      [$params.batch_consensus]
--batch_read_limit       read limit for predictions         [$params.batch_read_limit]
--batch_sketch_files     reference sketch file glob         [$params.batch_sketch_files]
--batch_read_files       reads file glob for prediction     [$params.batch_read_files]
--batch_genotype_file    genotype file matching sketches    [$params.batch_genotype_file]

Example:

    nextflow run esteinig/sketchy -profile local -entry batch_predict --prefix test \\
        --batch_consensus "-t 5 -c" \\
        --batch_read_limit 1000 \\
        --batch_sketch_files "*.msh" \\
        --batch_read_files "*.fastq" \\
        --batch_genotype_file "genotypes.tsv"

Notes:

    * batch consensus must be string quoted and specify the --top and --consensus parameters
    * globs with wildcard (sketch files, read files) must be string quoted e.g. "*.msh"


Resources:

--sketch_mem     memory for prediction process          [$params.batch_predict_mem]
--sketch_cpus    threads for prediction process         [$params.batch_predict_cpus]
--sketch_time    maximum time for prediction process    [$params.batch_predict_time]

"""
}
