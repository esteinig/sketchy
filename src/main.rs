extern crate dirs;
extern crate cute;
extern crate clap;
extern crate colored; 
extern crate indicatif;
extern crate serde_json;
extern crate prettytable;

mod sketchy;

use std::io::Error;
use clap::{Arg, App, AppSettings, SubCommand};

fn main() -> Result<(), Error> {


    let matches = App::new("sketchy")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::DisableHelpSubcommand)
        .version("0.5.0")
        .about("\nGenomic neighbor typing using MinHash\n")
        .subcommand(SubCommand::with_name("stream")
            .about("\ncompute sum of shared hashes from fasta/q stream")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Genomic neighbor typing DB [required]"))
            .arg(Arg::with_name("FASTX").short("f").long("fastx").takes_value(true).help("Fasta/q path or STDIN [-]"))
            .arg(Arg::with_name("READS").short("r").long("reads").takes_value(true).help("Limit reads to prediction [none]"))
            .arg(Arg::with_name("RANKS").short("n").long("ranks").takes_value(true).help("Ranked sum of shared hashes [10]"))
            .arg(Arg::with_name("STABILITY").short("s").long("stability").takes_value(true).help("Reads to stable breakpoint [100]"))
            .arg(Arg::with_name("THREADS").short("t").long("threads").takes_value(true).help("Maximum threads for Mash [4]"))
            .arg(Arg::with_name("PROGRESS").short("p").long("progress").takes_value(false).help("Progress bar on [false]"))
            .arg(Arg::with_name("RAW").short("w").long("raw").takes_value(false).help("Print raw sum of shared hashes [false]"))
        )
        .subcommand(SubCommand::with_name("predict")
            .about("\npredict genotypes from sum of shared hashes stream")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Genomic neighbor typing DB [required]"))
            .arg(Arg::with_name("LIMIT").short("l").long("limit").takes_value(true).help("Limit predicted rank genotypes per read [10]"))
            .arg(Arg::with_name("RAW").short("w").long("raw").takes_value(false).help("Print raw translated genotypes [false]"))
            .arg(Arg::with_name("PRETTY").short("p").long("pretty").takes_value(false).help("Pretty print on [false]"))
        )
        .subcommand(SubCommand::with_name("screen")
            .about("\nquery read set against database with mash screen")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Genomic neighbor typing DB [required]"))
            .arg(Arg::with_name("FASTX").short("f").long("fastx").takes_value(true).required(true).help("Fasta/q input path [required]"))
            .arg(Arg::with_name("LIMIT").short("l").long("limit").takes_value(true).help("Limit predicted genotype output [10]"))
            .arg(Arg::with_name("THREADS").short("t").long("threads").takes_value(true).help("Maximum threads for Mash [4]"))
            .arg(Arg::with_name("PRETTY").short("p").long("pretty").takes_value(false).help("Pretty print on [false]"))
        )
        .subcommand(SubCommand::with_name("dist")
            .about("\nquery read set against database with mash dist")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Genomic neighbor typing DB [required]"))
            .arg(Arg::with_name("FASTX").short("f").long("fastx").takes_value(true).required(true).help("Fasta/q input path [required]"))
            .arg(Arg::with_name("LIMIT").short("l").long("limit").takes_value(true).help("Limit predicted genotype output [10]"))
            .arg(Arg::with_name("THREADS").short("t").long("threads").takes_value(true).help("Maximum threads for Mash [4]"))
            .arg(Arg::with_name("PRETTY").short("p").long("pretty").takes_value(false).help("Pretty print on [false]"))
        )
        .subcommand(SubCommand::with_name("head")
            .about("\ndisplay genotype database header")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Genomic neighbor typing DB [required]"))
            .arg(Arg::with_name("PRETTY").short("p").long("pretty").takes_value(false).help("Pretty print on [false]"))
        )
        .subcommand(SubCommand::with_name("check")
            .about("\ncheck genotype database format")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Genomic neighbor typing DB [required]"))
        )
        .subcommand(SubCommand::with_name("cite")
            .about("\noutput citations for sketchy")
            .version("0.5.0")
        )
        .get_matches();
        
    if let Some(stream) = matches.subcommand_matches("stream") {
        
        let db: String = stream.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let fastx: String = stream.value_of("FASTX").unwrap_or("-").to_string();
        let reads: u32 = stream.value_of("READS").unwrap_or("0").parse::<u32>().unwrap();
        let ranks: usize = stream.value_of("RANKS").unwrap_or("10").parse::<usize>().unwrap();
        let threads: i32 = stream.value_of("THREADS").unwrap_or("4").parse::<i32>().unwrap();
        let stability: usize = stream.value_of("STABILITY").unwrap_or("100").parse::<usize>().unwrap();
        let progress: bool = stream.is_present("PROGRESS");
        let raw: bool = stream.is_present("RAW");

        let (sketch_msh, _, genotype_index, _) = sketchy::get_sketch_files(db);
        let (sketch_size, sketch_index): (usize, usize) = sketchy::get_sketch_info(&sketch_msh);

        sketchy::stream(fastx, sketch_msh, genotype_index, threads, reads, ranks, stability, progress, raw, sketch_index, sketch_size).map_err(
            |err| println!("{:?}", err)
        ).ok();
        
    }

    if let Some(screen) = matches.subcommand_matches("screen") {
        
        let db: String = screen.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();
        
        let fastx: String = screen.value_of("FASTX").unwrap_or_else(||
            clap::Error::with_description("Could not find input read file", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let threads: i32 = screen.value_of("THREADS").unwrap_or("4").parse::<i32>().unwrap();
        let limit: usize = screen.value_of("LIMIT").unwrap_or("10").parse::<usize>().unwrap();
        let pretty: bool = screen.is_present("PRETTY");

        let (sketch_msh, genotypes, _, _) = sketchy::get_sketch_files(db);

        sketchy::screen(fastx, sketch_msh, genotypes, threads, limit, pretty).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }
    
    if let Some(dist) = matches.subcommand_matches("dist") {
        
        let db: String = dist.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();
        
        let fastx: String = dist.value_of("FASTX").unwrap_or_else(||
            clap::Error::with_description("Could not find input read file", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let threads: i32 = dist.value_of("THREADS").unwrap_or("4").parse::<i32>().unwrap();
        let limit: usize = dist.value_of("LIMIT").unwrap_or("10").parse::<usize>().unwrap();
        let pretty: bool = dist.is_present("PRETTY");

        let (sketch_msh, genotypes, _, _) = sketchy::get_sketch_files(db);

        sketchy::dist(fastx, sketch_msh, genotypes, threads, limit, pretty).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    if let Some(predict) = matches.subcommand_matches("predict") {

        let db: String = predict.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let limit: usize = predict.value_of("LIMIT").unwrap_or("1").parse::<usize>().unwrap();
        let raw: bool = predict.is_present("RAW");
        let pretty: bool = predict.is_present("PRETTY");

        let (_, _, _, genotype_key) = sketchy::get_sketch_files(db);

        sketchy::predict(genotype_key, limit, raw, pretty).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    if let Some(head) = matches.subcommand_matches("head") {

        let db: String = head.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let pretty: bool = head.is_present("PRETTY");

        let (_, _, _, genotype_key) = sketchy::get_sketch_files(db);

        sketchy::display_header(genotype_key, pretty).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    if let Some(check) = matches.subcommand_matches("check") {

        let db: String = check.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        sketchy::get_sketch_files(db);
        println!("Ok");
    }

    if let Some(_cite) = matches.subcommand_matches("cite") {
                
        println!("\nPlease cite the following authors whose work we used for Sketchy:\n");
        
        println!("Ondov  et al. (2016) : https://doi.org/10.1186/s13059-016-0997-x");
        println!("Ondov  et al. (2019) : https://doi.org/10.1186/s13059-019-1841-x");
        println!("Brinda et al. (2020) : https://doi.org/10.1038/s41564-019-0656-6\n");
        
    }

    Ok(())

}
