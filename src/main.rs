extern crate dirs;
extern crate cute;
extern crate clap;
extern crate indicatif;
extern crate serde_json;
extern crate prettytable;

mod sketchy;

use std::env;
use std::io::Error;
use clap::{Arg, App, SubCommand};

fn main() -> Result<(), Error> {

    let user_home: String = dirs::home_dir().unwrap().to_str().unwrap_or("").to_string();
    let sketchy_home: String = format!("{}/.sketchy", user_home);
    let sketchy_path: String = env::var("SKETCHY_PATH").unwrap_or(sketchy_home).to_string();

    let matches = App::new("sketchy")
        .version("0.5.0")
        .about("\nNanopore lineage calling and genotyping of bacterial pathogens using Mash\n")
        .subcommand(SubCommand::with_name("stream")
            .about("\ncompute sum of ranked sum of shared hashes from fasta/q stream")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Reference sketch DB [required]"))
            .arg(Arg::with_name("FASTX").short("f").long("fastx").takes_value(true).help("Fasta/q path or STDIN [-]"))
            .arg(Arg::with_name("RANKS").short("r").long("ranks").takes_value(true).help("Ranked sum of shared hashes [10]"))
            .arg(Arg::with_name("STABILITY").short("s").long("stability").takes_value(true).help("Reads to stable breakpoint [100]"))
            .arg(Arg::with_name("THREADS").short("t").long("threads").takes_value(true).help("Maximum threads for Mash [4]"))
            .arg(Arg::with_name("PROGRESS").short("p").long("progress").takes_value(false).help("Progress bar on [false]"))
            .arg(Arg::with_name("RAW").short("w").long("raw").takes_value(false).help("Print raw sum of shared hashes [false]"))
        )
        .subcommand(SubCommand::with_name("predict")
            .about("\npredict genotypes from sum of shared hashes stream")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Reference sketch DB [required]"))
            .arg(Arg::with_name("LIMIT").short("l").long("limit").takes_value(true).help("Limit predicted rank genotypes per read [10]"))
            .arg(Arg::with_name("RAW").short("r").long("raw").takes_value(false).help("Print raw translated genotypes [false]"))
        )
        .subcommand(SubCommand::with_name("screen")
            .about("\nscreen read set against reference sketch with mash")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Reference sketch DB [required]"))
            .arg(Arg::with_name("FASTX").short("f").long("fastx").takes_value(true).required(true).help("Fasta/q input path [required]"))
            .arg(Arg::with_name("LIMIT").short("l").long("limit").takes_value(true).help("Limit predicted genotype output [10]"))
            .arg(Arg::with_name("THREADS").short("t").long("threads").takes_value(true).help("Maximum threads for Mash [4]"))
            .arg(Arg::with_name("PRETTY").short("p").long("pretty").takes_value(false).help("Pretty print on [false]"))
        )
        .subcommand(SubCommand::with_name("head")
            .about("\ndisplay genotype database header")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Reference sketch DB [required]"))
            .arg(Arg::with_name("PRETTY").short("p").long("pretty").takes_value(false).help("Pretty print on [false]"))
        )
        .subcommand(SubCommand::with_name("check")
            .about("\ncheck genotype database format")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("Reference sketch DB [required]"))
        )
        .get_matches();
        
    if let Some(stream) = matches.subcommand_matches("stream") {
        
        let db: String = stream.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let fastx: String = stream.value_of("FASTX").unwrap_or("-").to_string();
        let ranks: usize = stream.value_of("RANKS").unwrap_or("10").parse::<usize>().unwrap();
        let threads: i32 = stream.value_of("THREADS").unwrap_or("4").parse::<i32>().unwrap();
        let stability: usize = stream.value_of("STABILITY").unwrap_or("100").parse::<usize>().unwrap();
        let progress: bool = stream.is_present("PROGRESS");
        let raw: bool = stream.is_present("RAW");

        let (sketch_msh, _, genotype_index, _) = sketchy::get_sketch_files(db, &sketchy_path);
        let (sketch_size, sketch_index): (usize, usize) = sketchy::get_sketch_info(&sketch_msh);

        sketchy::stream(fastx, sketch_msh, genotype_index, threads, ranks, stability, progress, raw, sketch_index, sketch_size).map_err(
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

        let (sketch_msh, genotypes, _, genotype_key) = sketchy::get_sketch_files(db, &sketchy_path);

        sketchy::screen(fastx, sketch_msh, genotypes, genotype_key, threads, limit, pretty).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    if let Some(predict) = matches.subcommand_matches("predict") {

        let db: String = predict.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let limit: usize = predict.value_of("LIMIT").unwrap_or("1").parse::<usize>().unwrap();
        let raw: bool = predict.is_present("RAW");

        let (_, _, _, genotype_key) = sketchy::get_sketch_files(db, &sketchy_path);

        sketchy::predict(genotype_key, limit, raw).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    if let Some(head) = matches.subcommand_matches("head") {

        let db: String = head.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        let pretty: bool = head.is_present("PRETTY");

        let (_, _, _, genotype_key) = sketchy::get_sketch_files(db, &sketchy_path);

        sketchy::display_header(genotype_key, pretty).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    if let Some(check) = matches.subcommand_matches("check") {

        let db: String = check.value_of("DB").unwrap_or_else(||
            clap::Error::with_description("Please input a reference sketch database", clap::ErrorKind::InvalidValue).exit()
        ).to_string();

        sketchy::get_sketch_files(db, &sketchy_path);

    }

    Ok(())

}
