extern crate cute;
extern crate clap;
extern crate indicatif;
extern crate prettytable;

mod sketchy;

use std::io::Error;
use clap::{Arg, App, SubCommand};

fn main() -> Result<(), Error> {

    let matches = App::new("sketchy")
        .version("0.5.0")
        .about("\nNanopore lineage calling and genotyping of bacterial pathogens using Mash\n")
        .subcommand(SubCommand::with_name("stream")
            .about("\nsum of shared hashes from fasta/q on stdin")
            .version("0.5.0")
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("reference sketch db"))
            .arg(Arg::with_name("RANKS").short("r").long("ranks").takes_value(true).help("max ssh ranks per read"))
            .arg(Arg::with_name("STABILITY").short("s").long("stability").takes_value(true).help("reads to stable breakpointh"))
            .arg(Arg::with_name("THREADS").short("t").long("threads").takes_value(true).help("max threads for mash"))
            .arg(Arg::with_name("PROGRESS").short("p").long("progress").takes_value(false).help("progress bar on"))
        )
        .subcommand(SubCommand::with_name("screen")
            .about("\nscreen read set against reference sketch")
            .version("0.5.0")
            .arg(Arg::with_name("FASTX").short("f").long("fastx").takes_value(true).required(true).help("fasta/q input path"))
            .arg(Arg::with_name("DB").short("d").long("db").takes_value(true).required(true).help("reference sketch db"))
            .arg(Arg::with_name("LIMIT").short("l").long("limit").takes_value(true).help("limit ranked results"))
            .arg(Arg::with_name("THREADS").short("t").long("threads").takes_value(true).help("max threads for mash"))
            .arg(Arg::with_name("PRETTY").short("p").long("pretty").takes_value(false).help("pretty print on"))
        )
        .get_matches();
        
    if let Some(stream) = matches.subcommand_matches("stream") {
        
        let db: String = stream.value_of("db").unwrap_err("TEST").to_string();
        let ranks: usize = stream.value_of("ranks").unwrap_or("10").parse::<usize>().unwrap();
        let threads: i32 = stream.value_of("threads").unwrap_or("4").parse::<i32>().unwrap();
        let stability: usize = stream.value_of("stability").unwrap_or("100").parse::<usize>().unwrap();
        let progress: bool = stream.is_present("progress");

        let (sketch_msh, _, genotype_index, _) = sketchy::get_sketch_files(db);
        let (sketch_size, sketch_index): (usize, usize) = sketchy::get_sketch_info(&sketch_msh);
        

        sketchy::run(sketch_msh, genotype_index, threads, ranks, stability, progress, sketch_index, sketch_size).map_err(
            |err| println!("{:?}", err)
        ).ok();
        
    }

    if let Some(screen) = matches.subcommand_matches("screen") {
        
        let fastx: String = screen.value_of("fastx").unwrap().unwrap_err("TEST").to_string();
        let db: String = screen.value_of("db").unwrap().unwrap_err("TEST").to_string();
        let threads: i32 = screen.value_of("threads").unwrap_or("4").parse::<i32>().unwrap();
        let limit: usize = screen.value_of("limit").unwrap_or("10").parse::<usize>().unwrap();
        let pretty: bool = screen.is_present("pretty");

        let (sketch_msh, genotypes, _, _) = sketchy::get_sketch_files(db);

        sketchy::screen(fastx, sketch_msh, genotypes, threads, limit).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    Ok(())

}
