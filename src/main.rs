extern crate cute;
extern crate clap;
extern crate indicatif;
extern crate prettytable;

mod sketchy;

use std::io::Error;
use clap::{Arg, App, SubCommand};

fn main() -> Result<(), Error> {

    let matches = App::new("sketchy")
        .version("0.4.5")
        .about("\nNanopore lineage calling and genotyping of bacterial pathogens using Mash\n")
        .subcommand(SubCommand::with_name("stream")
            .about("\ncompute sum of shared hashes from fasta/q on stdin")
            .version("0.4.5")
            .arg(Arg::from_usage("-s, --sketch=[FILE] 'reference sketch db'"))
            .arg(Arg::from_usage("-r, --ranks=[INT] 'max ssh ranks per read'"))
            .arg(Arg::from_usage("-b, --stability=[INT] 'reads to stable breakpoint'"))
            .arg(Arg::from_usage("-t, --threads=[INT] 'max threads for mash'"))
            .arg(Arg::with_name("progress").short("p").long("progress").takes_value(false).help("progress bar on"))
        )
        .subcommand(SubCommand::with_name("screen")
            .about("\nscreen read set with mash against reference sketch")
            .version("0.4.5")
            .arg(Arg::from_usage("-f, --fastx=[FILE] 'fasta/q input path'"))
            .arg(Arg::from_usage("-s, --sketch=[FILE] 'reference sketch db'"))
            .arg(Arg::from_usage("-l, --limit=[INT] 'limit output to top ranking'"))
            .arg(Arg::from_usage("-t, --threads=[INT] 'max threads for mash'"))
            .arg(Arg::with_name("pretty").short("p").long("pretty").takes_value(false).help("pretty print on"))
        )
        .get_matches();
        
    if let Some(stream) = matches.subcommand_matches("stream") {
        
        let sketch: String = stream.value_of("sketch").unwrap().to_string();
        let ranks: usize = stream.value_of("ranks").unwrap().parse::<usize>().unwrap();
        let threads: i32 = stream.value_of("threads").unwrap().parse::<i32>().unwrap();
        let stability: usize = stream.value_of("stability").unwrap().parse::<usize>().unwrap();
        let progress: bool = stream.is_present("progress")

        let (sketch_msh, genotypes): (str, str) = sketchy::get_sketch_files(&sketch);
        let (sketch_size, sketch_index): (usize, usize) = sketchy::get_sketch_info(&sketch_msh);

        sketchy::run(sketch_msh, genotypes, threads, ranks, stability, sketch_index, sketch_size, progress).map_err(
            |err| println!("{:?}", err)
        ).ok();
        
    }

    if let Some(screen) = matches.subcommand_matches("screen") {
        
        let fastx: String = screen.value_of("fastx").unwrap().to_string();
        let sketch: String = screen.value_of("sketch").unwrap().to_string();
        let genotypes: String = screen.value_of("genotypes").unwrap().to_string();
        let threads: i32 = screen.value_of("threads").unwrap().parse::<i32>().unwrap();
        let limit: usize = screen.value_of("limit").unwrap().parse::<usize>().unwrap();
        let pretty: bool = screen.is_present("pretty")

        sketchy::screen(fastx, sketch, genotypes, threads, limit).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    Ok(())

}
