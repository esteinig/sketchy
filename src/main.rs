extern crate cute;
extern crate clap;
extern crate indicatif;

mod sketchy;

use std::io::Error;
use clap::{Arg, App, SubCommand};

fn main() -> Result<(), Error> {

    let matches = App::new("sketchy")
        .version("0.4.5")
        .about("\nNanopore lineage calling and genotyping of bacterial pathogens using Mash\n")
        .subcommand(SubCommand::with_name("compute")
            .about("\ncompute sum of shared hashes from fasta/q on stdin")
            .version("0.4.5")
            .arg(Arg::from_usage("-r, --ranks=[INT] 'max ranks per read'"))
            .arg(Arg::from_usage("-s, --sketch=[FILE] 'reference sketch'"))
            .arg(Arg::from_usage("-p, --progress=[INT] 'progress switch > 0'"))
            .arg(Arg::from_usage("-t, --threads=[INT] 'max threads for mash'"))
        )
        .subcommand(SubCommand::with_name("evaluate")
            .about("\nevaluate sum of shared hashes from sketchy compute on stdin")
            .version("0.4.5")
            .arg(Arg::from_usage("-f, --features=[FILE] 'genotype feature index'"))
            .arg(Arg::from_usage("-s, --stable=[INT] 'reads to stable breakpoint'"))
        )
        .subcommand(SubCommand::with_name("screen")
            .about("\nscreen reads with mash against reference sketch")
            .version("0.4.5")
            .arg(Arg::from_usage("-f, --fastx=[FILE] 'fasta/q input path'"))
            .arg(Arg::from_usage("-s, --sketch=[FILE] 'reference sketch'"))
            .arg(Arg::from_usage("-f, --features=[FILE] 'genotype feature index'"))
            .arg(Arg::from_usage("-t, --threads=[INT] 'max threads for mash'"))
        )
        .get_matches();
        
    if let Some(compute) = matches.subcommand_matches("compute") {
        
        let sketch: String = compute.value_of("sketch").unwrap().to_string();
        let ranks: usize = compute.value_of("ranks").unwrap().parse::<usize>().unwrap();
        let progress: usize = compute.value_of("progress").unwrap().parse::<usize>().unwrap();
        let threads: i32 = compute.value_of("threads").unwrap().parse::<i32>().unwrap();
        
        let (sketch_size, sketch_index): (usize, usize) = sketchy::get_sketch_info(&sketch);
        
        sketchy::run(sketch, threads, ranks, sketch_index, sketch_size, progress).map_err(
            |err| println!("{:?}", err)
        ).ok();
        
    }

    if let Some(evaluate) = matches.subcommand_matches("evaluate") {
        
        let features: String = evaluate.value_of("features").unwrap().to_string();
        let stable: usize = evaluate.value_of("stable").unwrap().parse::<usize>().unwrap();

        sketchy::evaluate(features, stable).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    if let Some(screen) = matches.subcommand_matches("screen") {
        
        let fastx: String = screen.value_of("fastx").unwrap().to_string();
        let sketch: String = screen.value_of("sketch").unwrap().to_string();
        let features: String = screen.value_of("features").unwrap().to_string();
        let threads: i32 = screen.value_of("threads").unwrap().parse::<i32>().unwrap();

        let (sketch_size, sketch_index): (usize, usize) = sketchy::get_sketch_info(&sketch);

        sketchy::screen(fastx, sketch, features, threads, sketch_size, sketch_index).map_err(
            |err| println!("{:?}", err)
        ).ok();

    }

    Ok(())

}
