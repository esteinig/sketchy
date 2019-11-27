extern crate clap;

mod sketchy;

use std::io::Error;
use clap::App;

fn main() -> Result<(), Error> {

    let matches = App::new("\nsketchy")
                          .version("0.4")
                          .about("\nNanopore lineage calling and genotyping of bacterial pathogens\n")
                          .args_from_usage(
                              "-i, --input=[FILE]  'Sets the input read list for Mash'
                               -s, --sketch=[FILE] 'Sets the sketch database to query'
                               -r, --ranks=[INT] 'Sets the ranked ssh scores to extract
                               -p, --procs=[INT] 'Sets the number of processors for Mash'")
                          .get_matches();
    
    let sketch: String = matches.value_of("sketch").unwrap().to_string();
    let reads: String = matches.value_of("input").unwrap().to_string();
    
    let procs: i32 = matches.value_of("procs").unwrap().parse::<i32>().unwrap();
    let ranks: usize = matches.value_of("ranks").unwrap().parse::<usize>().unwrap();
    

    // Get sketch_size and sketch_index length from Mash

    
    let sketch_size: usize = 1000;
    let sketch_index: usize = 38981;

    sketchy::run(sketch, reads, procs, ranks, sketch_index, sketch_size)

}
