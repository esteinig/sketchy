use anyhow::Result;
use structopt::StructOpt;
use crate::cli::Cli;
use crate::cli::Commands::{Sketch, Predict, Info, Shared, Check};
use crate::sketchy::Sketchy;

mod cli;
mod sketchy;

/// Sketchy application
///
/// Run the application from arguments provided
/// by the command line interface
///
/// Hash seed by default is 0; for hashing to 
/// replicate Mash, seed must be 42.
fn main() -> Result<()> {    

    let args = Cli::from_args();
    let sketchy = Sketchy::new();

    match args.commands {
        Sketch { input, output, sketch_size, kmer_size, scale, seed } => {
            sketchy.sketch(input, output, sketch_size, kmer_size, seed, scale)?;
        },
        Info { input, params } => {
            sketchy.info(input, params)?;
        },
        Shared { reference, query} => {
            sketchy.shared(reference, query)?;
        },
        Predict { input, reference, genotypes, top, limit, stream } => {
            sketchy.predict(input, reference, genotypes, top, limit, stream)?;
        }
        Check { input, genotypes } => {
            sketchy.check(input, genotypes)?;
        }
    }

    Ok(())
}