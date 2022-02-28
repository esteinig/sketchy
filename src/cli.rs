use std::ffi::OsStr;
use std::{ffi::OsString, path::PathBuf};
use structopt::StructOpt;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CliError {
    #[error("Scale parameter must be between 0 and 1")]
    InvalidScaleRange,
    #[error("Scale parameter must be a float")]
    InvalidScaleFloat,
}

/// Bacterial genomic neighbor typing using MinHash
#[derive(Debug, StructOpt)]
#[structopt(name = "sketchy")]
pub struct Cli {
    #[structopt(subcommand)]
    pub commands: Commands,
}

#[derive(Debug, StructOpt)]
pub enum Commands {
    /// Create a sketch from input sequences
    Sketch {
        /// Fast{a,q}.{gz,xz,bz}, stdin if not present
        #[structopt(short, long, parse(from_os_str), multiple = true)]
        input: Option<Vec<PathBuf>>,
        /// Output sketch file path.
        #[structopt(short, long, parse(from_os_str), required = true)]
        output: PathBuf,
        /// Sketch size.
        #[structopt(short, long, default_value = "1000")]
        sketch_size: usize,
        /// K-mer size.
        #[structopt(short = "k", long, default_value = "16")]
        kmer_size: u8,
        /// Hash scaler for finch format.
        #[structopt(
            short = "c", 
            long,
            parse(try_from_str = check_scale_limits),
            default_value = "0.001"
        )]
        scale: f64,
        /// Seed for hashing k-mers.
        #[structopt(short = "e", long, default_value = "0")]
        seed: u64,
    },
    /// List sketch genome order, sketch build parameters
    Info {
        /// Sketch file, format: Mash (.msh) or Finch (.fsh)
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists)
        )]
        input: PathBuf,
        /// Display the sketch build parameters.
        #[structopt(short, long)]
        params: bool,
    },

    /// Check match between sketch and genotype file
    Check {
        /// Sketch file, format: Mash (.msh) or Finch (.fsh)
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists)
        )]
        reference: PathBuf,
        /// Genotype file to validate with sketch file
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists)
        )]
        genotypes: PathBuf,
    },
    /// Compute shared hashes between two sketches
    Shared {
        /// Sketch file, format: Mash (.msh) or Finch (.fsh)
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists)
        )]
        reference: PathBuf,
        /// Sketch file, matching format: Mash (.msh) or Finch (.fsh)
        #[structopt(short, long)]
        query: PathBuf,
    },
    /// Predict genotypes from reads or read streams
    Predict {
        /// Fast{a,q}.{gz,xz,bz}, stdin if not present.
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists)
        )]
        input: Option<PathBuf>,
        /// Reference sketch, Mash (.msh) or Finch (.fsh)
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists)
        )]
        reference: PathBuf,
        /// Reference genotype table (.tsv)
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists)
        )]
        genotypes: PathBuf,
        /// Number of top ranked prediction to output
        #[structopt(short, long, default_value = "1")]
        top: usize,
        /// Number of reads to process, all reads default
        #[structopt(short, long, default_value = "0")]
        limit: usize,
        /// Sum of shared hashes per read output
        #[structopt(short, long)]
        stream: bool,
        /// Consensus prediction over top feature values
        #[structopt(short, long)]
        consensus: bool,
        /// Header added to output based on genotype file
        #[structopt(short = "H", long)]
        header: bool,
    },
}

fn check_scale_limits(scale: &str) -> Result<f64, CliError> {
    match scale.parse::<f64>() {
        Ok(x) => match x {
            x if (0.0..=1.0).contains(&x) => Ok(x),
            _ => Err(CliError::InvalidScaleRange),
        },
        _ => Err(CliError::InvalidScaleFloat),
    }
}

fn check_file_exists(file: &OsStr) -> Result<PathBuf, OsString> {
    let path = PathBuf::from(file);
    if path.exists() {
        Ok(path)
    } else {
        Err(OsString::from(format!("{:?} does not exist", path)))
    }
}
