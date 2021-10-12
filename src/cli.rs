use std::{ffi::OsString, path::PathBuf};
use structopt::StructOpt;
use thiserror::Error;
use std::ffi::OsStr;

#[derive(Error, Debug)]
pub enum CliError {
    #[error("Scale parameter must be between 0 and 1")]
    InvalidScaleRange,
    #[error("Scale parameter must be a float")]
    InvalidScaleFloat,
}

/// Predict genotypes based on cumulative shared hashes
#[derive(Debug, StructOpt)]
#[structopt(name = "sketchy")]
pub struct Cli {
    #[structopt(subcommand)]
    pub commands: Commands
}


#[derive(Debug, StructOpt)]
pub enum Commands {
    
    Sketch {
        /// Fast{a,q}.{gz,xz,bz}, stdin if not present
        #[structopt(
            short, 
            long, 
            parse(from_os_str),
            multiple = true
        )]
        input: Option<Vec<PathBuf>>,
        /// Output sketch filepath
        #[structopt(
            short, 
            long, 
            parse(from_os_str),
            required = true,
        )]
        output: PathBuf,
        /// Size of sketch .
        #[structopt(
            short, 
            long,
            default_value = "1000"
        )]
        sketch_size: usize,
        /// K-mer size
        #[structopt(
            short = "k", 
            long,
            default_value = "15"
        )]
        kmer_size: u8,
        /// Hash scaler for 'scaled' format.
        #[structopt(
            short = "c", 
            long, 
            parse(try_from_str = check_scale_limits),
            default_value = "0.001"
        )]
        scale: f64,
        /// Seed value for hashing, 42 replicates Mash.
        #[structopt(
            short = "e", 
            long, 
            default_value = "0"
        )]
        seed: u64,
    },
    Info {
        /// Sketch file, Mash (.msh) or Finch (.fsh)
        #[structopt(
            short, 
            long, 
            parse(try_from_os_str = check_file_exists)
        )]
        input: PathBuf,
        /// Display the sketch build parameters.
        #[structopt(short, long)]
        build: bool,
    },
    Check {
        /// Sketch file, Mash (.msh) or Finch (.fsh)
        #[structopt(
            short, 
            long, 
            parse(try_from_os_str = check_file_exists)
        )]
        input: PathBuf,
        /// Genotype file to validate with sketch file
        #[structopt(
            short, 
            long, 
            parse(try_from_os_str = check_file_exists)
        )]
        genotypes: PathBuf,
    },
    Shared {
        /// Sketch file, Mash (.msh) or Finch (.fsh)
        #[structopt(
            short, 
            long, 
            parse(try_from_os_str = check_file_exists)
        )]
        reference: PathBuf,
        /// Sketch file,  Mash (.msh) or Finch (.fsh)
        #[structopt(
            short, 
            long
        )]
        query: PathBuf,
    },
    Predict {
        /// Fast{a,q}.{gz,xz,bz}, stdin if not present.
        #[structopt(
            short, 
            long, 
            parse(try_from_os_str = check_file_exists)
        )]
        fastx: Option<PathBuf>,
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
        /// Number of top ranked prediction to output/
        #[structopt(
            short, 
            long,
            default_value = "1"
        )]
        top: usize,
        /// Number of reads to process, all reads by default/
        #[structopt(
            short, 
            long,
            default_value = "0"
        )]
        limit: usize,
        /// Sum of shared hashes per read output
        #[structopt(
            short, 
            long
        )]
        online: bool,
    }
}




fn check_scale_limits(scale: &str) -> Result<f64, CliError> {

    match scale.parse::<f64>() {
        Ok(x) => match x {
            x if x >= 0.0 && x <= 1.0 => Ok(x),
            _ => Err(CliError::InvalidScaleRange)
        },
        _ => Err(CliError::InvalidScaleFloat)
    }
}

fn check_file_exists(file: &OsStr) -> Result<PathBuf, OsString> {

    let path = PathBuf::from(file);
    if path.exists(){
        Ok(path)
    } else {
        Err(OsString::from(format!("{:?} does not exist", path)))
    }
    
}