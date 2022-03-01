use anyhow::Result;
use finch::serialization::{
    read_finch_file, read_mash_file, write_finch_file, write_mash_file, Sketch,
};
use finch::sketch_schemes::{KmerCount, SketchParams};
use finch::statistics::cardinality;
use needletail::{parse_fastx_file, parse_fastx_stdin, FastxReader};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SketchyError {
    #[error("reference sketch identifier {0} does not match genotype identifier {1} at line {2}")]
    InvalidIdentifier(String, String, String),
    #[error("reference sketch and genotype table must have the same length")]
    InvalidSize,
    #[error("reference sketch file must have Mash (.msh) or Finch (.fsh) extension")]
    InvalidExtension,
    #[error("reference ({0}) {1} ({2}) does not match query ({3}) {1} ({4})")]
    InvalidSketchMatch(String, String, String, String, String),
    #[error("consensus genotype could not be computed")]
    InvalidConsensusGenotype,
    #[error("--top must be an odd number when using --consensus")]
    InvalidConsensusTop,
    #[error("failed to open file")]
    IOError(#[from] std::io::Error),
    #[error("failed to open file with Finch")]
    FinchError(#[from] finch::errors::FinchError),
    #[error("failed to open genotype file or record with CSV")]
    GenotypeTableError(#[from] csv::Error),
    #[error("failed to open Fastx file or record with Needletail")]
    FastxError(#[from] needletail::errors::ParseError),
}

/// A `Struct` used for configuring some
/// params of the predict methods
#[derive(Debug, PartialEq)]
pub struct PredictConfig {
    pub top: usize,
    pub limit: usize,
    pub stream: bool,
    pub consensus: bool,
    pub header: bool,
}

/// A `Struct` used for building the formatted
/// Sketchy reference sketch
#[derive(Debug, PartialEq)]
pub struct Sketchy {}

impl Sketchy {
    /// Create a new reference sketch instance
    pub fn new() -> Self {
        Sketchy {}
    }

    /// Prediction method for Sketchy
    ///
    /// Implement documentation.
    pub fn predict(
        &self,
        fastx: Option<PathBuf>,
        reference: PathBuf,
        genotypes: PathBuf,
        config: PredictConfig,
    ) -> Result<(), SketchyError> {
        // First check that top_results is odd when consensus is true
        if config.consensus {
            let is_odd = matches!(config.top % 2, 1);
            if !is_odd {
                return Err(SketchyError::InvalidConsensusTop);
            }
        }

        let reference_sketches = self._read_sketch(reference)?;
        let reference_params = &reference_sketches[0].sketch_params;

        let mut min_scale = 0.;
        if let Some(ref_scale) = reference_params.hash_info().3 {
            min_scale = ref_scale;
        }

        let fastx_reader = match fastx {
            Some(file) => parse_fastx_file(file)?,
            None => parse_fastx_stdin()?,
        };

        let (geno, header) = self._read_genotypes(&genotypes)?;
        let geno_map = self._genotype_hashmap(geno)?;

        // Header is printed already here, regardless of streaming or direct prediction mode

        if config.header {
            println!("reads\tsketch_id\tshared_hashes\t{}", header);
        }

        if config.stream {
            self._sum_of_shared_hashes(
                fastx_reader,
                reference_params,
                &reference_sketches,
                geno_map,
                min_scale,
                config,
            )?;
        } else {
            self._shared_hashes(
                fastx_reader,
                reference_params,
                &reference_sketches,
                geno_map,
                min_scale,
                config,
            )?;
        }

        Ok(())
    }
    /// Sketch building method for Sketchy
    ///
    /// Implement documentation.
    pub fn sketch(
        &self,
        input: Option<Vec<PathBuf>>,
        output: PathBuf,
        sketch_size: usize,
        kmer_size: u8,
        seed: u64,
        scale: f64,
    ) -> Result<(), SketchyError> {
        let files = match input {
            None => {
                let stdin = std::io::stdin();
                let files = stdin
                    .lock()
                    .lines()
                    .map(|x| PathBuf::from(x.unwrap()))
                    .collect();
                files
            }
            Some(f) => f,
        };

        let sketch_params =
            self._get_sketch_params_from_extension(&output, sketch_size, kmer_size, scale, seed)?;

        let mut writer = File::create(&output)?;
        let sketches = self._sketch_files(&sketch_params, files)?;

        match &sketch_params {
            SketchParams::Mash { .. } => {
                write_mash_file(&mut writer, &sketches)?;
            }
            SketchParams::Scaled { .. } => {
                write_finch_file(&mut writer, &sketches)?;
            }
            SketchParams::AllCounts { .. } => unreachable!(), // no AllCounts from Sketchy::new()
        }

        Ok(())
    }
    /// Information method for Sketchy
    ///
    /// Given a sketch input, list the sketch identifiers, the size of the sequence the
    /// sketch was build from (bp) and the estimated  uniqueness of the sequence (bp)
    pub fn info(&self, input: PathBuf, build: bool) -> Result<(), SketchyError> {
        let sketches = self._read_sketch(input)?;

        if build {
            let rep_sketch = &sketches[0]; // assumes all sketch params same
            let rep_params = &rep_sketch.sketch_params;

            match rep_params {
                SketchParams::Scaled {
                    kmers_to_sketch,
                    kmer_length,
                    scale,
                    hash_seed,
                } => {
                    println! {"type=scaled sketch_size={:} kmer_size={:} scale={:} seed={:}", kmers_to_sketch, kmer_length, scale, hash_seed}
                }
                SketchParams::Mash {
                    kmers_to_sketch,
                    final_size: _,
                    kmer_length,
                    no_strict: _,
                    hash_seed,
                } => {
                    println! {"type=mash sketch_size={:} kmer_size={:} seed={:}", kmers_to_sketch, kmer_length, hash_seed}
                }
                SketchParams::AllCounts { .. } => unimplemented!(),
            };
        } else {
            for sketch in &sketches {
                let kmers = &sketch.hashes;
                if let Ok(c) = cardinality(kmers) {
                    println!("{} {} {}", &sketch.name, &sketch.seq_length, c);
                }
            }
        }
        Ok(())
    }
    /// Checks that sketches in the input refernce sketch and identifiers in
    /// the genotype table are in the same order and that the reference sketch
    /// collection and genotype table are of the same size
    pub fn check(&self, input: PathBuf, genotypes: PathBuf) -> Result<(), SketchyError> {
        let sketches = self._read_sketch(input)?;
        let (genotype_data, _) = self
            ._read_genotypes(&genotypes)
            .expect("Could not read genotype file");

        // Check for same order of identifiers in sketch and genotype table
        for (i, (sketch, genotype)) in sketches.iter().zip(&genotype_data).enumerate() {
            match sketch.name == genotype[0] {
                true => continue,
                false => SketchyError::InvalidIdentifier(
                    sketch.name.to_owned(),
                    genotype[0].to_owned(),
                    i.to_string(),
                ),
            };
        }
        // Check for same size of sketch collection and genotype table
        if sketches.len() != genotype_data.len() {
            Err(SketchyError::InvalidSize)
        } else {
            println!("ok");
            Ok(())
        }
    }
    /// Compute and print shared hashes between reference and query sketches
    pub fn shared(&self, reference: PathBuf, query: PathBuf) -> Result<(), SketchyError> {
        let reference_sketches = self._read_sketch(reference)?;
        let query_sketches = self._read_sketch(query)?;

        // Scale of sketches is inferred from first sketch in file --> might need to implement
        // an empty file check which is not implemented in the finch::readers
        let mut min_scale = 0.;
        if let Some(scale1) = &query_sketches[0].sketch_params.hash_info().3 {
            if let Some(scale2) = &reference_sketches[0].sketch_params.hash_info().3 {
                min_scale = f64::min(*scale1, *scale2);
            }
        }
        // Pairwise shared hashes computation
        for ref_sketch in &reference_sketches {
            for query_sketch in &query_sketches {
                let compatible = ref_sketch
                    .sketch_params
                    .check_compatibility(&query_sketch.sketch_params);
                let common = match compatible {
                    None => Ok(self._common_hashes(
                        &ref_sketch.hashes,
                        &query_sketch.hashes,
                        min_scale,
                    )),
                    Some(incomp) => Err(SketchyError::InvalidSketchMatch(
                        ref_sketch.name.to_owned(),
                        incomp.0.to_string(),
                        incomp.1,
                        query_sketch.name.to_owned(),
                        incomp.2,
                    )),
                };
                println!(
                    "{:} {:} {:}",
                    ref_sketch.name,
                    query_sketch.name,
                    common.unwrap()
                ); // unwrap should be safe here
            }
        }
        Ok(())
    }

    fn _shared_hashes(
        &self,
        mut fastx_reader: Box<dyn FastxReader>,
        reference_params: &SketchParams,
        reference_sketches: &[Sketch],
        geno_map: HashMap<String, Vec<String>>,
        min_scale: f64,
        config: PredictConfig,
    ) -> Result<(), SketchyError> {
        // Sketcher is created for all reads
        let mut sketcher = reference_params.create_sketcher();
        // Kmers are extracted for all reads
        let mut read = 0;
        while let Some(record) = fastx_reader.next() {
            sketcher.process(record?);
            read += 1;
            if read == config.limit {
                break;
            }
        }
        // Hashed kmers and counts are extracted across all reads
        let read_hashes = sketcher.to_vec();

        let mut result_vec = vec![];
        for ref_sketch in reference_sketches {
            // Shared hashes are computed for each ref sketch
            let shared_hashes = self._common_hashes(&ref_sketch.hashes, &read_hashes, min_scale);
            result_vec.push((&ref_sketch.name, shared_hashes, &geno_map[&ref_sketch.name]));
        }
        result_vec.sort_by(|a, b| b.1.cmp(&a.1));

        self._print_results(result_vec, read, config.top, config.consensus)?;

        Ok(())
    }

    fn _sum_of_shared_hashes(
        &self,
        mut fastx_reader: Box<dyn FastxReader>,
        reference_params: &SketchParams,
        reference_sketches: &[Sketch],
        geno_map: HashMap<String, Vec<String>>,
        min_scale: f64,
        config: PredictConfig,
    ) -> Result<(), SketchyError> {
        let mut sum_of_shared_hashes = vec![0; reference_sketches.len()];
        let mut read = 1;
        while let Some(record) = fastx_reader.next() {
            let mut result_vec = vec![];
            // At each read we create a new sketcher based on the reference sketch
            let mut sketcher = reference_params.create_sketcher();
            // Kmers are then extracted for the record and hashed
            sketcher.process(record?);
            // Hashed kmers and counts are extracted for this read
            let read_hashes = sketcher.to_vec();
            // With each read, we compute the shared hashes with the reference sketch
            for (i, ref_sketch) in reference_sketches.iter().enumerate() {
                let shared_hashes =
                    self._common_hashes(&ref_sketch.hashes, &read_hashes, min_scale);
                // Finally the sum of shared hashes are updated
                sum_of_shared_hashes[i] += shared_hashes;
                result_vec.push((
                    &ref_sketch.name,
                    sum_of_shared_hashes[i],
                    &geno_map[&ref_sketch.name],
                ));
            }
            result_vec.sort_by(|a, b| b.1.cmp(&a.1));
            self._print_results(result_vec, read, config.top, config.consensus)?;
            read += 1;
            if read == config.limit + 1 {
                break;
            }
        }
        Ok(())
    }

    fn _print_results(
        &self,
        result_vec: Vec<(&String, u64, &Vec<String>)>,
        read: usize,
        top_results: usize,
        consensus: bool,
    ) -> Result<(), SketchyError> {
        if consensus {
            // For consensus calling, ignore reference genome names and shared hashes
            // and for each genotype feature, gather the calls in a new vector
            let ngenotypes = result_vec[0].2.len();
            let mut genotype_features: Vec<Vec<&String>> = vec![Vec::new(); ngenotypes];

            for (_, _, genotypes) in result_vec[..top_results].iter() {
                for (j, genotype) in genotypes.iter().enumerate() {
                    genotype_features[j].push(genotype)
                }
            }
            // For each feature consensus vector, call the most frequent value
            // as the consensus. CLI implements a strict rule for only using
            // odd --top values when using consensus calling
            let mut consensus_genotype: Vec<String> = Vec::new();
            for consensus_feature in genotype_features.iter() {
                let mut counts = HashMap::new();
                for genotype in consensus_feature {
                    *counts.entry(genotype).or_insert(0) += 1;
                }
                let consensus_value = self._get_consensus_value(counts)?;
                consensus_genotype.push(consensus_value)
            }
            println!("{:}\t-\t-\t{:}", read, consensus_genotype.join("\t"));
        } else {
            // If not computing consensus simply iterate over top results and print to console
            for (name, shared_hashes, genotype) in result_vec[..top_results].iter() {
                println!(
                    "{:}\t{:}\t{:}\t{:}",
                    read,
                    name,
                    shared_hashes,
                    genotype.join("\t")
                );
            }
        }
        Ok(())
    }

    fn _get_consensus_value(
        &self,
        counts: HashMap<&&String, usize>,
    ) -> Result<String, SketchyError> {
        let consensus_value = counts.iter().max_by(|a, b| a.1.cmp(b.1)).map(|(k, _v)| k);
        match consensus_value {
            Some(value) => Ok(value.to_string()),
            None => Err(SketchyError::InvalidConsensusGenotype),
        }
    }
    /// Analogue of `finch::distance::raw_distance` reduced to extracting common hashes
    ///
    /// Assumes hashes are sorted - not sure if need to implement the check here, in particular
    /// because we implement a read-by-read shared hashes computation in the prediction method
    /// and this may increase cost in the end. Need to test.
    fn _common_hashes(
        &self,
        ref_hashes: &[KmerCount],
        query_hashes: &[KmerCount],
        min_scale: f64,
    ) -> u64 {
        let mut i: usize = 0;
        let mut j: usize = 0;
        let mut common: u64 = 0;
        while let (Some(query), Some(refer)) = (query_hashes.get(i), ref_hashes.get(j)) {
            match query.hash.cmp(&refer.hash) {
                Ordering::Less => i += 1,
                Ordering::Greater => j += 1,
                Ordering::Equal => {
                    common += 1;
                    i += 1;
                    j += 1;
                }
            }
        }
        // At this point we've exhausted one of the two sketches, but we may have
        // more counts in the other to compare if these were scaled sketches
        if min_scale > 0. {
            let max_hash = u64::max_value() / min_scale.recip() as u64;
            while query_hashes
                .get(i)
                .map(|kmer_count| kmer_count.hash < max_hash)
                .unwrap_or(false)
            {
                i += 1;
            }
            while ref_hashes
                .get(j)
                .map(|kmer_count| kmer_count.hash < max_hash)
                .unwrap_or(false)
            {
                j += 1;
            }
        }
        common
    }

    /// Analogous method to `finch::sketch_files` excluding filtering options
    ///
    /// Filtering excluded, as we are not interested in sketching read files
    /// but assembled reference genome sequences for the genotype database.
    fn _sketch_files(
        &self,
        sketch_params: &SketchParams,
        sequence_files: Vec<PathBuf>,
    ) -> Result<Vec<Sketch>, SketchyError> {
        sequence_files
            .par_iter()
            .map(|file| {
                let mut sketcher = sketch_params.create_sketcher();
                let mut fastx_reader = parse_fastx_file(file)?;

                while let Some(record) = fastx_reader.next() {
                    sketcher.process(record?);
                }

                let sketch_hashes = sketcher.to_vec();
                let (seq_length, num_valid_kmers) = sketcher.total_bases_and_kmers();

                Ok(Sketch {
                    name: file.file_name().unwrap().to_str().unwrap().to_string(),
                    seq_length,
                    num_valid_kmers,
                    comment: "".to_string(),
                    hashes: sketch_hashes,
                    filter_params: finch::filtering::FilterParams::default(), // no filter params
                    sketch_params: sketch_params.clone(),
                })
            })
            .collect()
    }

    /// Read a sketch file into a vector of sketches based on the extension of the file
    fn _read_sketch(&self, sketch_file: PathBuf) -> Result<Vec<Sketch>, SketchyError> {
        let sketch_ext = match sketch_file.extension() {
            None => Err(SketchyError::InvalidExtension),
            Some(os_str) => match os_str.to_str() {
                Some("msh") => Ok("msh"),
                Some("fsh") => Ok("fsh"),
                _ => Err(SketchyError::InvalidExtension),
            },
        };

        let mut reader = BufReader::new(File::open(&sketch_file)?);

        match sketch_ext {
            Ok("msh") => {
                let sketches = read_mash_file(&mut reader)?;
                // Fixing the sketch param object, as the parameters is not written to file for some reason [MASH]
                let sketches_with_params = sketches
                    .iter()
                    .map(|sketch| {
                        // Rethink if really necessary, currently used only to instantiate a sketcher in the main
                        // straming function and to output a summary statistic - so technically, we only need to add
                        // the sketch size to SketchParams::Mash in the first sketch after reading into these methods
                        // rather than fixing all of them, which may slow down with very large sketch collections
                        let mut new_sketch = sketch.clone();
                        new_sketch.sketch_params = SketchParams::Mash {
                            kmers_to_sketch: sketch.hashes.len(),
                            final_size: sketch.hashes.len(),
                            kmer_length: sketch.sketch_params.k(),
                            no_strict: false,
                            hash_seed: sketch.sketch_params.hash_info().2,
                        };
                        new_sketch
                    })
                    .collect();
                Ok(sketches_with_params)
            }
            Ok("fsh") => Ok(read_finch_file(&mut reader)?),
            _ => Err(SketchyError::InvalidExtension),
        }
    }

    fn _read_genotypes(
        &self,
        genotype_file: &Path,
    ) -> Result<(Vec<Vec<String>>, String), SketchyError> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(genotype_file)?;
        let mut genotypes: Vec<Vec<String>> = vec![];
        for result in reader.records() {
            let record = result?;
            let str_vec: Vec<String> = record.iter().map(|field| field.to_string()).collect();
            genotypes.push(str_vec);
        }
        let header: Vec<String> = reader
            .headers()?
            .into_iter()
            .map(|field| field.to_string())
            .collect();
        let header_str = header[1..].join("\t"); // exclude first column identifier
        Ok((genotypes, header_str))
    }

    fn _genotype_hashmap(
        &self,
        genotypes: Vec<Vec<String>>,
    ) -> Result<HashMap<String, Vec<String>>, SketchyError> {
        let genotype_map: HashMap<String, Vec<String>> = genotypes
            .iter()
            .map(|gvec| (gvec[0].to_owned(), gvec[1..].to_owned()))
            .collect();

        Ok(genotype_map)
    }

    fn _get_sketch_params_from_extension(
        &self,
        output: &Path,
        sketch_size: usize,
        kmer_size: u8,
        scale: f64,
        seed: u64,
    ) -> Result<SketchParams, SketchyError> {
        match output.extension() {
            None => Err(SketchyError::InvalidExtension),
            Some(os_str) => match os_str.to_str() {
                Some("msh") => Ok(SketchParams::Mash {
                    kmers_to_sketch: sketch_size,
                    final_size: sketch_size,
                    no_strict: false,
                    kmer_length: kmer_size,
                    hash_seed: seed,
                }),
                Some("fsh") => Ok(SketchParams::Scaled {
                    kmers_to_sketch: sketch_size,
                    kmer_length: kmer_size,
                    scale,
                    hash_seed: seed,
                }),
                _ => Err(SketchyError::InvalidExtension),
            },
        }
    }
}
