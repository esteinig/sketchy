/*

@date: November 2019
@authors: Eike Steinig, Michael Hall

Sketchy computes the sum of shared hashes from STDOUT of MASH

*/

use cute::c;
use std::fs::File;
use std::path::Path;
use std::cmp::Reverse;
use std::time::Instant;
use indicatif::ProgressBar;
use std::iter::FromIterator;
use std::collections::HashMap;
use std::process::{Command, Stdio};
use prettytable::{Table, Row, Cell};
use prettytable::{Attr, color};
use prettytable::format::{FormatBuilder};
use std::io::{BufRead, BufReader, Error, ErrorKind};
use serde_json::{Value};

pub fn stream(fastx: String, sketch: String, genotype_index: String, threads: i32, reads: usize, ranks: usize, stability: usize, progress: bool, raw: bool, index_size: usize, sketch_size: usize) -> Result<(), Error> {
    
    /* Sketchy core compute function for sum of shared hashes from MASH

    Arguments
    =========

    sketch:
        path to input sketch database file in Sketchy created with MASH

    reads:
        stream (-) of fasta/q reads for input to MASH

    ranks: 
        extract the <ranks> highest ranking sum of shared hashes of each read to STDOUT

    index_size: 
        size of sketch index, required for continuous parsing of reads from MASH
    
    sketch_size: 
        sketch size used to construct sketch in MASH, required for fast clipping tail of line output

    */


    if fastx != "-" {
        if !Path::new(&fastx).exists(){
            clap::Error::with_description("Could not detect FASTX", clap::ErrorKind::InvalidValue).exit();
        };
    };

    let mash_args = [
        "dist", "-p", &*format!("{}", threads), "-i", &*format!("{}", sketch), &*format!("{}", fastx)
    ];

    let mash_dist_stream = Command::new("mash") // system call to MASH   
        .args(&mash_args)
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(|| Error::new(ErrorKind::Other, "Could not capture standard output from MASH."))?;

    // SKETCHY - SUM OF SHARED HASHES

    let mash_reader = BufReader::new(mash_dist_stream);
    let tail_index: usize = sketch_size.to_string().len(); // <tail_index> to reach shared hashes
    
    let data_file = File::open(&genotype_index)?;
    let data_reader = BufReader::new(data_file);

    sum_of_shared_hashes(mash_reader, data_reader, tail_index, index_size, reads, ranks, stability, progress, raw).map_err(
        |err| println!("{:?}", err)
    ).ok();

    Ok(())
}

pub fn screen(fastx: String, sketch: String, genotypes: String, genotype_key: String, threads: i32, limit: usize, pretty: bool) -> Result<(), Error> {
    
    /* Sketchy screening of species-wide reference sketches using `mash screen` and genomic neighbor inference

    Arguments
    =========

    fastx:
        fasta/q reads input path for mash screen

    sketch:
        path to input sketch database file in Sketchy created with MASH
    
    features:
        prepared feature index for evaluation, numeric categorical feature columns, row order as sketch

    index_size: 
        size of sketch index, required for continuous parsing of reads from MASH
    
    sketch_size: 
        sketch size used to construct sketch in MASH, required for fast clipping tail of line output

    */


    let mash_args = [
        "screen", "-p", &*format!("{}", threads), "-w", &*format!("{}", sketch), &*format!("{}", fastx)
    ];

    let screen_out = Command::new("mash") // system call to MASH   
        .args(&mash_args)
        .stderr(Stdio::null())
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(
            || Error::new(ErrorKind::Other, "Could not capture standard output from MASH SCREEN")
        )?;
    
    let screen_sorted = Command::new("sort")
        .arg("-gr")
        .stdin(screen_out)
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(|| Error::new(ErrorKind::Other, "Could not capture standard output from SORT"))?;


    let reader = BufReader::new(screen_sorted);
    
    let mut table = Table::new();
    
    if !pretty {
        let raw = FormatBuilder::new().column_separator('\t').build();
        table.set_format(raw);
    } else {
        let header_row = get_header_row(genotype_key).unwrap();
        // table.add_row(header_row); fix at genotype translation update
    }

    for (_i, line) in reader.lines().enumerate() {

        if _i >= limit {
             break   
        }

        let line = line?;
        let values: Vec<&str> = line.split_whitespace().collect();   
                
        let _identity: &str = values[0];
        let _shared_hashes: &str = values[1];

        let _sketch_id: &str = values[4];

        let _name_values: Vec<&str> = _sketch_id.split("/").collect();
        let _name: &str = _name_values.last().expect("Failed to get name from sketch reference identifier");

        let _id_values: Vec<&str> = _name.split(".").collect();
        let _id: &str = _id_values.first().expect("Failed to get unique identifier from sketch reference file name");
        
        let contents = std::fs::read_to_string(&genotypes)?;

        let _grep_results = grep(&_id, &contents); 
        let _genotype_str = _grep_results[0]; // there is only ever one unique id

        let _genotype_values: Vec<&str> = _genotype_str.split("\t").collect();

        let _screen_rank: &str = &(_i+1).to_string();

        let mut screen_row = Row::new(vec![
            Cell::new(_screen_rank),
            Cell::new(_identity),
            Cell::new(_shared_hashes)
        ]);

        for x in _genotype_values.iter() {
            screen_row.add_cell(Cell::new(x).with_style(Attr::ForegroundColor(if x == &"R" { color::RED } else { color::WHITE } )));
        }; 
        
        table.add_row(screen_row);
            
    };

    table.printstd();

    Ok(())
}

pub fn dist(fastx: String, sketch: String, genotypes: String, genotype_key: String, threads: i32, limit: usize, pretty: bool) -> Result<(), Error> {
    
    /* Sketchy screening of species-wide reference sketches using `mash dist` and genomic neighbor inference

    Arguments
    =========

    fastx:
        fasta/q reads input path for mash screen

    sketch:
        path to input sketch database file in Sketchy created with MASH
    
    features:
        prepared feature index for evaluation, numeric categorical feature columns, row order as sketch

    index_size: 
        size of sketch index, required for continuous parsing of reads from MASH
    
    sketch_size: 
        sketch size used to construct sketch in MASH, required for fast clipping tail of line output

    */


    let mash_args = [
        "dist", "-p", &*format!("{}", threads), "-r", &*format!("{}", sketch), &*format!("{}", fastx)
    ];

    let dist_out = Command::new("mash") // system call to MASH   
        .args(&mash_args)
        .stderr(Stdio::null())
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(
            || Error::new(ErrorKind::Other, "Could not capture standard output from MASH SCREEN")
        )?;
    
    let dist_sorted = Command::new("sort")
        .arg("-k5nr")
        .stdin(dist_out)
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(|| Error::new(ErrorKind::Other, "Could not capture standard output from SORT"))?;


    let reader = BufReader::new(dist_sorted);
    
    let mut table = Table::new();
    
    if !pretty {
        let raw = FormatBuilder::new().column_separator('\t').build();
        table.set_format(raw);
    } else {
        let header_row = get_header_row(genotype_key).unwrap();
        // table.add_row(header_row); fix at genotype translation update
    }

    for (_i, line) in reader.lines().enumerate() {

        if _i >= limit {
             break   
        }

        let line = line?;
        let values: Vec<&str> = line.split_whitespace().collect();   

                
        let _sketch_id: &str = values[0];
        
        let _dist: &str = values[2];
        let _shared_hashes: &str = values[4];

        let _name_values: Vec<&str> = _sketch_id.split("/").collect();
        let _name: &str = _name_values.last().expect("Failed to get name from sketch reference identifier");

        let _id_values: Vec<&str> = _name.split(".").collect();
        let _id: &str = _id_values.first().expect("Failed to get unique identifier from sketch reference file name");
        
        let contents = std::fs::read_to_string(&genotypes)?;

        let _grep_results = grep(&_id, &contents); 
        let _genotype_str = _grep_results[0]; // there is only ever one unique id

        let _genotype_values: Vec<&str> = _genotype_str.split("\t").collect();

        let _screen_rank: &str = &(_i+1).to_string();

        let mut screen_row = Row::new(vec![
            Cell::new(_screen_rank),
            Cell::new(_dist),
            Cell::new(_shared_hashes)
        ]);

        for x in _genotype_values.iter() {
            screen_row.add_cell(Cell::new(x).with_style(Attr::ForegroundColor(if x == &"R" { color::RED } else { color::WHITE } )));
        }; 
        
        table.add_row(screen_row);
            
    };

    table.printstd();

    Ok(())
}

pub fn predict(genotype_key: String, limit: usize, raw: bool) -> Result<(), Error>{

    /* Predict the genotype using either top running total match (mode = total) or last highest ranked match (mode = last)  */

    let key_file = File::open(genotype_key)?;
    let reader = BufReader::new(key_file);

    let feature_translation: HashMap<usize, Value> = serde_json::from_reader(reader)?;
    
    let mut _read_tracker: Vec<String> = vec!["0".to_string()]; // read change tracker
    let mut read_prediction: HashMap<usize, Vec<String>> = HashMap::new();

    let stdin = std::io::stdin();
    let stdin_reader = BufReader::new(stdin);

    for (_i, line) in stdin_reader.lines().enumerate() {
        
        let line = line?;
        let content: Vec<String> = line.trim().split("\t").map(
            |x| x.parse::<String>().unwrap()
        ).collect();

        let read = &content[0];

        if !_read_tracker.contains(read) {

            // At the start of a new read:

            _read_tracker[0] = read.to_string(); // reset read tracker vec
            
            // prepare the variables for genotype reconstruction
            let _values: Vec<Vec<String>> = read_prediction.values().cloned().collect();
            let _lengths: Vec<usize> = _values.iter().map(|x| x.len()).collect();
            let _max_genotype_ranks: &usize = _lengths.iter().max().unwrap();
            let _keys: Vec<usize> = read_prediction.keys().cloned().collect();
            let _max_genotype_categories: &usize = _keys.iter().max().unwrap();
            
            // iterate over genotype ranusizeks ...
            for rank in 0..*_max_genotype_ranks {

                let mut genotype: Vec<String> = vec![]; // ... start a new genotype at this rank ...
                for i in 0..*_max_genotype_categories {  // ... iterate over genotype categories ...
                    let category = &read_prediction[&i];
                    let prediction = match category.get(rank) {  // ... get prediction for this category and rank ...
                        Some(value) => value,
                        None => category.last().unwrap()  // current mode: fill with higher ranked genotypes
                    };
                    genotype.push(prediction.to_string()); // ... add prediction to genotype
                }
                let genotype_str = genotype.join("\t");

                println!("{}\t{}", &read, &genotype_str);
                
                if rank+1 >= limit { // prevent negative 
                    break
                }
            }
            
            read_prediction.clear();

            
        }

        let feature_value: usize = content[2].parse::<usize>().unwrap();
        let feature_key: usize = content[1].parse::<usize>().unwrap();

        // translate features for raw output
        let feature_data = &feature_translation[&feature_key];
        let feature_name: String = feature_data["name"].as_str().unwrap().to_string();
        let feature_prediction: &String = &feature_data["values"][feature_value].as_str().unwrap().trim().to_string();
        
        read_prediction.entry(feature_key)
            .or_insert(vec![])
            .push(feature_prediction.to_string());
        
       
        // read, feature, feat_value, feat_rank, sssh_score, stable, preference_score
        if raw {
            println!("{} {} {} {} {} {} {}", read, feature_name, feature_prediction, &content[3], &content[4], &content[5], &content[6]);
        }

    }

    Ok(())

}

pub fn display_header(genotype_key: String, pretty: bool) -> Result<(), Error> {

    let mut table = Table::new();
    
    if !pretty {
        let raw = FormatBuilder::new().column_separator('\t').build();
        table.set_format(raw);
    }

    let header_row = get_header_row(genotype_key).unwrap();

    table.add_row(header_row);

    table.printstd();

    Ok(())
}

fn get_header_row(genotype_key: String) -> Result<Row, Error> {

    let key_file = File::open(genotype_key)?;
    let reader = BufReader::new(key_file);
    
    let feature_translation: HashMap<usize, Value> = serde_json::from_reader(reader)?;
    let mut keys: Vec<usize> = feature_translation.keys().cloned().collect();

    keys.sort();
    
    let mut header_row = Row::new(vec![]);
    for key in keys.iter() {
        header_row.add_cell(
            Cell::new(feature_translation[&key]["name"].as_str().unwrap())
        );
    }

    Ok(header_row)
}

fn sum_of_shared_hashes<R: BufRead>(
    reader: R, data_reader: BufReader<File>, tail_index: usize, index_size: usize, 
    reads: usize, ranks: usize, stability: usize, progress: bool, raw: bool
) -> Result<(), Error> {
    
    /* Sum of shared hashes core function */ 

    let bar = if progress {
        ProgressBar::new_spinner()
    } else {
        ProgressBar::hidden()
    };

    let start = Instant::now();
    
    // Ranked sums of sum of shared hashes (by feature) for evaluation

    let mut feature_data = vec![];
    for (_i, line) in data_reader.lines().enumerate() {
        let line = line?;
        let vec: Vec<i32> = line.trim()
            .split("\t").map(
                |x| x.parse::<i32>().unwrap()  // catch error here: non i32 value in genotype index (i32 for -1 missing data)
            )
            .collect(); 
        feature_data.push(vec);
    }
    
    // Init the feature map

    let number_features = feature_data[0].len();
    let mut sssh: HashMap<usize, HashMap<usize, usize>> = c!{
        fidx => HashMap::new(), for fidx in 0..number_features
        };

    let mut top_predictions: HashMap<usize, Vec<usize>> = c!{
        fidx => vec![], for fidx in 0..number_features
        };
     
    
    let mut idx: usize = 0;
    let mut read: u32 = 0; // max 4b reads in stream
    let mut sum_shared_hashes: Vec<u32> = vec![0; index_size]; // max 4b ssh score in stream
  
    reader.lines()  
        .filter_map(|line| line.ok())
        .for_each(|line| {
            
            let shared_hashes = get_shared_hashes(&line, tail_index).parse::<u32>().unwrap();

            sum_shared_hashes[idx] += shared_hashes; // update sum of shared hashes

            // at sketch index end: output ranked ssh + reset for next read
            if idx == index_size-1 {
                
                bar.tick();

                if read > reads {
                    break;  // break at read limit if defined, otherwise read limit = 0
                }

                // collect the index of the current sum of shared hashes
                let mut ssh_index: Vec<(usize, &u32)> = sum_shared_hashes.iter().enumerate().collect();
                // sort indexed and updated ssh scores  
                ssh_index.sort_by_key(|k| k.1);
                // collect the highest ranked ssh scores and their indices
                let ranked_ssh: Vec<(usize, &u32)> = ssh_index.drain(index_size-ranks..).collect();
                

                // write ranked ssh block for this read
                for (rank, (ix, ssh)) in ranked_ssh.iter().rev().enumerate() {                    
                    
                    if raw {
                        println!("{}\t{}\t{}\t{}", ix, ssh, rank, read); // ssh scores
                        continue;
                    } 
                    
                    let ssh: usize = **ssh as usize;

                    // feature evaluation block
                    if rank == 0 {

                        /* Start of new read in `sketchy compute` output block after <rank> rows
                        Since we report the sum at the start of each new read (at rank == 0)
                        we need to skip the first (non-summed) read, and report again when the
                        last read block is completed and no rank == 0 can be encountered ...  */
            
                        if read > 0 {
            
                            for (feature, fm) in sssh.iter(){
            
                                // Get a sorted feature map as vector
                                let sorted_feature_map = get_sorted_feature_map(&fm);
                                // Compute Brinda et al. (2019) preference score on SSSH
                                let preference_score = compute_preference_score_sssh(&sorted_feature_map);
            
                                // Add top feature value to feature vector map:
                                top_predictions.entry(*feature)
                                               .or_insert_with(Vec::new)  // prevented by init of vecmap above
                                               .push(*sorted_feature_map[0].0);
                                
                                let stable = evaluate_stability(&top_predictions[feature], stability);
            
                                for (feat_rank, (feat_value, sssh_score)) in sorted_feature_map.iter().enumerate() {
                                    println!(
                                        "{}\t{}\t{}\t{}\t{}\t{}\t{:.8}",
                                         read-1, feature, feat_value, feat_rank,
                                         sssh_score, stable, preference_score
                                    )
                                }
                            }
                            // Clear the feature sssh scores for next read
                            for (_, feature_map) in sssh.iter_mut(){
                                feature_map.clear()    
                            }
                        }
                    }

                    // This needs to be after, so that at each rank = 0 the purged feature map can be properly populated
                    let feature_row = &feature_data[*ix];

                    // Iterate mutable over feature keys
                    for (key, feature_map) in sssh.iter_mut() {
                        let feature_value = feature_row[*key];
                        // Add ssh score to feature value in feature map, or init with 0
                        let feature_value_sssh = feature_map.entry(feature_value as usize).or_insert(0);
                        *feature_value_sssh += ssh;
                    }

                }
                // sketch index end for this read
                idx = 0; read += 1;
            
            } else { idx += 1; } // in sketch index

        });

    let duration = start.elapsed();
    let msg = read.to_string() + " reads / " + &(duration.as_millis()/1000).to_string() + " seconds";

    bar.finish_with_message(&msg);
    
    Ok(())
}


#[test]
fn test_mash_dist() {
    
    /* Test dependency MASH DIST in $PATH */

    let _stdout = Command::new("mash")  
        .args(&["dist", "-h"])
        .output()
        .expect("Failed to run MASH DIST");

}


#[test]
fn test_grep_one_result() {
    let query = "duct";
    let contents = "\
        Rust:
        safe, fast, productive.
        Pick three.";
    assert_eq!(vec!["safe, fast, productive."], grep(query, contents));
}


pub fn grep<'a>(query: &str, contents: &'a str) -> Vec<&'a str> {
    let mut results = Vec::new();

    for line in contents.lines() {
        if line.contains(query) {
            results.push(line);
        }
    }

    results
}

pub fn get_sketch_files(db: String, sketchy_path: &String)  -> (String, String, String, String) {
    
    /* Get sketch files from database path and perform checks */

    let db_path = Path::new(&db);
    let db_name = db_path.file_name().unwrap().to_str().unwrap();
    
    let db_path = if !db_path.exists() {
        Path::new(&sketchy_path).join(db_name)
    }  else {
        Path::new(&sketchy_path).join(db_name)
    };
    
    if !db_path.exists() {
        clap::Error::with_description("Database sketch directory is not available", clap::ErrorKind::InvalidValue).exit();
    };

    let db_sketch = db_path.join(
        format!("{}.msh", db_name)
    );
    let db_genotypes = db_path.join(
        format!("{}.tsv", db_name)
    );
    let db_index = db_path.join(
        format!("{}.idx", db_name)
    );
    let db_key = db_path.join(
        format!("{}.key", db_name)
    );

    if !db_sketch.exists(){
        clap::Error::with_description("Database sketch is missing sketch file (.msh)", clap::ErrorKind::InvalidValue).exit();
    };
    if !db_genotypes.exists(){
        clap::Error::with_description("Database sketch is missing genotype file (.tsv)", clap::ErrorKind::InvalidValue).exit();
    };
    if !db_index.exists(){
        clap::Error::with_description("Database sketch is missing index file (.idx)", clap::ErrorKind::InvalidValue).exit();
    };
    if !db_key.exists(){
        clap::Error::with_description("Database sketch is missing key file (.key)", clap::ErrorKind::InvalidValue).exit();
    };

    (
        db_sketch.to_str().unwrap().to_string(),
        db_genotypes.to_str().unwrap().to_string(),
        db_index.to_str().unwrap().to_string(),
        db_key.to_str().unwrap().to_string()
    )
    
}


pub fn get_sketch_info(sketch: &String) -> (usize, usize) {
    
    /* Get sketch size and number of sketches from sketch file */

    let info = Command::new("mash")
        .args(&["info", "-H", &*format!("{}", sketch)])
        .output()
        .expect("Failed to run MASH INFO");

    let info_lines = std::str::from_utf8(&info.stdout).unwrap().lines();

    let mut return_values = vec![];    
    for (i, line) in info_lines.enumerate() {
        if i == 4 || i == 5 {
            let values: Vec<&str> = line.split_whitespace().collect();            
            let value: usize = values.last().unwrap().parse().unwrap();
            return_values.push(value);
        }
    }

    (return_values[0], return_values[1])
}

#[test]
fn test_mash_info() {
    
    /* Test dependency MASH INFO in $PATH */

    let _info = Command::new("mash")
        .args(&["info", "-h"])
        .output()
        .expect("Failed to run MASH INFO");

}

#[test]
fn test_get_sketch_info() {
    
    /* Test obtaining sketch size and number of sketches from sketch file */

    let sketch_file: String = String::from("src/data/test_mash.msh");
    let expected_values: (usize, usize) = (1000, 2);

    let return_values = get_sketch_info(&sketch_file);

    assert_eq!(return_values, expected_values);

}

fn get_shared_hashes(input: &str, tail_index: usize) -> String {

    /* Parse shared hashes from input line string

    Arguments
    =========

    input: 
        line str to extract sum of shared hashes from

    tail_index: 
        shared hashes char collection starts at tail_index+1 of reverse line
    

    Description
    ===========

    This function is critical for speed as splitting the line string appears slow, and the 
    reverse line string pattern much faster. It iterates over the reversed line to extract 
    the shared hashes defined as numerics after <tail_index> breaking on first non-numeric.

    Shared hashes start after (not including) <tail_index> = len(sketch_size as str) ("/") in 
    reversed line; insert to front of an empty string reconstructs the correct number of 
    shared hashes as string to return.

    */

    let mut output = String::new();
    for (i, c) in input.chars().rev().enumerate() {
        if i > tail_index  {
            if c.is_numeric()  {
                output.insert(0, c)
            } else {
                break;
            }
        }

    }
    output
}


#[test]
fn test_get_shared_hashes() {

    /* Test the get_shared_hashes function to extract shared hashes from line string */

    let sketch_size: usize = 1000;
    let default_tail_index: usize = sketch_size.to_string().len();

    let line_default = "ThisIs\tATestLine\t3/1000";
    assert_eq!(get_shared_hashes(line_default, default_tail_index), "3");

    let line_multiple_numeric_shared_hashes = "ThisIs\tATestLine\t42/1000";
    assert_eq!(get_shared_hashes(line_multiple_numeric_shared_hashes, default_tail_index), "42"); 

    let line_numeric_before_sep = "ThisIs\tATestLine0\t3/1000";
    assert_eq!(get_shared_hashes(line_numeric_before_sep, default_tail_index), "3");
    
    let line_without_preceding_sep = "ThisIs\tATestLine3/1000";
    assert_eq!(get_shared_hashes(line_without_preceding_sep, default_tail_index), "3");
    
    let line_with_newline = "ThisIs\tATestLine\t3/1000\n";
    assert_ne!(get_shared_hashes(line_with_newline, default_tail_index), "3"); 

    // wrong <tail_index> terminates at "/"
    assert_eq!(get_shared_hashes(line_default, 0), "100");

}

fn compute_preference_score_sssh(feature_map: &Vec<(&usize, &usize)>) -> f64 {
    if feature_map.len() < 2 {
        1.0
    } else {
        let primary = *feature_map[0].1 as f64;
        let secondary = *feature_map[1].1 as f64;

        if secondary == 0.0 {
            0.0
        } else {
            let numerator = 2.0 * primary;
            let denominator = primary + secondary;

            (numerator / denominator) - 1.0
        }
    }

}

fn evaluate_stability(top_vec: &Vec<usize>, breakpoint: usize) -> usize {

    let len = top_vec.len();

    if len < breakpoint {
        0
    } else {

        let mut last_predictions: Vec<&usize> = top_vec.iter().rev().take(breakpoint).collect();
        last_predictions.dedup();
        let uniques = last_predictions.len();

        if uniques == 1 {
            1
        } else {
            0
        }

    }

}

fn get_sorted_feature_map(fm: &HashMap<usize, usize>) -> Vec<(&usize, &usize)> {
        
    let mut feature_map = Vec::from_iter(fm);
    feature_map.sort_by_key(|k| Reverse(k.1));

    return feature_map
}

// #[test]
// fn test_ranked_sum_of_shared_hashes() {
    
//     /* Test general functionality of the compute sum of shared hashes function */

//     // https://www.reddit.com/r/rust/comments/40cte1/using_a_byte_vector_as_a_stdiobufread/

//     let sketch_size: usize = 1000;
//     let default_tail_index: usize = sketch_size.to_string().len();
    
//     let test_data_bytes = include_bytes!("data/test_mash.out");
//     let test_data_sliced: &[u8] = &test_data_bytes[..];

//     let reader = BufReader::new(test_data_sliced); // 2 reads
    
//     let ssh: Vec<u32> = ranked_sum_of_shared_hashes(reader, default_tail_index, 2, 1); // index_size = 2, ranks = 1

//     assert_eq!(ssh, vec![7, 77]);

// }


// fn compute_preference_score(primary: f64, secondary: f64) -> f64 {
    
//     if secondary == 0.0 {
//         0.0
//     } else {
//         let numerator = 2.0 * primary;
//         let denominator = primary + secondary;

//         (numerator / denominator) - 1.0
//     }

// }
