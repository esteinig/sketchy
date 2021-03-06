/*

@date: November 2019
@authors: Eike Steinig, Michael Hall

Sketchy computes the sum of shared hashes from STDOUT of MASH

*/

use std::fs::File;
use std::iter::FromIterator;
use std::cmp::Reverse;
use cute::c;
use indicatif::ProgressBar;
use std::collections::HashMap;
use std::process::{Command, Stdio};
use std::time::Instant;
use std::io::{BufRead, BufReader, Error, ErrorKind, stdin};

pub fn run(sketch: String, procs: i32, ranks: usize, index_size: usize, sketch_size: usize, progress: usize) -> Result<(), Error> {
    
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


    let mash_args = [
        "dist", "-p", &*format!("{}", procs), "-i", &*format!("{}", sketch), "-"
    ];

    let stdout = Command::new("mash") // system call to MASH   
        .args(&mash_args)
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(|| Error::new(ErrorKind::Other, "Could not capture standard output from MASH."))?;

    // SKETCHY - SUM OF SHARED HASHES

    let reader = BufReader::new(stdout);
    let tail_index: usize = sketch_size.to_string().len(); // <tail_index> to reach shared hashes

    compute_ssh(reader, tail_index, index_size, ranks, progress);

    Ok(())
}

#[test]
fn test_mash_dist() {
    
    /* Test dependency MASH DIST in $PATH */

    let _stdout = Command::new("mash")  
        .args(&["dist", "-h"])
        .output()
        .expect("Failed to run MASH DIST command.");

}


pub fn get_sketch_info(sketch: &String) -> (usize, usize) {
    
    /* Get sketch size and number of sketches from sketch file */

    let info = Command::new("mash")
        .args(&["info", "-H", &*format!("{}", sketch)])
        .output()
        .expect("Failed to run MASH INFO command.");

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
        .expect("Failed to run MASH INFO command.");

}

#[test]
fn test_get_sketch_info() {
    
    /* Test obtaining sketch size and number of sketches from sketch file */

    let sketch_file: String = String::from("src/data/test_mash.msh");
    let expected_values: (usize, usize) = (1000, 2);

    let return_values = get_sketch_info(&sketch_file);

    assert_eq!(return_values, expected_values);

}


fn compute_ssh<R: BufRead>(reader: R, tail_index: usize, index_size: usize, ranks: usize, progress: usize) -> Vec<u32> {
    
    /* Separated sum of shared hashes function for testing */ 

    let mut sum_shared_hashes: Vec<u32> = vec![0; index_size];

    let mut idx: usize = 0;
    let mut read: u32 = 0; // max 4b reads

    let bar = if progress > 0 {
        ProgressBar::new_spinner()
    } else {
        ProgressBar::hidden()
    };

    let start = Instant::now();
    
    reader.lines()  
        .filter_map(|line| line.ok())
        .for_each(|line| {
            
            let shared_hashes = get_shared_hashes(&line, tail_index).parse::<u32>().unwrap();

            sum_shared_hashes[idx] += shared_hashes; // update sum of shared hashes

            // at sketch index end: output ranked ssh + reset for next read
            if idx == index_size-1 {
                
                bar.tick();

                // collect the index of the current sum of shared hashes
                let mut ssh_index: Vec<(usize, &u32)> = sum_shared_hashes.iter().enumerate().collect();
                // sort indexed and updated ssh scores  
                ssh_index.sort_by_key(|k| k.1);
                // collect the highest ranked ssh scores and their indices
                let ranked_ssh: Vec<(usize, &u32)> = ssh_index.drain(index_size-ranks..).collect();

                // write ranked ssh block for this read
                for (rank, (ix, ssh)) in ranked_ssh.iter().rev().enumerate() {                    
                    println!("{}\t{}\t{}\t{}", ix, ssh, rank, read); // flush everytime?
                }

                idx = 0; read += 1;

            } else { idx += 1; }

        });

    let duration = start.elapsed();
    let msg = read.to_string() + " reads / " + &(duration.as_millis()/1000).to_string() + " seconds";

    bar.finish_with_message(&msg);
    
    sum_shared_hashes
}


#[test]
fn test_compute_ssh() {
    
    /* Test general functionality of the compute sum of shared hashes function */

    // https://www.reddit.com/r/rust/comments/40cte1/using_a_byte_vector_as_a_stdiobufread/

    let sketch_size: usize = 1000;
    let default_tail_index: usize = sketch_size.to_string().len();
    
    let test_data_bytes = include_bytes!("data/test_mash.out");
    let test_data_sliced: &[u8] = &test_data_bytes[..];

    let reader = BufReader::new(test_data_sliced); // 2 reads
    
    let ssh: Vec<u32> = compute_ssh(reader, default_tail_index, 2, 1); // index_size = 2, ranks = 1

    assert_eq!(ssh, vec![7, 77]);

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

pub fn evaluate(features: String, breakpoint: usize) -> Result<(), Error> {
    
    /* Compute the sum of ranked sums of shared hashes by feature 

    Creates an empty HashMap where keys are indices of columns in feature data
    and keys are HashMaps with indices of feature values as keys and current
    sum of ranked sums of shared hashes as values; maps are emptied on new read.
    
    Uses list comprehension macro: cute::c

    Arguments
    =========

    Compute the sum of ranked sum of shared hashes (sssh) by feature
    from the sum of shared hashes output from sketchy compute.

    Operates on STDIN : /dev/stdin

    features:
        prepared feature index for evaluation, numeric categorical 
        feature columns, row order as sketch

    breakpoint:
        denote a stable read breakpoint in the output stream after this many 
        continuous top feature value predictions, actual breakpoint is
        then: stable_breakpoint (in reads) - breakpoint (in reads)

    */

    let data_file = File::open(&features)?;
    let data_reader = BufReader::new(data_file);

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
    

    let mut read = 0;    
    for line in stdin().lock().lines() {
        // Each line in sketchy compute ssh output:
        let line_str = line?;
        let ssh_row: Vec<usize> = line_str.split("\t")
            .map(|x| x.parse::<usize>().unwrap())
            .collect();

        let idx  = ssh_row[0];
        let ssh  = ssh_row[1];
        let rank = ssh_row[2];
                
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
                    
                    let stable = evaluate_stability(&top_predictions[feature], breakpoint);

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
            read += 1;
        }
        // This needs to be after, so that at each rank = 0 the purged feature map can be properly populated
        let feature_row = &feature_data[idx];
        // Iterate mutable over feature keys
        for (key, feature_map) in sssh.iter_mut() {
            let feature_value = feature_row[*key];
            // Add ssh score to feature value in feature map, or init with 0
            let feature_value_sssh = feature_map.entry(feature_value as usize).or_insert(0);
            *feature_value_sssh += ssh;
        }
    }

    // ... this is done here:
    for (feature, fm) in sssh.iter(){
        let sorted_feature_map = get_sorted_feature_map(&fm);
        top_predictions.entry(*feature)
                        .or_insert_with(Vec::new) 
                        .push(*sorted_feature_map[0].0);
        let stable = evaluate_stability(&top_predictions[feature], breakpoint);
        let preference_score = compute_preference_score_sssh(&sorted_feature_map);
        for (feat_rank, (feat_value, sssh_score)) in sorted_feature_map.iter().enumerate() {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.8}",
                    read-1, feature, feat_value, feat_rank,
                    sssh_score, stable, preference_score
            )
        }
    }

    Ok(())

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
