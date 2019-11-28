//! @date: November 2019
//! @authors: Eike Steinig, Michael Hall
//!
//! Sketchy library module that computes the sum of shared hashes from STDOUT of MASH
//!
//! Questions:
//!
//!    - how best to output (currently println to stdout, but perhaps better to file)
//!    - how to best compile into a lib.so we could link into Python framework (plots etc.)
//!    - how are the index_size for vector init handled in that case, can they be read from MASH at runtime?
//!    - or do we need to compile a library for each reference sketch and just provide via the sketch download?
//!    - maybe keep a Rust library core with CLI
//!    - parallel implementation of online ssh compute on sequence runs <-- important

use std::io::{BufRead, BufReader, Error, ErrorKind};
use std::process::{Command, Stdio};

/// Sketchy core compute function for sum of shared hashes from `MASH`
///
/// Arguments
/// =========
///
/// sketch:
///     path to input sketch database file in Sketchy created with `MASH`
///
/// reads:
///     path to list of read files (`.fa`) for input to `MASH` (`mash dist -l mash.in`)
///
/// ranks:
///     extract the <ranks> highest ranking sum of shared hashes of each read to `STDOUT`
///
/// index_size:
///     size of sketch index, required for continuous parsing of reads from `MASH`
///
/// sketch_size:
///     sketch size used to construct sketch in `MASH`, required for fast clipping tail of line output
///
///
/// Description
/// ===========
///
/// Compute sum of shared hashes from system call `STDOUT` of `MASH`, pseudologic:
///
///     - system call to compute read file in `MASH` -> capture all reads in `STDOUT`
///     - iterate over lines (hits against sketch entries, ordered as in sketch)
///     - extract shared hashes from tail of line
///     - update sum of shared hashes by sketch index (line index)
///     - reset when lines reach end of index and new read starts
///
/// # Errors
/// todo: describe the circumstances under which this function errors, and what error(s) to expect.
///
/// # Examples
/// todo: provide code example
///
pub fn run(
    sketch: String,
    reads: String,
    procs: i32,
    ranks: usize,
    index_size: usize,
    sketch_size: usize,
) -> Result<(), Error> {
    let mash_args = [
        "dist",
        "-p",
        &*format!("{}", procs),
        "-l",
        &*format!("{}", sketch),
        &*format!("{}", reads),
    ];

    let stdout = Command::new("mash") // system call to MASH
        .args(&mash_args)
        .stdout(Stdio::piped())
        .spawn()?
        .stdout
        .ok_or_else(|| {
            Error::new(
                ErrorKind::Other,
                "Could not capture standard output from MASH.",
            )
        })?;

    // SKETCHY - SUM OF SHARED HASHES

    let reader = BufReader::new(stdout);
    let tail_index: usize = sketch_size.to_string().len(); // <tail_index> to reach shared hashes

    sum_of_shared_hashes(reader, tail_index, index_size, ranks);

    Ok(())
}

/// todo: add description of function
///
/// # Examples
/// todo: add example code snippet
fn sum_of_shared_hashes<R: BufRead>(
    reader: R,
    tail_index: usize,
    index_size: usize,
    ranks: usize,
) -> Vec<u32> {

    let mut sum_shared_hashes: Vec<u32> = vec![0; index_size];

    let mut idx: usize = 0;
    let mut read: u32 = 0; // max 4b reads

    reader
        .lines()
        .filter_map(|line| line.ok())
        .for_each(|line| {
            let shared_hashes = get_shared_hashes(&line, tail_index).parse::<u32>().unwrap();

            sum_shared_hashes[idx] += shared_hashes; // update sum of shared hashes

            // at sketch index end: output ranked ssh + reset for next read
            if idx == index_size - 1 {
                // collect the index of the current sum of shared hashes
                let mut ssh_index: Vec<(usize, &u32)> =
                    sum_shared_hashes.iter().enumerate().collect();
                // sort indexed and updated ssh scores
                ssh_index.sort_by_key(|k| k.1);
                // collect the highest ranked ssh scores and their indices
                let ranked_ssh: Vec<(usize, &u32)> =
                    ssh_index.drain(index_size - ranks..).collect();

                // write ranked ssh block for this read
                for (rank, (ix, ssh)) in ranked_ssh.iter().rev().enumerate() {
                    println!("{}\t{}\t{}\t{}", ix, ssh, rank, read); // flush everytime?
                }

                idx = 0;
                read += 1;
            } else {
                idx += 1;
            }
        });

    sum_shared_hashes
}


/// Parse shared hashes from input line string
///
/// Arguments
/// =========
///
/// input:
///     line str to extract sum of shared hashes from
///
/// tail_index:
///     shared hashes char collection starts at tail_index+1 of reverse line
///
///
/// Description
/// ===========
///
/// This function is critical for speed as splitting the line string appears slow, and the
/// reverse line string pattern much faster. It iterates over the reversed line to extract
/// the shared hashes defined as numerics after <tail_index> breaking on first non-numeric.
///
/// Shared hashes start after (not including) <tail_index> = len(sketch_size as str) ("/") in
/// reversed line; insert to front of an empty string reconstructs the correct number of
/// shared hashes as string to return.
///
/// # Examples
/// todo: add example code snippet
///
fn get_shared_hashes(input: &str, tail_index: usize) -> String {


    let mut output = String::new();
    for (i, c) in input.chars().rev().enumerate() {
        if i > tail_index {
            if c.is_numeric() {
                output.insert(0, c)
            } else {
                break;
            }
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_ssh() {
        /* Test general functionality of the compute sum of shared hashes function  */

        let sketch_size: usize = 1000;
        let default_tail_index: usize = sketch_size.to_string().len();

        let test_data_bytes = include_bytes!("data/test_mash.out");

        // https://www.reddit.com/r/rust/comments/40cte1/using_a_byte_vector_as_a_stdiobufread/
        let test_data_sliced: &[u8] = &test_data_bytes[..];

        let reader = BufReader::new(test_data_sliced); // 2 reads

        let ssh: Vec<u32> = sum_of_shared_hashes(reader, default_tail_index, 2, 1); // index_size = 2, ranks = 1

        assert_eq!(ssh, vec![7, 77]);
    }

    #[test]
    fn test_mash_dist() {
        /* Test dependency MASH in $PATH  */

        let _stdout = Command::new("mash") // system call to MASH
            .args(&["dist", "-h"])
            .output()
            .expect("Failed to run MASH dist command.");
    }

    #[test]
    fn test_compute_ssh() {
        /* Test general functionality of the compute sum of shared hashes function  */

        let sketch_size: usize = 1000;
        let default_tail_index: usize = sketch_size.to_string().len();

        let test_data_bytes = include_bytes!("data/test_mash.out");

        // https://www.reddit.com/r/rust/comments/40cte1/using_a_byte_vector_as_a_stdiobufread/
        let test_data_sliced: &[u8] = &test_data_bytes[..];

        let reader = BufReader::new(test_data_sliced); // 2 reads

        let ssh: Vec<u32> = sum_of_shared_hashes(reader, default_tail_index, 2, 1); // index_size = 2, ranks = 1

        assert_eq!(ssh, vec![7, 77]);
    }

    #[test]
    fn test_get_shared_hashes() {
        /* Test the get_shared_hashes function to extract shared hashes from line string */

        let sketch_size: usize = 1000;
        let default_tail_index: usize = sketch_size.to_string().len();

        let line_default = "ThisIs\tATestLine\t3/1000";
        assert_eq!(get_shared_hashes(line_default, default_tail_index), "3");

        let line_multiple_numeric_shared_hashes = "ThisIs\tATestLine\t42/1000";
        assert_eq!(
            get_shared_hashes(line_multiple_numeric_shared_hashes, default_tail_index),
            "42"
        );

        let line_numeric_before_sep = "ThisIs\tATestLine0\t3/1000";
        assert_eq!(
            get_shared_hashes(line_numeric_before_sep, default_tail_index),
            "3"
        );

        let line_without_preceding_sep = "ThisIs\tATestLine3/1000";
        assert_eq!(
            get_shared_hashes(line_without_preceding_sep, default_tail_index),
            "3"
        );

        let line_with_newline = "ThisIs\tATestLine\t3/1000\n";
        assert_ne!(
            get_shared_hashes(line_with_newline, default_tail_index),
            "3"
        );

        // wrong <tail_index> terminates at "/"
        assert_eq!(get_shared_hashes(line_default, 0), "100");
    }
}
