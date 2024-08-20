mod binary_ntt_cache;
mod binary_ntt_cache_build_test;
use binary_ntt_cache::WiEvalCache;
use lazy_static::lazy_static;
use std::fs;
use std::sync::Mutex;
mod binary_field16;
mod binary_ntt;
mod challenger;
mod merkle_tree;
pub mod prover;
mod utils;
pub mod verifier;

use std::io::Result;

fn main() -> Result<()> {
    let path = "example.txt";
    let contents = "Hello, world!";

    fs::write(path, contents)?;

    println!("File written successfully!");
    Ok(())
}
