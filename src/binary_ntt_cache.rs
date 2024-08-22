//! This module provides functionality for binary NTT (Number Theoretic Transform) caching.
//!
//! According to the framegraph of the vanilla version of Binius, the top 3 most time-consuming functions are:
//! 1. additive ntt: takes 60% of the total time
//! 2. computed tprimes: takes 16% of the total time(xor_along_axis takes 9% of the total time)
//! 3. inverse additive ntt: takes 7% of the total time
//!
//! So there might be 3 ways to optimize the performance:
//! 1. build a cache for Wi_eval(in this file)
//! 2. additive ntt function and inverse additive ntt from recursive to iterative
//! 3. build big mul cache(not work)

use super::binary_field16::BinaryFieldElement16 as B16;
use lazy_static::lazy_static;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::sync::Mutex;

lazy_static! {
    pub static ref WI_EVAL_CACHE: Mutex<WiEvalCache> = Mutex::new(load_or_build_wi_eval_cache());
}
const MAX_DIM: usize = 16;
const MAX_SIZE: usize = 1 << MAX_DIM;

fn load_or_build_wi_eval_cache() -> WiEvalCache {
    let cache_file = "wi_eval_cache.json";
    if let Ok(data) = fs::read_to_string(cache_file) {
        if let Ok(cache) = serde_json::from_str(&data) {
            return cache;
        }
    }

    let mut cache = WiEvalCache::new();
    cache.build_Wi_eval_cache();
    match serde_json::to_string(&cache) {
        Ok(data) => {
            println!("Writing cache to file");
            if let Err(e) = fs::write(cache_file, data) {
                eprintln!("Failed to write cache to file: {}", e);
            }
        }
        Err(e) => {
            eprintln!("Failed to serialize cache: {}", e);
        }
    }
    println!("wi_eval_cache built");
    cache
}

#[derive(Serialize, Deserialize)]
pub struct WiEvalCache {
    cache: Vec<HashMap<B16, B16>>,
}

impl WiEvalCache {
    pub fn new() -> Self {
        WiEvalCache { cache: vec![] }
    }

    pub fn build_Wi_eval_cache(&mut self) -> &mut Self {
        let mut Wi_eval_cache = vec![HashMap::new(); MAX_DIM];
        // for wi_eval_cache[0], for all key, value = key
        for pt in 0..MAX_SIZE {
            Wi_eval_cache[0].insert(B16::new(pt as u16), B16::new(pt as u16));
        }
        for dim in 1..MAX_DIM {
            let prev = Wi_eval_cache[dim - 1].clone();
            let prev_quot = Wi_eval_cache[dim - 1]
                .get(&B16::new(1 << dim))
                .cloned()
                .unwrap();
            let inv_quot = (prev_quot * (prev_quot + B16::new(1))).inv();
            // for each element in prev, get the value and conduct (value_prev_element * (value_prev_element + B16::new(1))) * inv_quot
            let mut result = HashMap::new();
            for (key, value) in prev.iter() {
                result.insert(
                    key.clone(),
                    (value.clone() * (value.clone() + B16::new(1))) * inv_quot,
                );
            }

            Wi_eval_cache[dim] = result;
        }
        self.cache = Wi_eval_cache;
        self
    }

    pub fn get_Wi_eval(&self, dim: usize, pt: u16) -> B16 {
        let coord = B16::new(pt);
        if dim == 0 {
            return coord;
        }
        self.cache[dim].get(&coord).cloned().unwrap()
    }
}

/** additive ntt: Converts a polynomial with coefficients into evaluations

in the Binius, it used in the extension of the rows. when we have transform the original row into coefficients,
    we can use the additive ntt to convert the extended row(row_length * EXPANSION_FACTOR) into evaluations

Args:
    vals: the coefficients of the polynomial
    start: the start index of the polynomial, reserved for the recursive call
    wi_eval_cache: the cache to store the evaluations

Returns:
    the evaluations of the polynomial

Appendix: page 4-5 of https://arxiv.org/pdf/1802.03932
 */
// original recursive version
// fn additive_ntt(vals: &Vec<B16>, start: usize) -> Vec<B16> {
//     if vals.len() == 1 {
//         return vec![vals[0]];
//     }
//     let halflen = vals.len() / 2;
//     let (L, R) = vals.split_at(halflen);
//     // coeff1 = W{i}(start), i = log2(halflen)
//     let coeff1 = {
//         let wi_eval_cache = WI_EVAL_CACHE.lock().unwrap();
//         wi_eval_cache.get_Wi_eval((halflen as f64).log2() as usize, start as u16)
//     };

//     // sub_input1 = L + R * coeff1
//     let sub_input1: Vec<_> = L
//         .iter()
//         .zip(R.iter())
//         .map(|(i, j)| *i + *j * coeff1)
//         .collect();
//     // sub_input2 = L + R
//     let sub_input2 = sub_input1
//         .iter()
//         .zip(R.iter())
//         .map(|(i, j)| *i + *j)
//         .collect();

//     // o = additive_ntt(sub_input1, start) + additive_ntt(sub_input2, start + halflen)
//     let mut o = additive_ntt(&sub_input1, start);
//     o.extend(additive_ntt(&sub_input2, start + halflen));
//     o
// }

// Optimized iterative version: save 46% of the time
fn additive_ntt(vals: &Vec<B16>, start: usize) -> Vec<B16> {
    let mut results = vals.clone();
    let size = results.len();
    let mut step = size;

    while step >= 2 {
        step >>= 1;
        let halflen = step;

        for i in (0..size).step_by(step * 2) {
            let coeff1 = {
                let wi_eval_cache = WI_EVAL_CACHE.lock().unwrap();
                wi_eval_cache.get_Wi_eval((halflen as f64).log2() as usize, (start + i) as u16)
            };

            for j in 0..halflen {
                let l = results[i + j];
                let r = results[i + j + halflen];
                let sub_input1 = l + r * coeff1;
                results[i + j] = sub_input1;
                results[i + j + halflen] = sub_input1 + r;
            }
        }
    }

    results
}

/** inverse additive ntt: Converts evaluations into a polynomial with coefficients

Args:
    vals: the evaluations of the polynomial
    start: the start index of the polynomial
    wi_eval_cache: the cache to store the evaluations

Returns:
    the coefficients of the polynomial
*/

// original recursive version
// fn inv_additive_ntt(vals: Vec<B16>, start: usize) -> Vec<B16> {
//     if vals.len() == 1 {
//         return vals;
//     }
//     let halflen = vals.len() / 2;
//     // L = inv_additive_ntt(vals[..halflen], start)
//     let L = inv_additive_ntt(vals[..halflen].to_vec(), start);
//     // R = inv_additive_ntt(vals[halflen..], start + halflen)
//     let R = inv_additive_ntt(vals[halflen..].to_vec(), start + halflen);
//     // coeff1 = W{i}(start), i = log2(halflen)
//     let coeff1 = {
//         let wi_eval_cache = WI_EVAL_CACHE.lock().unwrap();
//         wi_eval_cache.get_Wi_eval((halflen as f64).log2() as usize, start as u16)
//     };
//     // coeff2 = coeff1 + 1
//     let coeff2 = coeff1 + B16::new(1);
//     // o = [L * coeff2 + R * coeff1] + [L + R]
//     let mut o: Vec<_> = L
//         .iter()
//         .zip(R.iter())
//         .map(|(i, j)| *i * coeff2 + *j * coeff1)
//         .collect();
//     o.append(&mut L.iter().zip(R.iter()).map(|(i, j)| *i + *j).collect());
//     o
// }

// Optimized iterative version: save 15% of the time
fn inv_additive_ntt(vals: &Vec<B16>, start: usize) -> Vec<B16> {
    let size = vals.len();
    if size == 1 {
        return vals.clone();
    }

    let mut results = vals.clone();
    let mut step = 1;
    while step < size {
        let halflen = step;
        step <<= 1;

        for i in (0..size).step_by(step) {
            // 获取系数
            let coeff1 = {
                let wi_eval_cache = WI_EVAL_CACHE.lock().unwrap();
                wi_eval_cache.get_Wi_eval((halflen as f64).log2() as usize, (start + i) as u16)
            };
            let coeff2 = coeff1 + B16::new(1);

            for j in 0..halflen {
                let l = results[i + j];
                let r = results[i + j + halflen];
                let sub_input1 = l * coeff2 + r * coeff1;
                let sub_input2 = l + r;
                results[i + j] = sub_input1;
                results[i + j + halflen] = sub_input2;
            }
        }
    }

    results
}

/** Reed-Solomon extension, using the efficient algorithms above

the logic of the function is that:
    first use inv_additive_ntt to convert the row into coefficients,
    then expend the row by expansion_factor times(e,g, 2 times),
    then padding the row with 0s after the original row(e.g. expansion_factor - 1 times)
    then use the additive_ntt to convert the row into evaluations

Args:
    data: the coefficients of the polynomial, one row of the matrix before extension
    expansion_factor: the expansion factor

Returns:
    the coefficients of the extended polynomial
*/
// pub fn extend(data: &Vec<B16>, expansion_factor: usize) -> Vec<B16> {
//     let data = data;
//     let mut o = inv_additive_ntt(data.clone(), 0);
//     o.extend(vec![B16::new(0); data.len() * (expansion_factor - 1)]);
//     additive_ntt(&o, 0)
// }
pub fn extend(data: &Vec<B16>, expansion_factor: usize) -> Vec<B16> {
    // Avoid unnecessary clone by passing reference
    let mut o = inv_additive_ntt(data, 0);

    // Calculate the total length after expansion
    let total_len = data.len() * expansion_factor;

    // Pre-allocate the extended vector with the required capacity
    o.reserve(total_len - o.len());

    // Extend the vector with zeros
    o.extend((0..(total_len - o.len())).map(|_| B16::new(0)));

    additive_ntt(&o, 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_Wi_eval() {
        let mut wi_eval_cache = WiEvalCache::new();
        wi_eval_cache.build_Wi_eval_cache();
        let dim = 2;
        let pt = 4;
        let result = wi_eval_cache.get_Wi_eval(dim, pt);
        assert_eq!(result, B16::new(1));
    }

    #[test]
    fn test_cache_file_creation() {
        let cache_file = "wi_eval_cache.json";
        // Ensure the cache file is created
        assert!(std::path::Path::new(cache_file).exists());
    }

    #[test]
    fn test_additive_ntt() {
        let vals = vec![B16::new(1), B16::new(2), B16::new(3), B16::new(4)];
        let result = additive_ntt(&vals, 0);
        assert_eq!(
            result,
            vec![B16::new(1), B16::new(3), B16::new(9), B16::new(15)]
        );
    }

    #[test]
    fn test_inv_additive_ntt() {
        let vals = vec![B16::new(1), B16::new(3), B16::new(9), B16::new(15)];
        let result = inv_additive_ntt(&vals, 0);
        assert_eq!(
            result,
            vec![B16::new(1), B16::new(2), B16::new(3), B16::new(4)]
        );
    }

    #[test]
    fn test_extend() {
        let data = vec![B16::new(1), B16::new(3), B16::new(9), B16::new(15)];
        let expansion_factor = 2;
        let result = extend(&data, expansion_factor);
        assert_eq!(
            result,
            vec![
                B16::new(1),
                B16::new(3),
                B16::new(9),
                B16::new(15),
                B16::new(14),
                B16::new(15),
                B16::new(14),
                B16::new(11)
            ]
        );
    }
}
