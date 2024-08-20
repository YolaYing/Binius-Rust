//! This module provides functionality for binary NTT (Number Theoretic Transform) caching.
//!
//! According to the framegraph of the vanilla version of Binius, the top 3 most time-consuming functions are:
//! 1. additive ntt: takes 60% of the total time
//! 2. computed tprimes: takes 16% of the total time(xor_along_axis takes 9% of the total time)
//! 3. inverse additive ntt: takes 7% of the total time
//!
//! So there might be 3 ways to optimize the performance:
//! 1. build a cache for Wi_eval
//! 2. additive ntt function and inverse additive ntt from recursive to iterative
//! 3. build big mul cache

use crate::binary_field16::BinaryFieldElement16 as B16;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

const MAX_DIM: usize = 4;
const MAX_SIZE: usize = 1 << MAX_DIM;

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
}
