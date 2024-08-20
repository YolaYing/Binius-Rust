use super::binary_field16::BinaryFieldElement16 as B16;
use super::binary_ntt_cache::WiEvalCache;
use lazy_static::lazy_static;
use std::fs;
use std::sync::Mutex;

lazy_static! {
    static ref WI_EVAL_CACHE: Mutex<WiEvalCache> = Mutex::new(load_or_build_wi_eval_cache());
}

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
    cache
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_wi_eval() {
        let wi_eval_cache = WI_EVAL_CACHE.lock().unwrap();
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
}
