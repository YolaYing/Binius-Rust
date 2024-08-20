use crate::binary_field16::BinaryFieldElement16 as B16;
use std::collections::HashMap;

pub struct WiEvalCache {
    cache: Vec<HashMap<B16, B16>>,
}

impl WiEvalCache {
    pub fn new() -> Self {
        WiEvalCache { cache: vec![] }
    }
}
/** calculate Wi(pt)

Definition: Wi(x) is a polynomial that returns 0 on 0...2^i-1 and 1 on 2^i
    relies on the identity W{i+1}(X) = Wi(X) * (Wi(X) + Wi(2^i))

Way to calculate Wi(x):
    1. prev = W{i-1}(x)
    2. o(x) = W{i-1}(x)(W{i-1}(x) + 1)
    3. inv_quot = o(2^i).inv() = (W{i-1}(2^i)(W{i-1}(2^i) + 1)).inv()
    4. Wi(x) = inv_quot · o(x)

Similar logic as above, the way to calculate Wi(x) evaluated at pt is that:
    1. prev  = W{i-1}(pt)
    2. prev_quot = W{i-1}(2^i)
    3. Wi(pt) = W{i-1}(pt)*(W{i-1}(pt)+1) / W{i-1}(2^i)*(W{i-1}(2^i) + 1)
        explanation: W{i-1}(pt)*(W{i-1}(pt)+1) -> o(pt),  1 / W{i-1}(2^i)*(W{i-1}(2^i) + 1) -> inv_quot
    4. check if Wi(pt) ==? W{i-1}(pt)*(W{i-1}(pt)+ W{i-1}(2^i))

Args:
    dim: the dimension of the input, dim should >= pt^2, so the max dim is u32
    pt: the point to evaluate, size of pt should <= 1<<dim
    wi_eval_cache: the cache to store the evaluations
 */
pub fn get_Wi_eval(dim: usize, pt: u16, wi_eval_cache: &mut WiEvalCache) -> B16 {
    let coord = B16::new(pt);
    while wi_eval_cache.cache.len() <= dim {
        wi_eval_cache.cache.push(HashMap::new());
    }
    if dim == 0 {
        return coord;
    }
    if !wi_eval_cache.cache[dim].contains_key(&coord) {
        let prev = get_Wi_eval(dim - 1, pt, wi_eval_cache);
        let prev_quot = get_Wi_eval(dim - 1, 1 << dim as u16, wi_eval_cache);
        let result = (prev * (prev + B16::new(1))) / (prev_quot * (prev_quot + B16::new(1)));
        wi_eval_cache.cache[dim].insert(coord.clone(), result);
    }
    wi_eval_cache.cache[dim].get(&coord).cloned().unwrap()
}

/** additive ntt: Converts a polynomial with coefficients into evaluations

Args:
    vals: the coefficients of the polynomial
    start: the start index of the polynomial
    wi_eval_cache: the cache to store the evaluations

Returns:
    the evaluations of the polynomial

Appendix: page 4-5 of https://arxiv.org/pdf/1802.03932
 */
fn additive_ntt(vals: Vec<B16>, start: usize, wi_eval_cache: &mut WiEvalCache) -> Vec<B16> {
    if vals.len() == 1 {
        return vals;
    }
    let halflen = vals.len() / 2;
    let (L, R) = vals.split_at(halflen);
    let coeff1 = get_Wi_eval(
        (halflen as f64).log2() as usize,
        start as u16,
        wi_eval_cache,
    );
    let sub_input1: Vec<_> = L
        .iter()
        .zip(R.iter())
        .map(|(i, j)| *i + *j * coeff1)
        .collect();
    let sub_input2 = sub_input1
        .iter()
        .zip(R.iter())
        .map(|(i, j)| *i + *j)
        .collect();
    let mut o = additive_ntt(sub_input1, start, wi_eval_cache);
    o.extend(additive_ntt(sub_input2, start + halflen, wi_eval_cache));
    o
}

/** inverse additive ntt: Converts evaluations into a polynomial with coefficients

Args:
    vals: the evaluations of the polynomial
    start: the start index of the polynomial
    wi_eval_cache: the cache to store the evaluations

Returns:
    the coefficients of the polynomial
*/
fn inv_additive_ntt(vals: Vec<B16>, start: usize, wi_eval_cache: &mut WiEvalCache) -> Vec<B16> {
    if vals.len() == 1 {
        return vals;
    }
    let halflen = vals.len() / 2;
    let L = inv_additive_ntt(vals[..halflen].to_vec(), start, wi_eval_cache);
    let R = inv_additive_ntt(vals[halflen..].to_vec(), start + halflen, wi_eval_cache);
    let coeff1 = get_Wi_eval(
        (halflen as f64).log2() as usize,
        start as u16,
        wi_eval_cache,
    );
    let coeff2 = coeff1 + B16::new(1);
    let mut o: Vec<_> = L
        .iter()
        .zip(R.iter())
        .map(|(i, j)| *i * coeff2 + *j * coeff1)
        .collect();
    o.append(&mut L.iter().zip(R.iter()).map(|(i, j)| *i + *j).collect());
    o
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
pub fn extend(data: Vec<B16>, expansion_factor: usize) -> Vec<B16> {
    let data = data;
    let wi_eval_cache = &mut WiEvalCache::new();
    let mut o = inv_additive_ntt(data.clone(), 0, wi_eval_cache);
    o.extend(vec![B16::new(0); data.len() * (expansion_factor - 1)]);
    additive_ntt(o, 0, wi_eval_cache)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_Wi_eval() {
        let mut wi_eval_cache = WiEvalCache::new();
        let dim = 2;
        let pt = 4;
        let result = get_Wi_eval(dim, pt, &mut wi_eval_cache);
        assert_eq!(result, B16::new(1));
    }

    #[test]
    fn test_additive_ntt() {
        let mut wi_eval_cache = WiEvalCache::new();
        let vals = vec![B16::new(1), B16::new(2), B16::new(3), B16::new(4)];
        let start = 0;
        let result = additive_ntt(vals, start, &mut wi_eval_cache);
        assert_eq!(
            result,
            vec![B16::new(1), B16::new(3), B16::new(9), B16::new(15)]
        );
    }

    #[test]
    fn test_inv_additive_ntt() {
        let mut wi_eval_cache = WiEvalCache::new();
        let vals = vec![B16::new(1), B16::new(3), B16::new(9), B16::new(15)];
        let start = 0;
        let result = inv_additive_ntt(vals, start, &mut wi_eval_cache);
        assert_eq!(
            result,
            vec![B16::new(1), B16::new(2), B16::new(3), B16::new(4)]
        );
    }

    #[test]
    fn test_extend() {
        let data = vec![B16::new(1), B16::new(3), B16::new(9), B16::new(15)];
        let expansion_factor = 2;
        let result = extend(data, expansion_factor);
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
