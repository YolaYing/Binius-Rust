//! This module contains the prover implementation.
//! The prover is responsible for generating the proof.
//! The proof consists of the root of the Merkle tree, the evaluation, the t_prime, the columns, and the branches.
//!
//! The prover is implemented as a function that takes the evaluations and the evaluation point as input and returns the proof.
//! The prover consists of the following steps:
//! 1. Pack the evaluations into rows(the evaluations are nomarlly a list of u8s, and the rows are a list of lists of u16s).
//! 2. Extend the rows using the Fast-Fourier extension.
//! 3. Compute t_prime: linear combination of rows before extension.
//! 4. Pack columns into a Merkle tree.
//! 5. Get challenges.
//! 6. Compute evaluation.
//! 7. Return the proof.

// const EXPANSION_FACTOR: usize = 4;
// const NUM_CHALLENGES: usize = 2;
// const PACKING_FACTOR: usize = 16;
const EXPANSION_FACTOR: usize = 8;
const NUM_CHALLENGES: usize = 32;
const PACKING_FACTOR: usize = 16;

use crate::merkle_tree::get_branch;
use p3_util::log2_strict_usize;

use super::binary_field16_simd::{big_mul, uint16s_to_bits, BinaryFieldElement16};
use super::challenger::get_challenges;
use super::merkle_tree::{get_root, merkelize};
use super::utils::{
    choose_row_length_and_count, computed_tprimes, evaluation_tensor_product, extend_rows,
    pack_rows, transpose, transpose_bits, xor_along_axis,
};

pub struct Proof {
    pub root: Vec<u8>,
    pub evaluation_point: Vec<u128>,
    pub eval: Vec<u16>,
    pub t_prime: Vec<Vec<u16>>,
    pub columns: Vec<Vec<BinaryFieldElement16>>,
    pub branches: Vec<Vec<Vec<u8>>>,
}

pub fn prover(evaluations: &[u8], evalutaion_point: Vec<u128>) -> Proof {
    let log_evaluation_count = log2_strict_usize(evaluations.len() * 8);
    let (log_row_length, log_row_count, row_length, row_count) =
        choose_row_length_and_count(log_evaluation_count);

    // row packing, convert each rows into a list of BinaryFieldElement16s
    let rows = pack_rows(evaluations, row_count, row_length, PACKING_FACTOR);

    // Fast-Fourier extend the rows
    let extended_rows = extend_rows(&rows, EXPANSION_FACTOR);
    let extended_row_length = row_length * EXPANSION_FACTOR / PACKING_FACTOR;

    // Compute t_prime: linear combination of rows before extension
    let row_combination = evaluation_tensor_product(&evalutaion_point[log_row_length..].to_vec());
    assert_eq!(row_combination.len(), rows.len());
    let rows_as_bits_transpose =
        transpose_bits(rows.iter().map(|row| uint16s_to_bits(row)).collect());
    let t_prime = computed_tprimes(&rows_as_bits_transpose, &row_combination);

    // Pack columns into a Merkle tree
    let columns = transpose(&extended_rows);
    // packed_columns = [col.tobytes('C') for col in columns]
    let packed_columns = columns
        .iter()
        .map(|col| col.clone().into_iter().collect())
        .collect();
    let merkle_tree = merkelize(&packed_columns);
    let root = get_root(&merkle_tree);

    // Get challenges
    let challenges = get_challenges(&root, extended_row_length, NUM_CHALLENGES);

    // Compute evaluation
    let col_combination = evaluation_tensor_product(&evalutaion_point[..log_row_length].to_vec());
    // for each row in t_prime and each row in col_combination, use big_mul to multiply them
    let multi_result = t_prime
        .iter()
        .zip(col_combination.iter())
        .map(|(t_prime_row, col_combination_row)| big_mul(t_prime_row, col_combination_row))
        .collect::<Vec<Vec<u16>>>();
    let computed_eval = xor_along_axis(&multi_result, 0);

    Proof {
        root,
        evaluation_point: evalutaion_point,
        eval: computed_eval,
        t_prime,
        columns,
        // columns: challenges
        //     .iter()
        //     .map(|&c| columns[c as usize].clone())
        //     .collect(),
        branches: challenges
            .iter()
            .map(|c| get_branch(&merkle_tree, (*c).into()))
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prover() {
        let evaluations = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
        let evaluation_point = vec![1, 2, 3, 4, 5, 6, 7];
        let result = prover(&evaluations, evaluation_point);
        assert_eq!(
            result.root,
            vec![
                237, 216, 147, 34, 92, 247, 228, 153, 173, 169, 131, 177, 192, 61, 21, 78, 211, 52,
                232, 196, 45, 237, 95, 175, 114, 196, 110, 120, 213, 131, 64, 60
            ]
        );
        assert_eq!(result.evaluation_point, vec![1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(result.eval, vec![0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(result.t_prime[0], vec![1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(result.t_prime[1], vec![5, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(result.t_prime[2], vec![6, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(result.columns[0][0], BinaryFieldElement16::new(513));
        assert_eq!(
            result.branches[0][0],
            vec![
                93, 251, 171, 238, 223, 49, 139, 243, 60, 9, 39, 196, 61, 118, 48, 245, 27, 130,
                243, 81, 116, 3, 1, 53, 79, 163, 215, 252, 81, 240, 19, 46
            ]
        );
        assert_eq!(
            result.branches[0][1],
            vec![
                117, 231, 100, 194, 34, 254, 150, 92, 231, 217, 46, 180, 190, 247, 71, 93, 114, 98,
                57, 172, 82, 119, 48, 229, 64, 230, 118, 239, 184, 26, 154, 48
            ]
        );
    }

    #[test]
    fn test_big_data_prover() {
        let evaluations = vec![1; 1 << 20];
        let evaluation_point = vec![1; 23];
        let result = prover(&evaluations, evaluation_point);

        assert_eq!(
            result.root,
            vec![
                14, 137, 1, 182, 32, 73, 136, 127, 237, 218, 39, 11, 5, 243, 134, 95, 106, 158,
                189, 161, 93, 114, 169, 113, 24, 23, 215, 128, 16, 106, 56, 90
            ]
        );
        assert_eq!(result.evaluation_point.len(), 23);
        assert_eq!(result.eval, vec![0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(result.t_prime[0], vec![1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(result.columns.len(), 2048);
        assert_eq!(
            result.branches[7][4],
            vec![
                87, 16, 103, 115, 59, 231, 163, 189, 151, 96, 41, 109, 226, 231, 251, 42, 204, 154,
                35, 52, 8, 58, 252, 189, 51, 41, 4, 29, 30, 31, 212, 86
            ]
        );
    }
}
