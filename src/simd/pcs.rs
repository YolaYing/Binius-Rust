const EXPANSION_FACTOR: usize = 8;
const NUM_CHALLENGES: usize = 32;
const PACKING_FACTOR: usize = 16;

use std::str;

use crate::simd::binary_field16_simd::int_to_bigbin;

use super::merkle_tree::get_branch;
use p3_util::log2_strict_usize;

use super::binary_field16_simd::{big_mul, uint16_to_bit, uint16s_to_bits, BinaryFieldElement16};
use super::challenger::get_challenges;
use super::merkle_tree::{get_root, merkelize, verify_branch};
use super::utils::{
    choose_row_length_and_count, computed_tprimes, evaluation_tensor_product, extend_rows,
    multisubset, pack_row, pack_rows, transpose, transpose_3d, transpose_bits, xor_along_axis,
};

pub struct Commitment {
    pub root: Vec<u8>,
    pub packed_columns: Vec<Vec<u8>>,
    pub merkle_tree: Vec<Vec<u8>>,
    pub rows: Vec<Vec<BinaryFieldElement16>>,
    pub columns: Vec<Vec<BinaryFieldElement16>>,
}

pub struct Proof {
    pub evaluation_point: Vec<u128>,
    pub eval: u128,
    pub t_prime: Vec<u128>,
    pub columns: Vec<Vec<BinaryFieldElement16>>,
    pub branches: Vec<Vec<Vec<u8>>>,
}

pub fn commit(evaluations: &[u8]) -> Commitment {
    let log_evaluation_count = log2_strict_usize(evaluations.len() * 8);
    let (log_row_length, log_row_count, row_length, row_count) =
        choose_row_length_and_count(log_evaluation_count);

    // row packing, convert each rows into a list of BinaryFieldElement16s
    let rows = pack_rows(evaluations, row_count, row_length, PACKING_FACTOR);

    // Fast-Fourier extend the rows
    let extended_rows = extend_rows(&rows, EXPANSION_FACTOR);
    let extended_row_length = row_length * EXPANSION_FACTOR / PACKING_FACTOR;

    // Pack columns into a Merkle tree
    let columns = transpose(&extended_rows);
    // packed_columns = [col.tobytes('C') for col in columns]
    let packed_columns: Vec<Vec<u8>> = columns
        .iter()
        .map(|col| col.iter().copied().collect())
        .collect();
    let merkle_tree = merkelize(&packed_columns);
    let root = get_root(&merkle_tree);

    Commitment {
        root,
        packed_columns,
        merkle_tree,
        rows,
        columns,
    }
}

pub fn prove(commitment: &Commitment, evaluations: &[u8], evaluation_point: &Vec<u128>) -> Proof {
    let log_evaluation_count = log2_strict_usize(evaluations.len() * 8);
    let (log_row_length, log_row_count, row_length, row_count) =
        choose_row_length_and_count(log_evaluation_count);
    let extended_row_length = row_length * EXPANSION_FACTOR / PACKING_FACTOR;

    // Compute t_prime: linear combination of rows before extension
    let row_combination = evaluation_tensor_product(&evaluation_point[log_row_length..].to_vec());
    assert_eq!(row_combination.len(), commitment.rows.len());
    let rows_as_bits_transpose = transpose_bits(
        commitment
            .rows
            .iter()
            .map(|row| uint16s_to_bits(row))
            .collect(),
    );
    let t_prime = computed_tprimes(&rows_as_bits_transpose, &row_combination);

    // Get challenges
    let challenges = get_challenges(&commitment.root, extended_row_length, NUM_CHALLENGES);

    // Compute evaluation
    let col_combination = evaluation_tensor_product(&evaluation_point[..log_row_length].to_vec());
    // for each row in t_prime and each row in col_combination, use big_mul to multiply them
    let computed_eval = t_prime
        .iter()
        .zip(col_combination.iter())
        .fold(0u128, |acc, (&t, &c)| acc ^ big_mul(t, c));

    Proof {
        evaluation_point: evaluation_point.clone(),
        eval: computed_eval,
        t_prime,
        columns: challenges
            .iter()
            .map(|&c| commitment.columns[c as usize].clone())
            .collect(),
        branches: challenges
            .iter()
            .map(|c| get_branch(&commitment.merkle_tree, (*c).into()))
            .collect(),
    }
}

pub fn verifier(commitment: &Commitment, proof: &Proof, evaluation_point: &Vec<u128>) -> bool {
    let columns = &commitment.packed_columns;
    let evaluation_point = &proof.evaluation_point;
    let value = &proof.eval;
    let t_prime = &proof.t_prime;
    let root = &commitment.root;
    let branches = &proof.branches;

    // Compute the row length and row count of the grid. Should output same numbers as what prover gave
    let (log_row_length, log_row_count, row_length, row_count) =
        choose_row_length_and_count(evaluation_point.len());
    let extended_row_length = row_length * EXPANSION_FACTOR / PACKING_FACTOR;

    // Compute challenges. Should output the same as what prover computed
    let challenges = get_challenges(&root, extended_row_length, NUM_CHALLENGES);

    // Verify Merkle branches
    for i in 0..NUM_CHALLENGES {
        let challenge = challenges[i];
        let packed_column: Vec<u8> = columns[challenge as usize].clone().into_iter().collect();
        let branch = branches[i].clone();
        assert!(verify_branch(
            &root,
            challenge as usize,
            &packed_column,
            &branch
        ));
    }

    // Use the same Reed-Solomon code that the prover used to extend the rows,
    // but to extend t_prime. We do this separately for each bit of t_prime
    // try to use u8
    // let t_prime_bits: Vec<Vec<u8>> = t_prime
    //     .iter()
    //     .map(|row| row.to_le_bytes().to_vec())
    //     .collect();
    let t_prime_bits: Vec<Vec<u8>> = t_prime
        .iter()
        .map(|&row| (0..128).map(|i| ((row >> i) & 1) as u8).collect())
        .collect();
    // transpose the bits
    let t_prime_bits_transpose = transpose_bits(t_prime_bits);
    // pack the each row of t_prime_bits_transpose into a list of BinaryFieldElement16s
    let t_prime_columns: Vec<Vec<BinaryFieldElement16>> = t_prime_bits_transpose
        .iter()
        .map(|row| pack_row(row, t_prime_bits_transpose[0].len() * 8, PACKING_FACTOR))
        .collect();
    // extend the rows
    let extended_t_prime_columns = extend_rows(&t_prime_columns, EXPANSION_FACTOR);
    // Convert our FFT-extended t_prime rows into bits
    // step 1: use challenge to select columns, and convert to bits
    let extended_t_prime_columns_slices: Vec<Vec<Vec<BinaryFieldElement16>>> =
        extended_t_prime_columns
            .iter()
            .map(|row| challenges.iter().map(|&c| vec![row[c as usize]]).collect())
            .collect();
    let extended_t_prime_bits: Vec<Vec<Vec<u8>>> = extended_t_prime_columns_slices
        .iter()
        .map(|row| row.iter().map(|uint16| uint16s_to_bits(uint16)).collect())
        .collect();
    // step 2: transpose the bits
    let extended_t_prime_bits_transpose = transpose_3d(&extended_t_prime_bits, (1, 2, 0));

    // Here, we take advantage of the linearity of the code. A linear combination of the Reed-Solomon extension gives the same result as an extension of the linear combination.
    let row_combination = evaluation_tensor_product(&evaluation_point[log_row_length..].to_vec());
    let selected_columns = proof.columns.clone();
    // Each column is a vector of row_count uint16's. Convert each uint16 into bits
    let column_bits: Vec<Vec<Vec<u8>>> = selected_columns
        .iter()
        .map(|col| col.iter().map(|uint16| uint16_to_bit(uint16)).collect())
        .collect();
    // Take the same linear combination the prover used to compute t_prime, and apply it to the columns of bits.
    let transposed_column_bits = transpose_3d(&column_bits, (0, 2, 1));
    // tranform the row_combination into a list of uint16s
    let row_combination = row_combination.iter().map(|&x| int_to_bigbin(x)).collect();
    let computed_tprimes = multisubset(&row_combination, &transposed_column_bits);
    // Turn the computed tprimes into bits using uint16s_to_bits
    let computed_tprime_bits: Vec<Vec<Vec<u8>>> = computed_tprimes
        .iter()
        .map(|row| row.iter().map(|uint16| uint16s_to_bits(uint16)).collect())
        .collect();
    // The bits of the t_prime extension should equal the bits of the row linear combination of the column bits
    assert_eq!(computed_tprime_bits, extended_t_prime_bits_transpose);

    // Compute the evaluation
    let col_combination = evaluation_tensor_product(&evaluation_point[..log_row_length].to_vec());
    let computed_eval = t_prime
        .iter()
        .zip(col_combination.iter())
        .fold(0u128, |acc, (&t, &c)| acc ^ big_mul(t, c));
    assert_eq!(computed_eval, *value);
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_commit() {
        let evaluations = vec![1; 1 << 20];
        let result = commit(&evaluations);

        assert_eq!(
            result.root,
            vec![
                14, 137, 1, 182, 32, 73, 136, 127, 237, 218, 39, 11, 5, 243, 134, 95, 106, 158,
                189, 161, 93, 114, 169, 113, 24, 23, 215, 128, 16, 106, 56, 90
            ]
        );
    }

    #[test]
    fn test_prove() {
        let evaluations = vec![1; 1 << 20];
        let commitment = commit(&evaluations);
        let evaluation_point = vec![1; 23];
        let result = prove(&commitment, &evaluations, &evaluation_point);

        assert_eq!(result.evaluation_point.len(), 23);
        // assert_eq!(result.eval, vec![0, 0, 0, 0, 0, 0, 0, 0]);
        // assert_eq!(result.t_prime[0], vec![1, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(
            result.branches[7][4],
            vec![
                87, 16, 103, 115, 59, 231, 163, 189, 151, 96, 41, 109, 226, 231, 251, 42, 204, 154,
                35, 52, 8, 58, 252, 189, 51, 41, 4, 29, 30, 31, 212, 86
            ]
        );
    }

    #[test]
    fn test_verifier() {
        let evaluations = vec![1; 1 << 20];
        let commitment = commit(&evaluations);
        let evaluation_point = vec![1; 23];
        let proof = prove(&commitment, &evaluations, &evaluation_point);
        assert!(verifier(&commitment, &proof, &evaluation_point));
    }
}
