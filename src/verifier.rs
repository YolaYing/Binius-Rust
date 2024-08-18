const EXPANSION_FACTOR: usize = 4;
const NUM_CHALLENGES: usize = 2;
const PACKING_FACTOR: usize = 16;

use sha2::digest::typenum::uint;

use super::prover::Proof;
use super::utils::{
    choose_row_length_and_count, transpose_bits, extend_rows, 
    evaluation_tensor_product, computed_tprimes, pack_row, 
    xor_along_axis, transpose_3d, multisubset};
use super::challenger::get_challenges; 
use super::merkle_tree::verify_branch;
use super::binary_field16::{uint16s_to_bits,uint16_to_bit, big_mul, BinaryFieldElement16};

fn verify_optimized_binius_proof(proof:Proof) -> bool {

    let columns = proof.columns;
    let evaluation_point = proof.evaluation_point;
    let value = proof.eval;
    let t_prime = proof.t_prime;
    let root = proof.root;
    let branches = proof.branches;

    // Compute the row length and row count of the grid. Should output same numbers as what prover gave
    let (log_row_length, log_row_count, row_length, row_count) = choose_row_length_and_count(evaluation_point.len());
    let extended_row_length = row_length * EXPANSION_FACTOR / PACKING_FACTOR;
    
    // Compute challenges. Should output the same as what prover computed
    let challenges = get_challenges(&root, extended_row_length, NUM_CHALLENGES);

    // Verify Merkle branches
    for i in 0..NUM_CHALLENGES {
        let challenge = challenges[i];
        let packed_column:Vec<u8> = columns[challenge as usize]
            .clone()
            .into_iter()
            .collect();
        let branch = branches[i].clone();
        assert!(verify_branch(&root, challenge as usize, &packed_column, &branch));
    }

    // Use the same Reed-Solomon code that the prover used to extend the rows,
    // but to extend t_prime. We do this separately for each bit of t_prime
    // each row in t_prime is a list of uint16s, use uint16s_to_bits to convert it to a list of bits
    let t_prime_bits = t_prime.iter().map(|row| uint16s_to_bits(row)).collect();
    // transpose the bits
    let t_prime_bits_transpose = transpose_bits(t_prime_bits);
    // pack the each row of t_prime_bits_transpose into a list of BinaryFieldElement16s
    let t_prime_columns: Vec<Vec<BinaryFieldElement16>>= t_prime_bits_transpose.iter().map(|row| pack_row(row, t_prime_bits_transpose[0].len() * 8, PACKING_FACTOR)).collect();
    // extend the rows
    let extended_t_prime_columns = extend_rows(&t_prime_columns, EXPANSION_FACTOR);

    // Here, we take advantage of the linearity of the code. A linear combination of the Reed-Solomon extension gives the same result as an extension of the linear combination.
    let row_combination = evaluation_tensor_product(&evaluation_point[log_row_length..].to_vec());
    // Each column is a vector of row_count uint16's. Convert each uint16 into bits
    let column_bits: Vec<Vec<Vec<u8>>>= columns.iter().map(|col| col.iter().map(|uint16| uint16_to_bit(uint16)).collect()).collect();
    // Take the same linear combination the prover used to compute t_prime, and apply it to the columns of bits.
    let transposed_column_bits = transpose_3d(&column_bits, (0, 2, 1));
    let computed_tprimes = multisubset(&row_combination, &transposed_column_bits);
    // Turn the computed tprimes into bits using uint16s_to_bits
    let computed_tprime_bits: Vec<Vec<Vec<u8>>> = computed_tprimes.iter()
        .map(|row| row.iter()
            .map(|uint16| uint16s_to_bits(uint16))
            .collect())
        .collect();

    // Convert our FFT-extended t_prime rows into bits
    // step 1: use challenge to select columns, and convert to bits
    let extended_t_prime_columns_slices: Vec<Vec<Vec<BinaryFieldElement16>>> = extended_t_prime_columns.iter().map(|row| challenges.iter().map(|&c| vec![row[c as usize]]).collect()).collect();
    let extended_t_prime_bits: Vec<Vec<Vec<u8>>> = extended_t_prime_columns_slices.iter().map(|row| row.iter().map(|uint16| uint16s_to_bits(uint16)).collect()).collect();
    // step 2: transpose the bits
    let extended_t_prime_bits_transpose = transpose_3d(&extended_t_prime_bits, (1, 2, 0)); 
 
    // The bits of the t_prime extension should equal the bits of the row linear combination of the column bits
    assert_eq!(computed_tprime_bits, extended_t_prime_bits_transpose);

    // Compute the evaluation
    let col_combination = evaluation_tensor_product(&evaluation_point[..log_row_length].to_vec());
    let computed_eval = xor_along_axis(
        &t_prime.iter()
            .zip(col_combination.iter())
            .map(|(t_prime_row, col_combination_row)| big_mul(t_prime_row, col_combination_row))
            .collect::<Vec<Vec<u16>>>(),
        0
    );
    assert_eq!(computed_eval, value);
    true
   
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::prover::prover;
    use super::super::verifier::verify_optimized_binius_proof;

    #[test]
    fn test_verify_optimized_binius_proof() {
        let evaluations = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
        let evaluation_point = vec![1, 2, 3, 4, 5, 6, 7];
        let proof = prover(&evaluations, evaluation_point.clone());
        assert!(verify_optimized_binius_proof(proof));
    }
}