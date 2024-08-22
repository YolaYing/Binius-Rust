//! This module contains utility functions used by the prover and verifier.
//!
//! The functions in this module are used to perform various operations on the evaluations and the evaluations tensor product.
//! The operations include packing the evaluations into rows, extending the rows, computing the t_prime, and computing the evaluation.
//! In detail, the functions in this module are:
//! 1. choose_row_length_and_count: Choose the row length and row count based on the log of the evaluation count.
//! 2. pack_rows: Pack the evaluations into rows.
//! 3. extend_rows: Extend the rows using the Fast-Fourier extension.
//! 4. evaluation_tensor_product: Compute the tensor product of the evaluations.
//! 5. xor_along_axis: Perform XOR along rows or columns.
//! 6. transpose_bits: Transpose the matrix in the bit-level.
//! 7. transpose: Transpose the matrix
//! 8. computed_tprimes: Compute the t_prime.
//! 9. multisubset: Compute the multisubset sum.
//! 10. transpose_3d: Transpose the 3D matrix.

use super::binary_field16::{big_mul, int_to_bigbin, uint16s_to_bits};
// not use cache
// use super::binary_ntt::extend;
// use cache
use super::binary_ntt_cache::extend;
use crate::binary_field16::BinaryFieldElement16 as B16;
use crate::binary_ntt::WiEvalCache;
use std::convert::TryFrom;

/** transfrom the evaluations into a specific matrix

transform the evaluations into a matrix with row length = 16 bits

Args:
    evaluations: log(size), size = bits of evaluations

Returns:
    log(row length), log(row count), row length, row count
 */
pub fn choose_row_length_and_count(log_evaluation_count: usize) -> (usize, usize, usize, usize) {
    let log_row_length = (log_evaluation_count + 2) / 2;
    let log_row_count = (log_evaluation_count - 1) / 2;
    let row_length = 1 << log_row_length;
    let row_count = 1 << log_row_count;
    (log_row_length, log_row_count, row_length, row_count)
}

/** row packing

perform packing for each row, packing every 16 bits into a unit16, so each row is a list of uint16s
    and the 16 is controlled by packing_factor, and to make later calculation easier, we use BinaryFieldElement16s to represent the unit16s

Args:
    evaluations: the evaluations
    row_count: number of rows
    row_length: the number of bits in a row
    packing_factor: the number of bits in a unit16, control by the packing_factor

Returns:
    a list of rows, each row is a list of BinaryFieldElement16s
 */
pub fn pack_rows(
    evaluations: &[u8],
    row_count: usize,
    row_length: usize,
    packing_factor: usize,
) -> Vec<Vec<B16>> {
    let mut rows = Vec::with_capacity(row_count);
    let mut packed_row_length = row_length / packing_factor;

    // use B16 to represent the unit16s
    for i in 0..row_count {
        let mut packed_row = Vec::with_capacity(packed_row_length);

        for j in 0..packed_row_length {
            // let flipped: Vec<u8>= evaluations[i * row_length /8+ j * packing_factor/8..i * row_length/8 +(j + 1) * packing_factor/8].iter().map(|&byte|byte.reverse_bits()).collect();
            // packed_row.push(B16::new(u16::from_le_bytes(flipped.try_into().unwrap())));
            packed_row.push(B16::new(u16::from_le_bytes(
                evaluations[i * row_length / 8 + j * packing_factor / 8
                    ..i * row_length / 8 + (j + 1) * packing_factor / 8]
                    .try_into()
                    .unwrap(),
            )));
        }
        rows.push(packed_row);
    }
    rows
}

// similar logic as above, but return type is Vec<B16> instead of Vec<Vec<B16>>
// and the inputs are all Vec<u8>
pub fn pack_row(evaluations: &[u8], row_length: usize, packing_factor: usize) -> Vec<B16> {
    let mut packed_row = Vec::with_capacity(row_length / packing_factor);
    for j in 0..row_length / packing_factor {
        let flipped: Vec<u8> = evaluations[j * packing_factor / 8..(j + 1) * packing_factor / 8]
            .iter()
            .map(|&byte| byte.reverse_bits())
            .collect();
        packed_row.push(B16::new(u16::from_le_bytes(flipped.try_into().unwrap())));
    }
    packed_row
}

/** Fast-Fourier extend the rows

Reed-Solomon extension, using the binary-FFT algorithms to extend the rows

Args:
    rows: the packed rows, each row is a list of uint16s
    expansion_factor: EXPANSION_FACTOR, after extension, the row length will be row_length * EXPANSION_FACTOR

Returns:
    the extended rows, each row is a list of uint16s

 */
// Optimized implementation, rows use reference to avoid use row.to_vec(), save 0.75% running time
pub fn extend_rows(rows: &Vec<Vec<B16>>, expansion_factor: usize) -> Vec<Vec<B16>> {
    // use extend function from binary_ntt.rs to extend each row and get the extended rows
    rows.iter()
        .map(|row| extend(row, expansion_factor))
        .collect()
}

/** calculate the tensor product of evaluations

all possible results of walking through pt and at each step taking either coord or 1-coord

Args:
    evaluation_point: the evaluation point, a list of uint128s

Returns:
    field element: the result of the tensor product, a 2^k-long vector of Vector(u16)

 */
// Original implementation
// pub fn evaluation_tensor_product(eval_point: &Vec<u128>) -> Vec<Vec<u16>> {
//     // int_to_bigbin's return type is Vec<u16>
//     let mut o = vec![int_to_bigbin(1)];

//     for coord in eval_point {
//         // for all the element in o, conduct big_mul for the element and &int_to_bigbin(coord)
//         let o_times_coord: Vec<Vec<u16>> = o
//             .iter()
//             .map(|x| big_mul(x, &int_to_bigbin(*coord)))
//             .collect();
//         o = o
//             .iter()
//             .zip(o_times_coord.iter())
//             .map(|(x, y)| x.iter().zip(y.iter()).map(|(a, b)| a ^ b).collect())
//             .collect();
//         o.extend_from_slice(&o_times_coord);
//     }
//     o
// }
// Optimized implementation：save 1.66% prover time, 2.27% verifier time
pub fn evaluation_tensor_product(eval_point: &Vec<u128>) -> Vec<Vec<u16>> {
    // int_to_bigbin's return type is Vec<u16>
    let mut o = vec![int_to_bigbin(1)];

    for coord in eval_point {
        let int_bin = int_to_bigbin(*coord);
        // optimization trick: pre-allocate the o_times_coord vector
        let mut o_times_coord = Vec::with_capacity(o.len());

        // for all the element in o, conduct big_mul for the element and &int_to_bigbin(coord)
        // optimization trick: avoid using big_mul(x, &int_bin) directly in the map closure
        //      and use &o to avoid use o.to_vec() to avoid unnecessary memory allocation
        for x in &o {
            o_times_coord.push(big_mul(x, &int_bin));
        }

        let mut new_o = Vec::with_capacity(o.len() * 2);
        for (x, y) in o.iter().zip(o_times_coord.iter()) {
            let mut combined = Vec::with_capacity(x.len());
            for (a, b) in x.iter().zip(y.iter()) {
                combined.push(a ^ b);
            }
            new_o.push(combined);
        }
        new_o.extend(o_times_coord);
        o = new_o;
    }
    o
}

/** XOR along axis

XOR along rows or columns, if axis = 0, then XOR along rows, if axis = 1, then XOR along columns
    for example:
        we have a matrix like [[1, 2, 3], [4, 5, 6]]
        if axis = 0, then the result will be [1^4, 2^5, 3^6] = [5, 7, 5]
        if axis = 1, then the result will be [1^2^3, 4^5^6] = [0, 7]

Args:
    values: the matrix
    axis: the axis, 0 or 1

Returns:
    the result of XOR along the axis
*/
// Original implementation(takes 15% running time)
// pub fn xor_along_axis(values: &[Vec<u16>], axis: usize) -> Vec<u16> {
//     let mut result: Vec<u16>;

//     if axis == 0 {
//         // XOR along rows (axis=0)
//         result = values[0].clone();
//         for row in values.iter().skip(1) {
//             for (res, val) in result.iter_mut().zip(row.iter()) {
//                 *res ^= val;
//             }
//         }
//     } else if axis == 1 {
//         // XOR along columns (axis=1)
//         result = values.iter().map(|row| row[0]).collect();
//         for i in 1..values[0].len() {
//             for (j, val) in result.iter_mut().enumerate() {
//                 *val ^= values[j][i];
//             }
//         }
//     } else {
//         panic!("Unsupported axis");
//     }

//     result
// }

// Optimized implementation(takes 0.5% running time)
pub fn xor_along_axis(values: &[Vec<u16>], axis: usize) -> Vec<u16> {
    let (rows, _cols) = (values.len(), values[0].len());
    // optimization trick: pre-allocate the result vector
    let mut result = vec![0u16; rows];

    match axis {
        0 => {
            // XOR along rows (axis=0)
            result = values[0].clone();
            for row in &values[1..] {
                for (res, &val) in result.iter_mut().zip(row.iter()) {
                    *res ^= val;
                }
            }
        }
        1 => {
            // XOR along columns (axis=1)
            // optimized trick: cache friendly iteration, iterate over rows first
            for row in 0..rows {
                for (_col, val) in values[row].iter().enumerate() {
                    result[row] ^= val;
                }
            }
        }
        _ => panic!("Unsupported axis"),
    }

    result
}

fn xor_along_axis_4d(values: &Vec<Vec<Vec<Vec<u16>>>>, axis: usize) -> Vec<Vec<Vec<u16>>> {
    let mut result: Vec<Vec<Vec<u16>>> = Vec::new();
    if axis == 0 {
        for i in 0..values[0].len() {
            let mut row = Vec::new();
            for j in 0..values[0][0].len() {
                let mut col = Vec::new();
                for k in 0..values[0][0][0].len() {
                    let mut res = values[0][i][j][k];
                    for l in 1..values.len() {
                        res ^= values[l][i][j][k];
                    }
                    col.push(res);
                }
                row.push(col);
            }
            result.push(row);
        }
    } else if axis == 1 {
        for i in 0..values.len() {
            let mut row = Vec::new();
            for j in 0..values[0][0].len() {
                let mut col = Vec::new();
                for k in 0..values[0][0][0].len() {
                    let mut res = values[i][0][j][k];
                    for l in 1..values[0].len() {
                        res ^= values[i][l][j][k];
                    }
                    col.push(res);
                }
                row.push(col);
            }
            result.push(row);
        }
    } else if axis == 2 {
        for i in 0..values.len() {
            let mut row = Vec::new();
            for j in 0..values[0].len() {
                let mut col = Vec::new();
                for k in 0..values[0][0][0].len() {
                    let mut res = values[i][j][0][k];
                    for l in 1..values[0][0].len() {
                        res ^= values[i][j][l][k];
                    }
                    col.push(res);
                }
                row.push(col);
            }
            result.push(row);
        }
    } else if axis == 3 {
        for i in 0..values.len() {
            let mut row = Vec::new();
            for j in 0..values[0].len() {
                let mut col = Vec::new();
                for k in 0..values[0][0].len() {
                    let mut res = values[i][j][k][0];
                    for l in 1..values[0][0][0].len() {
                        res ^= values[i][j][k][l];
                    }
                    col.push(res);
                }
                row.push(col);
            }
            result.push(row);
        }
    } else {
        panic!("Unsupported axis");
    }
    result
}

/** transpose the bits

ragarding the input as bits, transpose the bits

Args:
    input: the input, a list of list of u8, representing the bits

Returns:
    the output, a list of list of u8, representing the transposed bits
 */
pub fn transpose_bits(input: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut output = vec![vec![0u8; (input.len() + 7) / 8]; input[0].len()];
    for i in 0..input.len() {
        for j in 0..input[0].len() {
            //
            output[j][i / 8] |= (input[i][j] as u8) << ((input.len() - 1 - i) % 8);
        }
    }
    output
}
/** transpose the matrix

different from the transpose_bits, this function transpose the matrix

Args:
    input: the input, a list of list of B16

Returns:
    the output, a transposed list of list of B16
*/
pub fn transpose(input: &Vec<Vec<B16>>) -> Vec<Vec<B16>> {
    let mut output = vec![vec![B16::new(0); input.len()]; input[0].len()];
    for i in 0..input.len() {
        for j in 0..input[0].len() {
            output[j][i] = input[i][j];
        }
    }
    output
}

/** compute the t'


Args:
    rows_as_bits_transpose: the transposed bits, assume rows's shape is (m,n unit16s), then the shape of rows_as_bits_transpose is (n * 16 bits, m)
    row_combination: the row combination, a 2^k-long vector of Vector(u16), len(row_combination) = len(rows) = m, shape of row_combination is (m, 2^k)

Returns:
    the t', a 2D array, shape is (m, 2^k)
*/
// Original implementation
// pub fn computed_tprimes(
//     rows_as_bits_transpose: &Vec<Vec<u8>>,
//     row_combination: &Vec<Vec<u16>>,
// ) -> Vec<Vec<u16>> {
//     // create a 2D array, shape is (m, 2^k), and all the elements are 0
//     let mut t_prime = vec![vec![0u16; row_combination[0].len()]; rows_as_bits_transpose.len()];

//     for j in 0..row_combination[0].len() {
//         // each column of row_combination as comb, and use comb to multiply the bits of each row in rows_as_bits_transpose
//         let multi_res: Vec<Vec<u16>> = rows_as_bits_transpose
//             .iter()
//             .map(|row| {
//                 (0..row_combination.len())
//                     .map(|i| {
//                         // get the i-th bit of the row
//                         let bit = (row[i / 8] >> (7 - i % 8)) & 1;
//                         // multiply the bit with the i-th row and j-th column of row_combination
//                         bit as u16 * row_combination[i][j]
//                     })
//                     .collect()
//             })
//             .collect();

//         // XOR along axis 1
//         let xor_res = xor_along_axis(&multi_res, 1);

//         for (i, res) in xor_res.iter().enumerate() {
//             t_prime[i][j] ^= res;
//         }
//     }

//     t_prime
// }

// Optimized implementation: save 5% prover time, 4% verifier time
pub fn computed_tprimes(
    rows_as_bits_transpose: &Vec<Vec<u8>>,
    row_combination: &Vec<Vec<u16>>,
) -> Vec<Vec<u16>> {
    let m = rows_as_bits_transpose.len();
    let num_bits = rows_as_bits_transpose[0].len() * 8;
    let k = row_combination[0].len();

    // optimization trick: pre-allocate the t_prime vector
    let mut t_prime = vec![vec![0u16; k]; m];
    // optimization trick: pre-allocate the multi_res vector
    let mut multi_res = vec![vec![0u16; num_bits]; m];

    // for each column of row_combination as comb, so we use j to iterate the columns
    for j in 0..k {
        // for each row in rows_as_bits_transpose, so we use i to iterate the rows
        for i in 0..m {
            for bit_pos in 0..num_bits {
                let byte_index = bit_pos / 8;
                let bit_index = 7 - (bit_pos % 8);
                let bit = (rows_as_bits_transpose[i][byte_index] >> bit_index) & 1;
                multi_res[i][bit_pos] = bit as u16 * row_combination[bit_pos][j];
            }
        }

        let xor_res = xor_along_axis(&multi_res, 1);

        for (i, res) in xor_res.iter().enumerate() {
            t_prime[i][j] ^= res;
        }
    }

    t_prime
}

/** transpose the 3D matrix

similar to np.transpose(column_bits, (0,2,1)) in python,
    swaps the rows and columns of each 2D array within the 3D array, effectively turning rows into columns and vice versa.
 */
pub fn transpose_3d(matrix: &Vec<Vec<Vec<u8>>>, order: (usize, usize, usize)) -> Vec<Vec<Vec<u8>>> {
    let dim0 = matrix.len(); // Number of matrices
    let dim1 = matrix[0].len(); // Number of rows in each matrix
    let dim2 = matrix[0][0].len(); // Number of columns in each matrix

    let (new_dim0, new_dim1, new_dim2) = match order {
        (0, 2, 1) => (dim0, dim2, dim1),
        (1, 2, 0) => (dim1, dim2, dim0),
        _ => panic!("Unsupported transpose order"),
    };

    // Initialize the transposed 3D matrix with zeros
    let mut transposed = vec![vec![vec![0; new_dim2]; new_dim1]; new_dim0];

    for i in 0..dim0 {
        for j in 0..dim1 {
            for k in 0..dim2 {
                match order {
                    (0, 2, 1) => transposed[i][k][j] = matrix[i][j][k],
                    (1, 2, 0) => transposed[j][k][i] = matrix[i][j][k],
                    _ => unreachable!(),
                }
            }
        }
    }

    transposed
}

/** Mutisubset sum

Given a list of N objects, and a list of length-N bitvectors representing subsets of those objects,
    compute the xor-sum of each subset. Uses the main subroutine of Pippenger-style algorithms, see: https://ethresear.ch/t/7238

Args:
    values: the values(row_combination, Vec<Vec<u16>)
    bits: the bits(transposed_column_bits, Vec<Vec<Vec<u8>>)
*/
pub fn multisubset(values: &Vec<Vec<u16>>, bits: &Vec<Vec<Vec<u8>>>) -> Vec<Vec<Vec<u16>>> {
    let GROUPING = 4;
    let mut subsets = vec![vec![vec![0u16; values[0].len()]; 16]; values.len() / GROUPING];

    for i in 0..GROUPING {
        for j in (0..values.len()).step_by(GROUPING) {
            subsets[j / GROUPING][1 << i] = values[j + i].clone();
        }
    }

    // generate the subsets
    let mut top_p_of_2 = 2;
    for i in 3..1 << GROUPING {
        if (i & (i - 1)) == 0 {
            top_p_of_2 = i;
        } else {
            for j in (0..values.len()).step_by(GROUPING) {
                for k in 0..values[0].len() {
                    subsets[j / GROUPING][i][k] = subsets[j / GROUPING][top_p_of_2][k]
                        ^ subsets[j / GROUPING][i - top_p_of_2][k];
                }
            }
        }
    }

    // use bits to generate the index_columns, and then use the index_columns to select the elements from subsets
    let index_columns: Vec<Vec<Vec<u8>>> = bits
        .iter()
        .map(|matrix| {
            matrix
                .iter()
                .map(|row| {
                    row.chunks(4)
                        .map(|chunk| chunk.iter().rev().fold(0, |acc, &bit| (acc << 1) | bit))
                        .collect()
                })
                .collect()
        })
        .collect();

    // use the index_columns to select the elements from subsets
    let selected_elements: Vec<Vec<Vec<Vec<u16>>>> = index_columns
        .iter()
        .map(|outer| {
            outer
                .iter()
                .map(|inner| {
                    inner
                        .iter()
                        .enumerate()
                        .map(|(i, &index)| subsets[i][index as usize].clone())
                        .collect()
                })
                .collect()
        })
        .collect();

    // XOR along axis 3
    let o = xor_along_axis_4d(&selected_elements, 2);
    o
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;

    #[test]
    fn test_choose_row_length_and_count() {
        let (log_row_length, log_row_count, row_length, row_count) = choose_row_length_and_count(6);
        assert_eq!(log_row_length, 4);
        assert_eq!(log_row_count, 2);
        assert_eq!(row_length, 16);
        assert_eq!(row_count, 4);
    }

    #[test]
    fn test_extend() {
        let rows = vec![
            vec![B16::new(1), B16::new(3)],
            vec![B16::new(9), B16::new(15)],
        ];
        let extended_rows = extend_rows(&rows, 2);
        assert_eq!(
            extended_rows[0],
            vec![B16::new(1), B16::new(3), B16::new(2), B16::new(0)]
        );
        assert_eq!(
            extended_rows[1],
            vec![B16::new(9), B16::new(15), B16::new(2), B16::new(4)]
        );
    }

    #[test]
    fn test_evaluation_tensor_product() {
        let eval_point = vec![2, 5];
        let result = evaluation_tensor_product(&eval_point);
        // i need to compare the equality of two vectors of Vec(Vec(u16)), the following code is not working
        // assert_eq!(result, vec![int_to_bigbin(12), int_to_bigbin(8), int_to_bigbin(15), int_to_bigbin(10)]);
        assert_eq!(result[0], int_to_bigbin(12));
        assert_eq!(result[1], int_to_bigbin(8));
        assert_eq!(result[2], int_to_bigbin(15));
        assert_eq!(result[3], int_to_bigbin(10));
    }

    #[test]
    fn test_xor_along_axis() {
        let values = vec![vec![1, 2, 3], vec![4, 5, 6]];
        let result = xor_along_axis(&values, 0);
        assert_eq!(result, vec![5, 7, 5]);
        let result = xor_along_axis(&values, 1);
        assert_eq!(result, vec![0, 7]);
    }

    #[test]
    fn test_transpose_bits() {
        let data = vec![
            vec![B16::new(1), B16::new(3)],
            vec![B16::new(9), B16::new(15)],
        ];
        // use uint16s_to_bits to convert the data into bits row by row
        let input = data.iter().map(|row| uint16s_to_bits(row)).collect();
        let output = transpose_bits(input);
        assert_eq!(output[0], [3]);
    }

    #[test]
    fn test_transpose() {
        let data = vec![
            vec![B16::new(1), B16::new(3)],
            vec![B16::new(9), B16::new(15)],
        ];
        let output = transpose(&data);
        assert_eq!(output[0], [B16::new(1), B16::new(9)]);
        assert_eq!(output[1], [B16::new(3), B16::new(15)]);
    }

    #[test]
    fn test_computed_tprimes() {
        let eval_point = vec![2, 5];
        let rows = vec![
            vec![B16::new(1), B16::new(3)],
            vec![B16::new(9), B16::new(15)],
            vec![B16::new(2), B16::new(4)],
            vec![B16::new(0), B16::new(0)],
        ];

        let rows_as_bits_transpose =
            transpose_bits(rows.iter().map(|row| uint16s_to_bits(row)).collect());
        let row_combination = evaluation_tensor_product(&eval_point);
        let result = computed_tprimes(&rows_as_bits_transpose, &row_combination);

        assert_eq!(result[0], [4, 0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_pack_row() {
        // data =  [1 1 0 1 0 0 0 0 0 0 1 0 1 0 0 0]
        let data = vec![0b11010000, 0b00101000];
        let result = pack_row(&data, 16, 16);
        // check if =  [5131]
        assert_eq!(result, [B16::new(5131)]);
    }

    #[test]
    fn test_pack_rows() {
        let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
        let result = pack_rows(&data, 8, 16, 16);
        assert_eq!(result[0], [B16::new(513)]);
        assert_eq!(result[1], [B16::new(1027)]);
        assert_eq!(result[2], [B16::new(1541)]);
    }
}
