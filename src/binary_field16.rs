//! This module defines the `BinaryFieldElement16` struct and its methods.
//! The `BinaryFieldElement16` struct is a wrapper around `u16` and is used in the binary field.
//! The binary field consists of two elements, 0 and 1, with operations defined as XOR, AND, and OR.
//!
//! The `BinaryFieldElement16` struct implements the following traits:
//! - `Add`, `Sub`, `Mul`, `Div`, and `Neg` for arithmetic operations.
//! - `FromIterator` to convert a vector of `BinaryFieldElement16` into a vector of `u8`.
//! - `BigMul` to multiply two large binary numbers.
//! - `ToU16` to convert a `BinaryFieldElement16` into a `u16`.
//!
//! Additionally, the `BinaryFieldElement16` struct provides the following functions:
//! - `int_to_bigbin`: Converts a 128-bit integer into a length-8 vector of `u16`.
//! - `uint16s_to_bits`: Converts a vector of `u16` into bits.
//! - `uint16_to_bit`: Converts a `BinaryFieldElement16` into bits.
//! - `bin_mul`: Multiplies two binary numbers in the binary tower field.
//! - `big_mul`: Multiplies two large binary numbers.
//! - `mul_by_Xi`: Multiplies a large binary number by `Xi`.

use serde::{Deserialize, Serialize};
use std::ops::{Add, Div, Mul, Neg, Sub};

/**
A binary field element：a wrapper of u64
 */
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct BinaryFieldElement16 {
    pub value: u16,
}

impl BinaryFieldElement16 {
    pub fn new(value: u16) -> Self {
        BinaryFieldElement16 { value }
    }

    /** Get the bit length of the element

    find the smallest power of 2 that is greater than the element, and count the number zeros before the first 1
        then return 2^count

    Returns:
        usize: the bit length of the element
     */
    fn bit_length(&self) -> u16 {
        16 - self.value.leading_zeros() as u16
    }

    /** Get the inverse of the element

    inverse = element^(2^(bit_length(element) - 2))

    Returns:
        BinaryFieldElement: the inverse of the element

    */
    pub fn inv(&self) -> Self {
        // L = 1 << (self.value.bit_length() - 1).bit_length()
        // return self ** (2**L - 2)
        let l = 1 << (16 - (self.bit_length() - 1).leading_zeros());
        self.pow(2u16.pow(l as u32) - 2)
    }

    /** Get the power of the element

    power = element^(exp), and it is calculated recursively, using the following rules:
        1. if exp = 0, return 1
        2. if exp = 1, return element
        3. if exp = 2, return element * element
        4. if exp is even, return (element^(exp/2))^2
        5. if exp is odd, return element * (element^(exp - 1))

    Args:
        exp (u16): the exponent, important: exp is not binary field element, it is u16

     */
    fn pow(&self, exp: u16) -> Self {
        if exp == 0 {
            BinaryFieldElement16::new(1)
        } else if exp == 1 {
            *self
        } else if exp == 2 {
            *self * *self
        } else {
            self.pow(exp % 2) * self.pow(exp / 2).pow(2)
        }
    }
}

/** Implement the Add trait for BinaryFieldElement

   The addition of two binary field elements is the XOR of the two elements

   Args:
       other (BinaryFieldElement): the other element to add

   Returns:
       BinaryFieldElement: the sum of the two elements

*/
impl Add for BinaryFieldElement16 {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        // let result = if USE_CACHE && self.value < 256 && other.value < 256 {
        //     unsafe { ADDCACHE[self.value as usize][other.value as usize].unwrap() }
        // } else {
        //     BinaryFieldElement::new(self.value ^ other.value)
        // };
        let result = BinaryFieldElement16::new(self.value ^ other.value);
        result
    }
}

/** Implement the Sub trait for BinaryFieldElement

   The subtraction of two binary field elements is the same as the addition

   Args:
       other (BinaryFieldElement): the other element to subtract

   Returns:
       BinaryFieldElement: the difference of the two elements

*/
impl Sub for BinaryFieldElement16 {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        self + other
    }
}

/** Implement the Neg trait for BinaryFieldElement

   The negation of a binary field element is the element itself

   Returns:
       BinaryFieldElement: the negation of the element

*/
impl Neg for BinaryFieldElement16 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self
    }
}

/** Implement the Mul trait for BinaryFieldElement

   The multiplication of two binary field elements is calculated using the Karatsuba algorithm(implemented in binmul)

   Args:
       other (BinaryFieldElement): the other element to multiply

   Returns:
       BinaryFieldElement: the product of the two elements

*/
impl Mul for BinaryFieldElement16 {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        // let result = if USE_CACHE && self.value < 256 && other.value < 256 {
        //     unsafe { MULCACHE[self.value as usize][other.value as usize].unwrap() }
        // } else {
        //     BinaryFieldElement::new(binmul(self.value, other.value, None))
        // };
        let result = BinaryFieldElement16::new(bin_mul(self.value, other.value, None));
        result
    }
}

/** Implement the Div trait for BinaryFieldElement

   The division of two binary field elements is the multiplication of the first element and the inverse of the second element

   Args:
       other (BinaryFieldElement): the other element to divide

   Returns:
       BinaryFieldElement: the quotient of the two elements

*/
impl Div for BinaryFieldElement16 {
    type Output = Self;

    fn div(self, other: Self) -> Self::Output {
        self * other.inv()
    }
}
/** used in backed_colunms step

Convert a vector of BinaryFieldElement16 into a vector of u8

Args:
    data: the vector of BinaryFieldElement16

Returns:
    Vec<u8>: the vector of u8
*/

impl FromIterator<BinaryFieldElement16> for Vec<u8> {
    fn from_iter<I: IntoIterator<Item = BinaryFieldElement16>>(iter: I) -> Self {
        let mut vec = Vec::new();
        for element in iter {
            let high_byte = (element.value >> 8) as u8;
            let low_byte = (element.value & 0xFF) as u8;
            vec.push(low_byte);
            vec.push(high_byte);
        }
        vec
    }
}

/** Multiply v1 * v2 in the binary tower field

   The multiplication of two binary field elements is calculated using the Karatsuba algorithm

   Args:
       v1 (u16): the first element, important: v1 is not binary field element, it is u16
       v2 (u16): the second element, important: v2 is not binary field element, it is u16
       length (Option<usize>): the length of the elements

   Returns:
       u16: the product of the two elements

   Appendix:
   See https://blog.lambdaclass.com/snarks-on-binary-fields-binius/ for introduction to how binary tower fields work
*/
pub fn bin_mul(v1: u16, v2: u16, length: Option<usize>) -> u16 {
    // if USE_CACHE && v1 < 256 && v2 < 256 && unsafe { RAWMULCACHE[v1 as usize][v2 as usize].is_some() } {
    //     return unsafe { RAWMULCACHE[v1 as usize][v2 as usize].unwrap() };
    // }
    if v1 < 2 || v2 < 2 {
        return v1 * v2;
    }

    let length = match length {
        Some(l) => l,
        None => {
            // get the maxium of v1 and v2, and find the smallest power of 2 that is greater than the maxium
            let max_v = v1.max(v2);
            // find the number of zeros before the max_v
            // leading_zeros() returns the number of leading zeros in the binary representation of the number, the overall length of the binary representation is depends on the input number
            // for example, 8 as u16 is 1000 in binary, so the leading_zeros() is 12
            // the return value of leading_zeros() is type u32
            let bit_length = 16 - max_v.leading_zeros();
            // the type of bit_length is u32, so the overall length of the binary representation is 32
            let adjusted_bit_length = 32 - (bit_length - 1).leading_zeros();
            // the length of the elements is 2^adjusted_bit_length
            1 << adjusted_bit_length
        }
    };

    let halflen = length / 2;
    let quarterlen = length / 4;
    let halfmask = (1 << halflen) - 1;

    let (l1, r1) = (v1 & halfmask, v1 >> halflen);
    let (l2, r2) = (v2 & halfmask, v2 >> halflen);

    // # Optimized special case (used to compute R1R2_high), sec III of
    // https://ieeexplore.ieee.org/document/612935
    if (l1, r1) == (0, 1) {
        let out_r = bin_mul(1 << quarterlen, r2, Some(halflen)) ^ l2;
        return r2 ^ (out_r << halflen);
    }

    // x_{i+1}^2 reduces to 1 + x_{i+1} * x_i
    // Uses Karatsuba to only require three sub-multiplications for each input
    // halving (R1R2_high doesn't count because of the above optimization)
    let l1l2 = bin_mul(l1, l2, Some(halflen));
    let r1r2 = bin_mul(r1, r2, Some(halflen));
    let r1r2_high = bin_mul(1 << quarterlen, r1r2, Some(halflen));
    let z3 = bin_mul(l1 ^ r1, l2 ^ r2, Some(halflen));
    l1l2 ^ r1r2 ^ ((z3 ^ l1l2 ^ r1r2 ^ r1r2_high) << halflen)
}

/** Multiplies together two list of binary number, using the Karatsuba algorithm

different from the function in binary_field.rs, this function is used to compute two big binary numbers
    and the type of each binary number is Vec<u16>

Args:
    x1: the first big binary number, the type is Vec<u16>
    x2: the second big binary number, the type is Vec<u16>

Returns:
    Vec<u16>: the product of the two big binary numbers
 */
trait BigMul {
    fn big_mul(self, other: Self) -> Vec<u16>;
}

impl BigMul for Vec<u16> {
    fn big_mul(self, other: Vec<u16>) -> Vec<u16> {
        big_mul_impl(&self, &other)
    }
}

impl<'a> BigMul for &'a Vec<u16> {
    fn big_mul(self, other: &'a Vec<u16>) -> Vec<u16> {
        big_mul_impl(self, other)
    }
}

// fn big_mul_impl(x1: &Vec<u16>, x2: &Vec<u16>) -> Vec<u16> {
//     let n = x1.len();
//     if n == 1 {
//         return vec![bin_mul(x1[0], x2[0], None)];
//     }
//     let l1 = x1[..n / 2].to_vec();
//     let l2 = x2[..n / 2].to_vec();
//     let r1 = x1[n / 2..].to_vec();
//     let r2 = x2[n / 2..].to_vec();
//     let l1l2 = big_mul_impl(&l1, &l2);
//     let r1r2 = big_mul_impl(&r1, &r2);
//     let r1r2_high = mul_by_Xi(r1r2.clone(), n / 2);
//     let z3 = big_mul_impl(
//         &l1.iter().zip(r1.iter()).map(|(a, b)| a ^ b).collect(),
//         &l2.iter().zip(r2.iter()).map(|(a, b)| a ^ b).collect(),
//     );
//     let part1 = l1l2
//         .iter()
//         .zip(r1r2.iter())
//         .map(|(a, b)| a ^ b)
//         .collect::<Vec<u16>>();
//     let part2 = z3
//         .iter()
//         .zip(l1l2.iter())
//         .zip(r1r2.iter())
//         .zip(r1r2_high.iter())
//         .map(|(((a, b), c), d)| a ^ b ^ c ^ d)
//         .collect::<Vec<u16>>();
//     let mut result = Vec::new();
//     result.extend_from_slice(&part1);
//     result.extend_from_slice(&part2);
//     result
// }
fn big_mul_impl(x1: &Vec<u16>, x2: &Vec<u16>) -> Vec<u16> {
    let n = x1.len();
    if n == 1 {
        return vec![bin_mul(x1[0], x2[0], None)];
    }

    let (l1, r1) = x1.split_at(n / 2);
    let (l2, r2) = x2.split_at(n / 2);

    let l1l2 = big_mul_impl(&l1.to_vec(), &l2.to_vec());
    let r1r2 = big_mul_impl(&r1.to_vec(), &r2.to_vec());
    let r1r2_high = mul_by_Xi(&r1r2.clone(), n / 2);

    let z3 = big_mul_impl(
        &l1.iter().zip(r1.iter()).map(|(a, b)| a ^ b).collect(),
        &l2.iter().zip(r2.iter()).map(|(a, b)| a ^ b).collect(),
    );

    let part1 = l1l2
        .iter()
        .zip(r1r2.iter())
        .map(|(a, b)| a ^ b)
        .collect::<Vec<u16>>();

    let part2 = z3
        .iter()
        .zip(l1l2.iter())
        .zip(r1r2.iter())
        .zip(r1r2_high.iter())
        .map(|(((a, b), c), d)| a ^ b ^ c ^ d)
        .collect::<Vec<u16>>();

    let mut result = Vec::with_capacity(part1.len() + part2.len());
    result.extend_from_slice(&part1);
    result.extend_from_slice(&part2);
    result
}

pub fn big_mul<T: BigMul>(x1: T, x2: T) -> Vec<u16> {
    x1.big_mul(x2)
}

/** Multiply a big binary number by Xi

multiply the right half of the big binary number by Xi, and then XOR the result with the left half

Args:
    x: the big binary number, the type is Vec<u16>
    n: the length of the big binary number

Returns:
    Vec<u16>: the product of the big binary number and Xi
 */
// pub fn mul_by_Xi(x: &Vec<u16>, n: usize) -> Vec<u16> {
//     if x.len() == 1 {
//         return vec![bin_mul(x[0], 256, None)];
//     }
//     let l = x[..n / 2].to_vec();
//     let r = x[n / 2..].to_vec();
//     let out_r = mul_by_Xi(&r, n / 2)
//         .iter()
//         .zip(l.iter())
//         .map(|(a, b)| a ^ b)
//         .collect::<Vec<u16>>();
//     let mut result = Vec::new();
//     result.extend_from_slice(&r);
//     result.extend_from_slice(&out_r);
//     result
// }
pub fn mul_by_Xi(x: &Vec<u16>, n: usize) -> Vec<u16> {
    if x.len() == 1 {
        return vec![bin_mul(x[0], 256, None)];
    }

    let (l, r) = x.split_at(n / 2);

    let out_r = mul_by_Xi(&r.to_vec(), n / 2)
        .iter()
        .zip(l.iter())
        .map(|(a, b)| a ^ b)
        .collect::<Vec<u16>>();

    let mut result = Vec::with_capacity(x.len());
    result.extend_from_slice(r);
    result.extend_from_slice(&out_r);
    result
}

/** Convert a 128-bit integer into a length-8 vector of uint16's

right shift the integer by 16 bits each time, and take the last 16 bits as the uint16

Args:
    value: the 128-bit integer

Returns:
    field element: a length-8 vector of uint16's

 */
pub fn int_to_bigbin(x: u128) -> Vec<u16> {
    let mut result = Vec::new();
    for k in 0..8 {
        result.push(((x >> (k * 16)) & 65535) as u16);
    }
    result
}

/** Convert a vector of uint16's into bits

right shift the uint16 by 1 bit each time, and take the last bit as the bit

Args:
    data: the vector of uint16's

Returns:
    Vec<u8>: the bits
*/
pub trait ToU16 {
    fn to_u16(&self) -> u16;
}

impl ToU16 for u16 {
    fn to_u16(&self) -> u16 {
        *self
    }
}

impl ToU16 for BinaryFieldElement16 {
    fn to_u16(&self) -> u16 {
        self.value
    }
}
pub fn uint16s_to_bits<T: ToU16>(data: &Vec<T>) -> Vec<u8> {
    let mut result = Vec::with_capacity(data.len() * 16);

    for value in data {
        // Extract each bit from the 16-bit value
        let value_u16 = value.to_u16();
        for i in 0..16 {
            result.push(((value_u16 >> i) & 1) as u8);
        }
    }
    result
}

pub fn uint16_to_bit(value: &BinaryFieldElement16) -> Vec<u8> {
    let mut result = Vec::with_capacity(16);
    for i in 0..16 {
        result.push(((value.value >> i) & 1) as u8);
    }
    result
}

/** Implement the Serialize trait for BinaryFieldElement

Serialize the element as a string

Args:
    serializer (S): the serializer

Returns:
    S::Ok: the serialized element

*/
impl Serialize for BinaryFieldElement16 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&format!("{:X}", self.value))
    }
}

/** Implement the Deserialize trait for BinaryFieldElement

Deserialize the element from a string

Args:
    deserializer (D): the deserializer

Returns:
    Result<BinaryFieldElement, D::Error>: the deserialized element

*/
impl<'de> Deserialize<'de> for BinaryFieldElement16 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        let value = u16::from_str_radix(&s, 16).map_err(serde::de::Error::custom)?;
        Ok(BinaryFieldElement16 { value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bin_mul() {
        assert_eq!(bin_mul(3, 5, None), 15);
        assert_eq!(bin_mul(7, 11, None), 4);
        assert_eq!(bin_mul(8, 2, None), 12);
        assert_eq!(bin_mul(32147, 48725, None), 43100);
    }

    #[test]
    fn test_binary_field_element_add() {
        let a = BinaryFieldElement16::new(8);
        let b = BinaryFieldElement16::new(5);
        assert_eq!(a + b, BinaryFieldElement16::new(13));
    }

    #[test]
    fn test_binary_field_element_sub() {
        let a = BinaryFieldElement16::new(8);
        let b = BinaryFieldElement16::new(5);
        assert_eq!(a - b, BinaryFieldElement16::new(13));
    }

    #[test]
    fn test_binary_field_element_mul() {
        let a = BinaryFieldElement16::new(8);
        let b = BinaryFieldElement16::new(5);
        assert_eq!(a * b, BinaryFieldElement16::new(6));
    }

    #[test]
    fn test_binary_field_element_div() {
        let a = BinaryFieldElement16::new(0);
        let b = BinaryFieldElement16::new(1);
        assert_eq!(a / b, BinaryFieldElement16::new(0));
    }

    #[test]
    fn test_binary_field_element_inv() {
        let a = BinaryFieldElement16::new(1);
        assert_eq!(a.inv(), BinaryFieldElement16::new(1));
    }

    #[test]
    fn test_binary_field_element_pow() {
        let a = BinaryFieldElement16::new(2);
        assert_eq!(a.pow(3), BinaryFieldElement16::new(1));
    }

    #[test]
    fn test_big_mul() {
        // big_mul(int_to_bigbin(3**29), int_to_bigbin(5**29))= [46732 49627 26993 63626 14101 27237 21150     0]
        let a = int_to_bigbin(3u128.pow(29));
        let b = int_to_bigbin(5u128.pow(29));
        let result = big_mul(a, b);
        assert_eq!(
            result,
            vec![46732, 49627, 26993, 63626, 14101, 27237, 21150, 0]
        );
    }

    #[test]
    fn test_uint16s_to_bits() {
        let data = vec![BinaryFieldElement16::new(1u16)];
        let result = uint16s_to_bits(&data);
        assert_eq!(result, vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        // test on [[1,3]],check result != [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        let data = vec![
            BinaryFieldElement16::new(1u16),
            BinaryFieldElement16::new(3u16),
        ];
        let result = uint16s_to_bits(&data);
        assert_eq!(
            result,
            vec![
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );
    }
}
