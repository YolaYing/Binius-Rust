// This module leverage simd instrcution to accelerate the multiplication of binary field elements.
// This module can be well supported by the CPU with simd instruction set.
// Original Binius further leverages the CLMUL（Carry-less Multiplication）instruction set
// to accelerate the multiplication of binary field elements, which is not supported by Apple M series chips.

use core::{arch::aarch64::*, mem};
use rayon::{result, str};
use serde::{Deserialize, Serialize};
use std::{
    ops::{Add, Div, Mul, Neg, Sub},
    slice::RSplit,
};

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
// trait BigMul {
//     fn big_mul(self, other: Self) -> Vec<u16>;
// }

// impl BigMul for Vec<u16> {
//     fn big_mul(self, other: Vec<u16>) -> Vec<u16> {
//         big_mul_impl(&self, &other)
//     }
// }

// impl<'a> BigMul for &'a Vec<u16> {
//     fn big_mul(self, other: &'a Vec<u16>) -> Vec<u16> {
//         big_mul_impl(self, other)
//     }
// }

// original implementation
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

// optimized implementation
// fn big_mul_impl(x1: &Vec<u16>, x2: &Vec<u16>) -> Vec<u16> {
//     let n = x1.len();
//     if n == 1 {
//         return vec![bin_mul(x1[0], x2[0], None)];
//     }

//     let (l1, r1) = x1.split_at(n / 2);
//     let (l2, r2) = x2.split_at(n / 2);

//     let l1l2 = big_mul_impl(&l1.to_vec(), &l2.to_vec());
//     let r1r2 = big_mul_impl(&r1.to_vec(), &r2.to_vec());
//     let r1r2_high = mul_by_Xi(&r1r2.clone(), n / 2);

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

//     let mut result = Vec::with_capacity(part1.len() + part2.len());
//     result.extend_from_slice(&part1);
//     result.extend_from_slice(&part2);
//     result
// }

// pub fn big_mul<T: BigMul>(x1: T, x2: T) -> Vec<u16> {
//     x1.big_mul(x2)
// }

// Montgomery (modular) multiplication is a method that allows computing such multiplications faster. Instead of dividing the product and subtracting
// $n$  multiple times, it adds multiples of
// $n$  to cancel out the lower bits and then just discards the lower bits.

pub fn big_mul(x1: u128, x2: u128) -> u128 {
    montgomery_multiply(x1, x2)
}

#[inline]
fn montgomery_multiply(a: u128, b: u128) -> u128 {
    unsafe {
        let h = vreinterpretq_u8_p128(a);
        let y = vreinterpretq_u8_p128(b);
        let (h, m, l) = karatsuba1(h, y);
        let (h, l) = karatsuba2(h, m, l);
        vreinterpretq_p128_u8(mont_reduce(h, l))
    }
}

/// Karatsuba decomposition for `x*y`.
#[inline]
unsafe fn karatsuba1(x: uint8x16_t, y: uint8x16_t) -> (uint8x16_t, uint8x16_t, uint8x16_t) {
    // First Karatsuba step: decompose x and y.
    //
    // (x1*y0 + x0*y1) = (x1+x0) * (y1+x0) + (x1*y1) + (x0*y0)
    //        M                                 H         L
    //
    // m = x.hi^x.lo * y.hi^y.lo
    // let z3 = big_mul_impl(
    //     &l1.iter().zip(r1.iter()).map(|(a, b)| a ^ b).collect(),
    //     &l2.iter().zip(r2.iter()).map(|(a, b)| a ^ b).collect(),
    // );
    let m = pmull(
        veorq_u8(x, vextq_u8(x, x, 8)), // x.hi^x.lo
        veorq_u8(y, vextq_u8(y, y, 8)), // y.hi^y.lo
    );
    // let l1l2 = big_mul_impl(&l1.to_vec(), &l2.to_vec());
    let h = pmull2(x, y); // h = x.hi * y.hi
                          // r1r2 = big_mul_impl(&r1.to_vec(), &r2.to_vec());
    let l = pmull(x, y); // l = x.lo * y.lo
    (h, m, l)
}

/// Karatsuba combine.
#[inline]
unsafe fn karatsuba2(h: uint8x16_t, m: uint8x16_t, l: uint8x16_t) -> (uint8x16_t, uint8x16_t) {
    // Second Karatsuba step: combine into a 2n-bit product.
    //
    // m0 ^= l0 ^ h0 // = m0^(l0^h0)
    // m1 ^= l1 ^ h1 // = m1^(l1^h1)
    // l1 ^= m0      // = l1^(m0^l0^h0)
    // h0 ^= l0 ^ m1 // = h0^(l0^m1^l1^h1)
    // h1 ^= l1      // = h1^(l1^m0^l0^h0)
    let t = {
        //   {m0, m1} ^ {l1, h0}
        // = {m0^l1, m1^h0}
        let t0 = veorq_u8(m, vextq_u8(l, h, 8));

        //   {h0, h1} ^ {l0, l1}
        // = {h0^l0, h1^l1}
        let t1 = veorq_u8(h, l);

        //   {m0^l1, m1^h0} ^ {h0^l0, h1^l1}
        // = {m0^l1^h0^l0, m1^h0^h1^l1}
        veorq_u8(t0, t1)
    };

    // {m0^l1^h0^l0, l0}
    let x01 = vextq_u8(
        vextq_u8(l, l, 8), // {l1, l0}
        t,
        8,
    );

    // {h1, m1^h0^h1^l1}
    let x23 = vextq_u8(
        t,
        vextq_u8(h, h, 8), // {h1, h0}
        8,
    );

    (x23, x01)
}

#[inline]
unsafe fn mont_reduce(x23: uint8x16_t, x01: uint8x16_t) -> uint8x16_t {
    // Perform the Montgomery reduction over the 256-bit X.
    //    [A1:A0] = X0 • poly
    //    [B1:B0] = [X0 ⊕ A1 : X1 ⊕ A0]
    //    [C1:C0] = B0 • poly
    //    [D1:D0] = [B0 ⊕ C1 : B1 ⊕ C0]
    // Output: [D1 ⊕ X3 : D0 ⊕ X2]
    let poly = vreinterpretq_u8_p128(1 << 127 | 1 << 126 | 1 << 121 | 1 << 63 | 1 << 62 | 1 << 57);
    let a = pmull(x01, poly);
    let b = veorq_u8(x01, vextq_u8(a, a, 8));
    let c = pmull2(b, poly);
    veorq_u8(x23, veorq_u8(c, b))
}

/// Multiplies the low bits in `a` and `b`.
#[inline]
unsafe fn pmull(a: uint8x16_t, b: uint8x16_t) -> uint8x16_t {
    mem::transmute(vmull_p64(
        vgetq_lane_u64(vreinterpretq_u64_u8(a), 0),
        vgetq_lane_u64(vreinterpretq_u64_u8(b), 0),
    ))
}

/// Multiplies the high bits in `a` and `b`.
#[inline]
unsafe fn pmull2(a: uint8x16_t, b: uint8x16_t) -> uint8x16_t {
    mem::transmute(vmull_p64(
        vgetq_lane_u64(vreinterpretq_u64_u8(a), 1),
        vgetq_lane_u64(vreinterpretq_u64_u8(b), 1),
    ))
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
    let mut result = Vec::with_capacity(8);
    for k in 0..8 {
        result.push(((x >> (k * 16)) & 65535) as u16);
    }
    result
}

/** Convert a length-8 vector of uint16's into a 128-bit integer
*/
pub fn bigbin_to_int(x: &Vec<u16>) -> u128 {
    x.iter()
        .enumerate()
        .fold(0, |acc, (i, &v)| acc | ((v as u128) << (i * 16)))
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
// original implementation
// pub fn uint16s_to_bits<T: ToU16>(data: &Vec<T>) -> Vec<u8> {
//     let mut result = Vec::with_capacity(data.len() * 16);

//     for value in data {
//         // Extract each bit from the 16-bit value
//         let value_u16 = value.to_u16();
//         for i in 0..16 {
//             result.push(((value_u16 >> i) & 1) as u8);
//         }
//     }
//     result
// }

// optimized implementation
pub fn uint16s_to_bits<T: ToU16>(data: &Vec<T>) -> Vec<u8> {
    let len = data.len() * 16;
    let mut result = Vec::with_capacity(len);

    // optimize trick: use unsafe code directly to avoid the overhead of bounds checking
    unsafe {
        result.set_len(len);
        let mut index = 0;
        for value in data {
            let value_u16 = value.to_u16();
            for i in 0..16 {
                *result.get_unchecked_mut(index) = ((value_u16 >> i) & 1) as u8;
                index += 1;
            }
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
// pub fn uint16_to_bit(value: &BinaryFieldElement16) -> Vec<u8> {
//     let mut result = Vec::with_capacity(2);
//     let byte1 = (value.value & 0xFF) as u8;
//     let byte2 = ((value.value >> 8) & 0xFF) as u8;

//     result.push(byte1);
//     result.push(byte2);
//     result
// }

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
    fn test_int_to_bigbin() {
        // int_to_bigbin(0x1234567890abcdef)= [52719 37035 22136  4660     0     0     0     0]
        let result = int_to_bigbin(0x1234567890abcdef);
        assert_eq!(result, vec![52719, 37035, 22136, 4660, 0, 0, 0, 0]);
    }

    #[test]
    fn test_bigbin_to_int() {
        // bigbin_to_int([52719 37035 22136  4660     0     0     0     0])= 0x1234567890abcdef
        let data = vec![52719, 37035, 22136, 4660, 0, 0, 0, 0];
        let result = bigbin_to_int(&data);
        assert_eq!(result, 0x1234567890abcdef);
    }

    // #[test]
    // fn test_big_mul() {
    //     // big_mul(int_to_bigbin(3**29), int_to_bigbin(5**29))= [46732 49627 26993 63626 14101 27237 21150     0]
    //     // let a = int_to_bigbin(3u128.pow(29));
    //     // let b = int_to_bigbin(5u128.pow(29));
    //     let a = 3u128.pow(29);
    //     let b = 5u128.pow(29);
    //     let result = big_mul(a, b);
    //     assert_eq!(
    //         result,
    //         vec![46732, 49627, 26993, 63626, 14101, 27237, 21150, 0]
    //     );
    // }

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
