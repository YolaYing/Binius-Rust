//! This module provides functionality for Merkle trees.
//!
//! The module provide the following functions:
//! 1. hash: hash a byte array using SHA256
//! 2. merkelize: build a Merkle tree from the inputs
//! 3. get_root: return the root of the Merkle tree
//! 4. get_branch: get the branch of the Merkle tree
//! 5. verify_branch: verify the Merkle branch

use sha2::{Digest, Sha256};

pub fn hash(x: &[u8]) -> Vec<u8> {
    let mut hasher = Sha256::new();
    hasher.update(x);
    hasher.finalize().to_vec()
}

/** Build a Merkle tree from the inputs

where o[i] is the parent node of o[2i] and o[2i+1], the second half of o is the original data, and o[1] is the root

Args:
    vals: the original data, should be packed_column

Returns:
    the Merkle tree
*/
pub fn merkelize(vals: &Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    assert_eq!(vals.len() & (vals.len() - 1), 0);
    let mut o = vec![vec![]; vals.len() * 2];
    for (i, x) in vals.iter().enumerate() {
        o[vals.len() + i] = hash(x);
    }
    for i in (1..vals.len()).rev() {
        let mut combined = o[i * 2].clone();
        combined.extend(o[i * 2 + 1].clone());
        o[i] = hash(&combined);
    }
    o
}

/** return the root of the Merkle tree

Args:
    tree: the Merkle tree

Returns:
    the root of the Merkle tree, o[1](the first element of the tree is None)
*/
pub fn get_root(tree: &Vec<Vec<u8>>) -> Vec<u8> {
    tree[1].clone()
}

/** Get the branch of the Merkle tree

the Merkle tree hash path from the leaf to the root, the branch is the sibling of the path

Args:
    tree: the Merkle tree
    pos: the position of the leaf

Returns:
    the hash path of the Merkle tree
 */
pub fn get_branch(tree: &Vec<Vec<u8>>, pos: usize) -> Vec<Vec<u8>> {
    let offset_pos = pos + tree.len() / 2;
    let branch_length = (tree.len() as f64).log2() as usize - 1;
    let mut branch = vec![];
    for i in 0..branch_length {
        branch.push(tree[(offset_pos >> i) ^ 1].clone());
    }
    branch
}

// # Verify that Merkle branch (requires only the root, not the tree)
// def verify_branch(root, pos, val, branch):
//     x = hash(val)
//     for b in branch:
//         if pos & 1:
//             x = hash(b + x)
//         else:
//             x = hash(x + b)
//         pos //= 2
//     return x == root
pub fn verify_branch(root: &[u8], pos: usize, val: &[u8], branch: &Vec<Vec<u8>>) -> bool {
    let mut x = hash(val);
    let mut pos = pos;
    for b in branch {
        if pos & 1 == 1 {
            x = hash(&[b.as_slice(), x.as_slice()].concat());
        } else {
            x = hash(&[x.as_slice(), b.as_slice()].concat());
        }
        pos /= 2;
    }
    x == root
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash() {
        let x = vec![1, 2, 3];
        let result = hash(&x);
        assert_eq!(
            result,
            vec![
                0x03, 0x90, 0x58, 0xc6, 0xf2, 0xc0, 0xcb, 0x49, 0x2c, 0x53, 0x3b, 0x0a, 0x4d, 0x14,
                0xef, 0x77, 0xcc, 0x0f, 0x78, 0xab, 0xcc, 0xce, 0xd5, 0x28, 0x7d, 0x84, 0xa1, 0xa2,
                0x01, 0x1c, 0xfb, 0x81
            ]
        );
    }

    #[test]
    fn test_merkelize() {
        let vals = vec![vec![1, 2], vec![3, 4]];
        let result = merkelize(&vals);
        assert_eq!(result[0], Vec::<u8>::new());
        assert_eq!(
            result[1],
            vec![
                0xbe, 0xd3, 0xd3, 0x3a, 0x81, 0x02, 0x6f, 0x7b, 0xe9, 0x3a, 0xef, 0xad, 0x44, 0xc5,
                0x89, 0x1c, 0x27, 0xfc, 0x82, 0x65, 0xaa, 0x15, 0x27, 0x9a, 0x58, 0xe2, 0x87, 0x74,
                0x4b, 0x7c, 0x77, 0x53
            ]
        );
    }

    #[test]
    fn test_get_root() {
        let tree = vec![
            vec![],
            vec![
                0xbe, 0xd3, 0xd3, 0x3a, 0x81, 0x02, 0x6f, 0x7b, 0xe9, 0x3a, 0xef, 0xad, 0x44, 0xc5,
                0x89, 0x1c, 0x27, 0xfc, 0x82, 0x65, 0xaa, 0x15, 0x27, 0x9a, 0x58, 0xe2, 0x87, 0x74,
                0x4b, 0x7c, 0x77, 0x53,
            ],
        ];
        let result = get_root(&tree);
        assert_eq!(
            result,
            vec![
                0xbe, 0xd3, 0xd3, 0x3a, 0x81, 0x02, 0x6f, 0x7b, 0xe9, 0x3a, 0xef, 0xad, 0x44, 0xc5,
                0x89, 0x1c, 0x27, 0xfc, 0x82, 0x65, 0xaa, 0x15, 0x27, 0x9a, 0x58, 0xe2, 0x87, 0x74,
                0x4b, 0x7c, 0x77, 0x53
            ]
        );
    }

    #[test]
    fn test_verify_branch() {
        let vals = vec![vec![1, 2], vec![3, 4]];
        let tree = merkelize(&vals);
        let pos = 1;
        let branch = get_branch(&tree, pos);
        let result = verify_branch(&tree[1], pos, &vals[1], &branch);
        assert_eq!(result, true);
    }
}
