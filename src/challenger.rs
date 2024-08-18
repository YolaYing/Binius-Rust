use super::merkle_tree::hash;

/** Get challenges from the root of the Merkle tree

Args:
    root: the root of the Merkle tree
    extended_row_length: the length of the extended row
    num_challenges: the number of challenges

Returns:
    Vec<u16>: the challenges, indexes of the columns
*/
pub fn get_challenges(root: &[u8], extended_row_length: usize, num_challenges: usize) -> Vec<u16> {
    let mut o = vec![];
    for i in 0..num_challenges {
        let mut bytes = root.to_vec();
        bytes.push(i as u8);
        let hash = hash(&bytes);
        let challenge =
            u16::from_le_bytes(hash[0..2].try_into().unwrap()) % extended_row_length as u16;
        o.push(challenge);
    }
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_challenges() {
        let root = vec![1, 2, 3, 4];
        let extended_row_length = 8;
        let num_challenges = 2;
        let result = get_challenges(&root, extended_row_length, num_challenges);
        assert_eq!(result, vec![6, 0]);
    }
}
