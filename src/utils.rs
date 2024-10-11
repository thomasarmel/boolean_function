use num_bigint::BigUint;
use num_traits::{One, Zero};

#[inline]
pub(crate) fn fast_binary_dot_product(a: u32, b: u32) -> u32 {
    (a & b).count_ones()
}

fn boolean_gaussian_elimination(matrix: &mut Vec<Vec<bool>>) -> Vec<Vec<bool>> {
    let m = matrix.len();
    let n = matrix[0].len();

    let mut row = 0;

    // Gaussian elimination over Boolean algebra
    for col in 0..n {
        if row >= m {
            break;
        }

        // Find the pivot
        if !matrix[row][col] {
            let mut found = false;
            for i in row + 1..m {
                if matrix[i][col] {
                    matrix.swap(row, i);
                    found = true;
                    break;
                }
            }
            if !found {
                continue;
            }
        }

        // Eliminate the rest of the column
        for i in 0..m {
            if i != row && matrix[i][col] {
                for j in col..n {
                    matrix[i][j] ^= matrix[row][j]; // XOR for Boolean addition (mod 2)
                }
            }
        }
        row += 1;
    }

    matrix.clone()
}

pub(crate) fn left_kernel_boolean(matrix: &Vec<Vec<bool>>) -> Vec<Vec<bool>> {
    if matrix.is_empty() || matrix[0].is_empty() {
        return Vec::new();
    }
    let mut matrix_transpose = matrix.iter().cloned().zip(0..).fold(
        vec![vec![false; matrix.len()]; matrix[0].len()],
        |mut acc, (row, i)| {
            for (j, &val) in row.iter().enumerate() {
                acc[j][i] = val;
            }
            acc
        },
    );

    let m = matrix_transpose.len();
    let n = matrix_transpose[0].len();

    // Apply Gaussian elimination
    let rref = boolean_gaussian_elimination(&mut matrix_transpose);

    // Identify free variables (columns without a leading 1)
    let mut free_variables = Vec::new();
    let mut pivot_columns = vec![false; n];

    let mut row = 0;
    for col in 0..n {
        if row < m && rref[row][col] {
            pivot_columns[col] = true;
            row += 1;
        } else {
            free_variables.push(col);
        }
    }

    // Construct the left kernel from free variables
    let mut left_kernel = Vec::new();

    for &free_col in &free_variables {
        let mut kernel_vector = vec![false; n];
        kernel_vector[free_col] = true;

        // Back-substitute to fill in the rest of the kernel vector
        for row in (0..m).rev() {
            let pivot_col = rref[row].iter().position(|&x| x).unwrap_or(n);
            if pivot_col < n {
                kernel_vector[pivot_col] = rref[row]
                    .iter()
                    .enumerate()
                    .skip(pivot_col + 1)
                    .any(|(j, &val)| val && kernel_vector[j]);
            }
        }

        left_kernel.push(kernel_vector);
    }

    left_kernel
}

pub(crate) fn fast_anf_transform_biguint(truth_table: &BigUint, variables_count: usize) -> BigUint {
    let u0 = BigUint::zero();
    let u1 = BigUint::one();

    let mut blocksize = 1usize;
    let mut anf_form = truth_table.clone();
    for _ in 0..variables_count {
        let mut source = 0usize;
        while source < (1 << variables_count) {
            let target = source + blocksize.clone();
            for i in 0..blocksize {
                let f_target_i: bool = ((anf_form.clone() >> (target + i)) & u1.clone()) != u0;
                let f_source_i: bool = ((anf_form.clone() >> (source + i)) & u1.clone()) != u0;
                let f_target_i_xor_f_source_i = f_target_i ^ f_source_i;
                if f_target_i_xor_f_source_i {
                    anf_form = anf_form | (u1.clone() << (target + i));
                } else {
                    // set (target + i) bit of final_f to 0 with not operation
                    anf_form.set_bit((target + i) as u64, false);
                }
            }
            source = source + (blocksize << 1);
        }
        blocksize = blocksize << 1;
    }
    anf_form
}

#[allow(dead_code)] // maybe useless, but I keep it for the beauty of the code
pub(crate) fn walsh_matrix(dim: usize) -> Vec<Vec<i8>> {
    (0usize..(1 << dim))
        .map(|y| {
            (0..(1 << dim))
                .map(|x| if (x & y).count_ones() & 1 == 0 { 1 } else { -1 })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::utils::left_kernel_boolean;

    #[test]
    fn test_left_kernel_boolean() {
        let matrix = vec![
            vec![true, false, false, false],
            vec![true, false, false, false],
            vec![false, false, false, false],
            vec![false, false, false, false],
            vec![false, false, false, false],
            vec![false, false, false, false],
            vec![false, false, false, false],
        ];

        let left_kernel = left_kernel_boolean(&matrix);

        assert_eq!(
            left_kernel,
            [
                [true, true, false, false, false, false, false],
                [false, false, true, false, false, false, false],
                [false, false, false, true, false, false, false],
                [false, false, false, false, true, false, false],
                [false, false, false, false, false, true, false],
                [false, false, false, false, false, false, true]
            ]
        );
    }

    #[test]
    fn test_walsh_matrix() {
        assert_eq!(super::walsh_matrix(0), vec![vec![1]]);
        assert_eq!(super::walsh_matrix(1), vec![vec![1, 1], vec![1, -1]]);
        assert_eq!(
            super::walsh_matrix(2),
            vec![
                vec![1, 1, 1, 1],
                vec![1, -1, 1, -1],
                vec![1, 1, -1, -1],
                vec![1, -1, -1, 1]
            ]
        );
        assert_eq!(
            super::walsh_matrix(3),
            vec![
                vec![1, 1, 1, 1, 1, 1, 1, 1],
                vec![1, -1, 1, -1, 1, -1, 1, -1],
                vec![1, 1, -1, -1, 1, 1, -1, -1],
                vec![1, -1, -1, 1, 1, -1, -1, 1],
                vec![1, 1, 1, 1, -1, -1, -1, -1],
                vec![1, -1, 1, -1, -1, 1, -1, 1],
                vec![1, 1, -1, -1, -1, -1, 1, 1],
                vec![1, -1, -1, 1, -1, 1, 1, -1]
            ]
        );
    }
}
