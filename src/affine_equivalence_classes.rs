//! Operations on affine equivalence classes of boolean functions
//!
//! $g$ and $f$ $n$-variable Boolean functions are affine equivalent if $\forall x \in \mathbb{F}_2, g(x) = f(Dx + a) + bx + c$
//! for some D $\in \mathcal M_n(\mathbb{F}_2)$ invertible matrix, $a, b \in \mathbb{F}^n_2$ vectors and $c \in \mathbb{F}_2$
//!
//! All the equivalence classes representatives have been computed by Joanne Elizabeth Fuller, in [her thesis](https://eprints.qut.edu.au/15828/1/Joanne_Fuller_Thesis.pdf).

use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::Zero;
use crate::{utils, BooleanFunctionImpl, BooleanFunctionType, SmallBooleanFunction};

/// Representatives of all affine equivalence classes of boolean functions with 3 variables
pub const BOOLEAN_FUNCTIONS_3_VAR_AFFINE_EQ_CLASSES: [SmallBooleanFunction; 3] = [
    SmallBooleanFunction::from_truth_table_unchecked(0xaa, 3),
    SmallBooleanFunction::from_truth_table_unchecked(0xab, 3),
    SmallBooleanFunction::from_truth_table_unchecked(0xac, 3),
];

/// Representatives of all affine equivalence classes of boolean functions with 4 variables.
/// The last class is the class of bent functions.
pub const BOOLEAN_FUNCTIONS_4_VAR_AFFINE_EQ_CLASSES: [SmallBooleanFunction; 8] = [
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xab55, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xbb55, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xaba5, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xaaff, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xaba4, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xab12, 4),
    SmallBooleanFunction::from_truth_table_unchecked(0xac90, 4), // bent
];

/// Representatives of all affine equivalence classes of boolean functions with 5 variables
pub const BOOLEAN_FUNCTIONS_5_VAR_AFFINE_EQ_CLASSES: [SmallBooleanFunction; 48] = [
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55aa55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55ab55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa55bb55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa5dbb55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaaddbb55, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaa5dbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x2a5dbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaaddbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x2a5dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x6a5dbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x2addbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xa8ddbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaeddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x0a5dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8addda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xa8dd9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88ddbb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88ddbb11, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8c5dda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xa89d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8eddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xaefdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x025dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88ddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88dd9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xceddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x0eddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x425dbf51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8cddda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88dddb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x289d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x86fdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x88dddb71, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xcefdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x0efdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x288d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8cfdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8cdddb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x8ccdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x289d9b41, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x488ddb51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xccfdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x688d9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x288d9b41, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x288d1b41, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0xdcfdda51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x68ad9b51, 5),
    SmallBooleanFunction::from_truth_table_unchecked(0x688ddb51, 5),
];

/// Affine-equivalence factors between 2 $n$-variable Boolean functions.
///
/// Check the definition of affine-equivalence [here](crate::affine_equivalence_classes).
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AffineEquivalenceFactors {
    /// Invertible Boolean matrix $\in \mathcal M_n(\mathbb{F}_2)$, each `u32` element of the `Vec<u32>` data structure is a column of the matrix.
    pub d: Vec<u32>,
    /// Boolean vector $\in \mathbb{F}^n_2$
    pub a: u32,
    /// Boolean vector $\in \mathbb{F}^n_2$
    pub b: u32,
    /// Boolean
    pub c: bool,
}

/// Check if two Boolean functions are affine-equivalent, and returns equivalence factors if so.
///
/// Check the definition of affine-equivalence [here](crate::affine_equivalence_classes).
///
/// The computation could take time depending on the variable count and internal loops heuristic.
///
/// # Parameters
/// - `f`: A reference to the first BooleanFunction
/// - `g`: A reference to the second BooleanFunction
///
/// # Returns
/// The affine equivalence factors if `f` and `g` are affine-equivalent, `None` otherwise.
///
/// # Note
/// If `f` anf `g` don't have the same input variable count, the function will return `None`.
///
/// # Example
/// ```rust
/// use boolean_function::{BooleanFunction, SmallBooleanFunction};
/// use boolean_function::affine_equivalence_classes::compute_affine_equivalence;
/// // 4-variable bent functions equivalence class
/// let f = BooleanFunction::from_hex_string_truth_table("ac90").unwrap();
/// let g = SmallBooleanFunction::from_truth_table(0xdbd4, 4).unwrap();
///
/// // f and g are affine-equivalent
/// assert!(compute_affine_equivalence(&f, &g).is_some());
pub fn compute_affine_equivalence(f: &impl BooleanFunctionImpl, g: &impl BooleanFunctionImpl) -> Option<AffineEquivalenceFactors> {
    if f.variables_count() != g.variables_count() || f.variables_count() == 0 {
        return None;
    }
    let g_0_connecting = g.get_1_local_neighbor(0);

    let possible_a = (0u32..=f.get_max_input_value())
        .into_iter()
        .filter(|a| {
            let f_a_connecting = f.get_1_local_neighbor(*a);
            f_a_connecting.algebraic_degree() == g_0_connecting.algebraic_degree()
            && f_a_connecting.absolute_walsh_hadamard_spectrum() == g_0_connecting.absolute_walsh_hadamard_spectrum()
            && f_a_connecting.absolute_autocorrelation() == g_0_connecting.absolute_autocorrelation()
        });

    let possible_d_columns = (0..f.variables_count())
        .into_iter()
        .map(|column| {
            let g_j_connecting = g.get_1_local_neighbor(1 << column);
            (0u32..=f.get_max_input_value())
                .into_iter()
                .filter(|i| {
                    let f_i_connecting = f.get_1_local_neighbor(*i);
                    f_i_connecting.algebraic_degree() == g_j_connecting.algebraic_degree()
                    && f_i_connecting.absolute_walsh_hadamard_spectrum() == g_j_connecting.absolute_walsh_hadamard_spectrum()
                    && f_i_connecting.absolute_autocorrelation() == g_j_connecting.absolute_autocorrelation()
                })
                .collect::<Vec<_>>()
        });
    let all_matrix_candidates_iter = possible_d_columns.multi_cartesian_product();

    for a_candidate in possible_a {
        for d_a_matrix_candidate in all_matrix_candidates_iter.clone() {
            for b_candidate in 0..=f.get_max_input_value() {
                for c_candidate in [true, false] {
                    let mut test_tt_small = 0u64;
                    let mut test_tt_big = BigUint::zero();
                    for x in 0..=f.get_max_input_value() {
                        let mut new_x = 0u32;
                        for (index, d_a_col_candidate) in d_a_matrix_candidate.iter().enumerate() {
                            let d_col = d_a_col_candidate ^ a_candidate;
                            let x_bit = utils::fast_binary_dot_product(d_col, x) & 1;
                            new_x |= x_bit << index;
                        }
                        let new_x_a = new_x ^ a_candidate;
                        let f_eval_new_x_a = f.compute_cellular_automata_rule(new_x_a);
                        let computed_bit_value = f_eval_new_x_a ^ ((utils::fast_binary_dot_product(x, b_candidate) & 1) != 0) ^ c_candidate;
                        match f.get_boolean_function_type() {
                            BooleanFunctionType::Small => {
                                test_tt_small |= (computed_bit_value as u64) << x;
                            }
                            BooleanFunctionType::Big => {
                                test_tt_big.set_bit(x as u64, computed_bit_value);
                            }
                        }
                    }
                    let same_tt = match f.get_boolean_function_type() {
                        BooleanFunctionType::Small => test_tt_small == g.try_u64_truth_table().unwrap(),
                        BooleanFunctionType::Big => test_tt_big == g.biguint_truth_table()
                    };
                    if same_tt {
                        return Some(AffineEquivalenceFactors {
                            d: d_a_matrix_candidate.iter()
                                .map(|column| column ^ a_candidate)
                                .collect(),
                            a: a_candidate,
                            b: b_candidate,
                            c: c_candidate,
                        });
                    }
                }
            }
        }
    }
    None
}


#[cfg(test)]
mod tests {
    use crate::affine_equivalence_classes::compute_affine_equivalence;
    use crate::{BooleanFunction, SmallBooleanFunction};

    #[test]
    fn test_compute_affine_equivalence_not_same_var_count() {
        let f = BooleanFunction::from_hex_string_truth_table("aa").unwrap();
        let g = BooleanFunction::from_hex_string_truth_table("aaaa").unwrap();
        assert_eq!(compute_affine_equivalence(&f, &g), None);
    }

    #[test]
    fn test_compute_affine_equivalence() {
        let f = BooleanFunction::from_hex_string_truth_table("1234").unwrap();
        let g = BooleanFunction::from_hex_string_truth_table("1234").unwrap();
        let factors = compute_affine_equivalence(&f, &g);
        assert!(factors.is_some());
        let factors = factors.unwrap();
        assert_eq!(factors.d, [1, 2, 4, 8]); // Identity matrix
        assert_eq!(factors.a, 0);
        assert_eq!(factors.b, 0);
        assert_eq!(factors.c, false);

        let f = BooleanFunction::from_hex_string_truth_table("1234").unwrap();
        let g = BooleanFunction::from_hex_string_truth_table("edcb").unwrap(); // f = !g
        let factors = compute_affine_equivalence(&f, &g);
        assert!(factors.is_some());
        let factors = factors.unwrap();
        assert_eq!(factors.d, [1, 2, 4, 8]); // Identity matrix
        assert_eq!(factors.a, 0);
        assert_eq!(factors.b, 0);
        assert_eq!(factors.c, true);

        let f = BooleanFunction::from_hex_string_truth_table("aaaa").unwrap();
        let g = BooleanFunction::from_hex_string_truth_table("5555").unwrap();
        let factors = compute_affine_equivalence(&f, &g);
        assert!(factors.is_some());
        let factors = factors.unwrap();
        assert_eq!(factors.d, [0, 0, 0, 0]); // Identity matrix
        assert_eq!(factors.a, 0);
        assert_eq!(factors.b, 1);
        assert_eq!(factors.c, true);

        let f = BooleanFunction::from_hex_string_truth_table("aa55").unwrap();
        let g = BooleanFunction::from_hex_string_truth_table("ab55").unwrap();
        let factors = compute_affine_equivalence(&f, &g);
        assert!(factors.is_none());

        // Class of 4-variable bent functions
        let f = BooleanFunction::from_hex_string_truth_table("ac90").unwrap();
        let g = BooleanFunction::from_hex_string_truth_table("dbd4").unwrap();
        let factors = compute_affine_equivalence(&f, &g);
        assert!(factors.is_some());
        let factors = factors.unwrap();
        assert_eq!(factors.d, [1, 2, 6, 8]);
        assert_eq!(factors.a, 0);
        assert_eq!(factors.b, 10);
        assert_eq!(factors.c, false);
    }

    #[test]
    fn test_compute_affine_equivalence_small() {
        let f = SmallBooleanFunction::from_truth_table(0x1234, 4).unwrap();
        let g = BooleanFunction::from_hex_string_truth_table("1234").unwrap();
        let factors = compute_affine_equivalence(&f, &g);
        assert!(factors.is_some());
        let factors = factors.unwrap();
        assert_eq!(factors.d, [1, 2, 4, 8]); // Identity matrix
        assert_eq!(factors.a, 0);
        assert_eq!(factors.b, 0);
        assert_eq!(factors.c, false);
    }

    /*#[test]
    fn test_compute_affine_equivalence_big() {
        let f = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("ffffffffffffffffffffffffffffffff", 16).unwrap(), 7);
        let g = BooleanFunction::from_hex_string_truth_table("00000000000000000000000000000000").unwrap();
        let factors = compute_affine_equivalence(&f, &g);
        assert!(factors.is_some());
        let factors = factors.unwrap();
        assert_eq!(factors.d, [0, 0, 0, 0, 0, 0, 0]); // Identity matrix
        assert_eq!(factors.a, 0);
        assert_eq!(factors.b, 0);
        assert_eq!(factors.c, true);
    }*/
}