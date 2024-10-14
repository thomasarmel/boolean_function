use crate::anf_polynom::AnfPolynomial;
#[cfg(not(feature = "unsafe_disable_safety_checks"))]
use crate::boolean_function_error::XOR_DIFFERENT_VAR_COUNT_PANIC_MSG;
use crate::iterator::BooleanFunctionIterator;
use crate::utils::left_kernel_boolean;
#[cfg(not(feature = "unsafe_disable_safety_checks"))]
use crate::BooleanFunctionError::TooBigTruthTableForVarCount;
use crate::BooleanFunctionError::TooBigVariableCount;
use crate::{BooleanFunction, BooleanFunctionError, BooleanFunctionImpl, BooleanFunctionType};
use fast_boolean_anf_transform::fast_bool_anf_transform_unsigned;
use itertools::{enumerate, Itertools};
use num_bigint::BigUint;
use num_integer::binomial;
use std::ops::{BitXor, BitXorAssign, Not};

/// Struct representing a boolean function with a big truth table.
///
/// The struct internally stores the truth table as an u64, meaning that the maximum number of variables is 6.
///
/// For more than 6 variables, use the [crate::BigBooleanFunction] struct, or the [crate::BooleanFunction] type to store both small and big boolean functions.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct SmallBooleanFunction {
    variables_count: usize,
    truth_table: u64,
}

impl SmallBooleanFunction {
    /// Creates a new [SmallBooleanFunction] from a truth table and the number of variables.
    ///
    /// # Parameters
    /// - `truth_table` - The truth table of the Boolean function, where the lower bit represents the output of the Boolean function for the input 0.
    /// - `variables_count` - The number of variables of the Boolean function.
    ///
    /// # Returns
    /// A [SmallBooleanFunction] instance from the truth table and the number of variables or an error if:
    /// - The number of variables is greater than 6.
    /// - The truth table is too big for the number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
    pub fn from_truth_table(
        truth_table: u64,
        variables_count: usize,
    ) -> Result<Self, BooleanFunctionError> {
        if variables_count > 6 {
            return Err(TooBigVariableCount(6));
        }
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if truth_table > (u64::MAX >> (64 - (1 << variables_count))) {
            return Err(TooBigTruthTableForVarCount);
        }
        Ok(SmallBooleanFunction {
            variables_count,
            truth_table,
        })
    }

    pub(crate) const fn from_truth_table_unchecked(
        truth_table: u64,
        variables_count: usize,
    ) -> Self {
        SmallBooleanFunction {
            variables_count,
            truth_table,
        }
    }

    /// Returns the truth table of the Boolean function, as an [u64].
    pub fn get_truth_table_u64(&self) -> u64 {
        self.truth_table
    }

    /// Computes the [derivative](crate::BooleanFunctionImpl::derivative) of the Boolean function for a given direction.
    ///
    /// # Parameters
    /// * `direction` - The direction of the derivative.
    ///
    /// # Returns
    /// The derivative of the Boolean function for the given direction, or an error if the direction is greater than the maximum input value and the `unsafe_disable_safety_checks` feature is not enabled.
    pub fn derivative_inner(&self, direction: u32) -> Result<Self, BooleanFunctionError> {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        {
            let max_input_value = self.get_max_input_value();
            if direction > max_input_value {
                return Err(BooleanFunctionError::TooBigDerivativeDirection(
                    max_input_value,
                ));
            }
        }
        let mut derivative_truth_table: u64 = 0;
        for x in 0..=self.get_max_input_value() {
            if self.compute_cellular_automata_rule(x)
                ^ self.compute_cellular_automata_rule(x ^ direction)
            {
                derivative_truth_table |= 1 << x;
            }
        }
        Ok(SmallBooleanFunction {
            variables_count: self.variables_count,
            truth_table: derivative_truth_table,
        })
    }

    /// Computes the [reverse](crate::BooleanFunctionImpl::reverse) of the Boolean function.
    ///
    /// # Returns
    /// The reverse of the Boolean function.
    pub fn reverse_inner(&self) -> Self {
        Self {
            variables_count: self.variables_count,
            truth_table: !self.truth_table & (u64::MAX >> (64 - (1 << self.variables_count))),
        }
    }

    /// Computes the [annihilator](crate::BooleanFunctionImpl::annihilator) of the Boolean function for a given maximum degree.
    ///
    /// # Parameters
    /// * `max_degree` - The maximum degree of the wished annihilator.
    ///
    /// # Returns
    /// A tuple containing the annihilator function, its degree and the dimension of the annihilator vector space, or `None` no annihilator was found.
    pub fn annihilator_inner(
        &self,
        max_degree: usize,
    ) -> Option<(SmallBooleanFunction, usize, usize)> {
        if self.truth_table == 0 {
            let max_possible_function_tt = u64::MAX >> (64 - (1 << self.variables_count));
            let dim_annihilator_vec_space = (0..=max_degree)
                .map(|i| binomial(self.variables_count as u64, i as u64))
                .sum::<u64>() as usize;
            return Some((
                Self::from_truth_table(max_possible_function_tt, self.variables_count).unwrap(),
                0,
                dim_annihilator_vec_space,
            ));
        }

        let truth_table_non_zero_positions = (0u32..(1 << self.variables_count))
            .filter(|bit_pos| self.truth_table & (1u64 << bit_pos) != 0)
            .collect::<Vec<u32>>();

        let matrix_out_len = (0..=max_degree)
            .map(|i| binomial(self.variables_count as u64, i as u64))
            .sum::<u64>() as usize;
        let matrix_in_len = truth_table_non_zero_positions.len();
        let mut matrix: Vec<Vec<bool>> = vec![vec![false; matrix_in_len]; matrix_out_len];

        let mut r = [AnfPolynomial::from_anf_small(0x1, self.variables_count)].to_vec();

        for i in 1..=max_degree {
            for comb in (0..self.variables_count).combinations(i) {
                let mut bit_index = 0;
                for monomial in comb {
                    bit_index |= 1 << monomial;
                }
                let anf = 1u64 << bit_index;
                r.push(AnfPolynomial::from_anf_small(anf, self.variables_count));
            }
        }

        for (i, m) in enumerate(r.iter()) {
            let truth_table = fast_bool_anf_transform_unsigned(
                m.get_polynomial_small().unwrap(),
                self.variables_count,
            );
            for (j, v) in enumerate(truth_table_non_zero_positions.iter()) {
                matrix[i][j] = truth_table & (1u64 << v) != 0;
            }
        }

        let left_kernel = left_kernel_boolean(&matrix);

        if left_kernel.is_empty() {
            return None;
        }

        let annihilator_anf = enumerate(r.iter())
            .filter(|(i, _)| left_kernel[0][*i])
            .map(|(_, v)| v.get_polynomial_small().unwrap())
            .sum();

        let annihilator_function = Self::from_truth_table(
            fast_bool_anf_transform_unsigned(annihilator_anf, self.variables_count),
            self.variables_count,
        )
        .unwrap();

        let annihilator_degree =
            AnfPolynomial::from_anf_small(annihilator_anf, self.variables_count).get_degree();

        Some((annihilator_function, annihilator_degree, left_kernel.len()))
    }

    /// Computes a [SmallBooleanFunction] from [Walsh-Hadamard values](crate::BooleanFunctionImpl::walsh_hadamard_values), by applying the inverse Walsh-Hadamard transform.
    ///
    /// # Parameters
    /// * `walsh_values` - The Walsh-Hadamard values of the Boolean function.
    ///
    /// # Returns
    /// The Boolean function created from the Walsh-Hadamard values list, or an error if:
    /// - the list length is less than 4
    /// - the list length is not a power of 2.
    /// - the number of variables is greater than 6 (meaning that the list length is greater than 64).
    pub fn from_walsh_hadamard_values(walsh_values: &[i32]) -> Result<Self, BooleanFunctionError> {
        let walsh_values_count = walsh_values.len();
        if walsh_values_count < 4 || walsh_values_count.count_ones() != 1 {
            return Err(BooleanFunctionError::InvalidWalshValuesCount(
                walsh_values_count,
            ));
        }
        let num_variables = walsh_values_count.trailing_zeros() as usize;
        if num_variables > 6 {
            return Err(TooBigVariableCount(6));
        }
        let mut truth_table = 0u64;
        for i in 0..(1 << num_variables) {
            let value = walsh_values
                .iter()
                .enumerate()
                .map(|(w, walsh_value)| {
                    walsh_value * (if (w & i).count_ones() & 1 == 0 { 1 } else { -1 })
                })
                .sum::<i32>()
                < 0;
            if value {
                truth_table |= 1 << i;
            }
        }
        Ok(Self {
            variables_count: num_variables,
            truth_table,
        })
    }

    /// Computes a [SmallBooleanFunction] from [Walsh-Fourier values](crate::BooleanFunctionImpl::walsh_fourier_values), by applying the inverse Walsh-Fourier transform.
    ///
    /// # Parameters
    /// * `walsh_values` - The Walsh-Hadamard values of the Boolean function.
    ///
    /// # Returns
    /// The Boolean function created from the Walsh-Hadamard values list, or an error if:
    /// - the list length is less than 4
    /// - the list length is not a power of 2.
    /// - the number of variables is greater than 6 (meaning that the list length is greater than 64).
    pub fn from_walsh_fourier_values(walsh_values: &[i32]) -> Result<Self, BooleanFunctionError> {
        let walsh_values_count = walsh_values.len();
        if walsh_values_count < 4 || walsh_values_count.count_ones() != 1 {
            return Err(BooleanFunctionError::InvalidWalshValuesCount(
                walsh_values_count,
            ));
        }
        let num_variables = walsh_values_count.trailing_zeros() as usize;
        if num_variables > 6 {
            return Err(TooBigVariableCount(6));
        }
        let mut truth_table = 0u64;
        for i in 0..(1 << num_variables) {
            let value = walsh_values
                .iter()
                .enumerate()
                .map(|(w, walsh_value)| {
                    walsh_value * (if (w & i).count_ones() & 1 == 0 { 1 } else { -1 })
                })
                .sum::<i32>()
                != 0;
            if value {
                truth_table |= 1 << i;
            }
        }
        Ok(Self {
            variables_count: num_variables,
            truth_table,
        })
    }
}

impl BooleanFunctionImpl for SmallBooleanFunction {
    #[inline]
    fn variables_count(&self) -> usize {
        self.variables_count
    }

    fn get_boolean_function_type(&self) -> BooleanFunctionType {
        BooleanFunctionType::Small
    }

    fn is_balanced(&self) -> bool {
        let expected_set_number: u32 = 1 << (self.variables_count - 1);
        self.truth_table.count_ones() == expected_set_number
    }

    #[inline]
    fn compute_cellular_automata_rule(&self, input_bits: u32) -> bool {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        {
            let max_input_value = self.get_max_input_value();
            if input_bits > max_input_value {
                panic!("Input bits must be less or equal than {}", max_input_value);
            }
        }
        (self.truth_table & (1 << input_bits)) != 0
    }

    fn derivative(&self, direction: u32) -> Result<BooleanFunction, BooleanFunctionError> {
        Ok((self.derivative_inner(direction)?).into())
    }

    fn is_linear(&self) -> bool {
        let max_input_value = self.get_max_input_value();
        [self.truth_table, self.reverse_inner().truth_table]
            .iter()
            .any(|rule| {
                let mut equivalent_xor_function: u64 = 0;
                for i in 0..=max_input_value {
                    let mut equivalent_xor_function_eval_i = false;
                    for j in 0..self.variables_count {
                        if *rule & (1 << (1 << j)) != 0 {
                            equivalent_xor_function_eval_i ^= (i & (1 << j)) == 0;
                        }
                    }
                    equivalent_xor_function |= (equivalent_xor_function_eval_i as u64) << i;
                }
                *rule == equivalent_xor_function
                    || *rule
                        == (Self {
                            variables_count: self.variables_count,
                            truth_table: equivalent_xor_function,
                        })
                        .reverse_inner()
                        .truth_table
            })
    }

    fn reverse(&self) -> BooleanFunction {
        (self.reverse_inner()).into()
    }

    fn algebraic_normal_form(&self) -> AnfPolynomial {
        let anf_form = fast_bool_anf_transform_unsigned(self.truth_table, self.variables_count);
        AnfPolynomial::from_anf_small(anf_form, self.variables_count)
    }

    fn annihilator(&self, max_degree: usize) -> Option<(BooleanFunction, usize, usize)> {
        let annihilator = self.annihilator_inner(max_degree)?;
        Some(((annihilator.0).into(), annihilator.1, annihilator.2))
    }

    fn iter(&self) -> BooleanFunctionIterator {
        BooleanFunctionIterator::new((self.clone()).into())
    }

    fn printable_hex_truth_table(&self) -> String {
        format!("{:01$x}", self.truth_table, 1 << (self.variables_count - 2))
    }

    fn biguint_truth_table(&self) -> BigUint {
        BigUint::from(self.truth_table)
    }
}

/// In-place XOR operator for Boolean functions truth tables.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitXorAssign for SmallBooleanFunction {
    fn bitxor_assign(&mut self, rhs: Self) {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if self.variables_count != rhs.variables_count {
            panic!("{}", XOR_DIFFERENT_VAR_COUNT_PANIC_MSG);
        }
        self.truth_table ^= rhs.truth_table;
    }
}

/// XOR operator for Boolean functions truth tables.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitXor for SmallBooleanFunction {
    type Output = Self;

    fn bitxor(mut self, rhs: Self) -> Self::Output {
        self ^= rhs;
        self
    }
}

/// NOT operator for Boolean functions.
///
/// This is equivalent to the [crate::BooleanFunctionImpl::reverse] operation: it reverses each output of the Boolean function.
impl Not for SmallBooleanFunction {
    type Output = Self;

    fn not(self) -> Self::Output {
        self.reverse_inner()
    }
}

#[cfg(test)]
mod tests {
    use crate::{BooleanFunctionImpl, SmallBooleanFunction};

    #[test]
    fn test_from_truth_table() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xaa55aa55);

        let boolean_function = SmallBooleanFunction::from_truth_table(0xffffffffffffffff, 7);
        assert!(boolean_function.is_err());
    }

    #[test]
    fn test_variables_count() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert_eq!(boolean_function.variables_count(), 5);
    }

    #[test]
    fn test_is_balanced() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert!(boolean_function.is_balanced());

        let boolean_function = SmallBooleanFunction::from_truth_table(0x17ffe80, 5).unwrap();
        assert!(boolean_function.is_balanced());

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaaaaaaaa, 5).unwrap();
        assert!(boolean_function.is_balanced());

        let boolean_function = SmallBooleanFunction::from_truth_table(0xabce1234, 5).unwrap();
        assert!(!boolean_function.is_balanced());
    }

    #[test]
    fn test_compute_cellular_automata_rule() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xabce1234, 5).unwrap();
        assert_eq!(boolean_function.compute_cellular_automata_rule(0), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(1), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(4), true);
        assert_eq!(boolean_function.compute_cellular_automata_rule(8), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(23), true);
    }

    #[test]
    fn test_derivative() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        let derivative = boolean_function.derivative_inner(1).unwrap();
        assert_eq!(derivative.variables_count(), 5);
        assert_eq!(derivative.get_truth_table_u64(), 0xffffffff);

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaaddbb55, 5).unwrap();
        let derivative = boolean_function.derivative_inner(1).unwrap();
        assert_eq!(derivative.get_truth_table_u64(), 0xff33ccff);

        let derivative = boolean_function.derivative_inner(2).unwrap();
        assert_eq!(derivative.get_truth_table_u64(), 0x00aa5500);

        let derivative = boolean_function.derivative_inner(3).unwrap();
        assert_eq!(derivative.get_truth_table_u64(), 0xff6666ff);

        let derivative = boolean_function.derivative_inner(32);
        assert!(derivative.is_err());
    }

    #[test]
    fn test_printable_hex_truth_table() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "aa55aa55");

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55, 5).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "0000aa55");

        let boolean_function = SmallBooleanFunction::from_truth_table(0, 5).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00000000");
    }

    #[test]
    fn test_algebraic_normal_form() {
        let boolean_function = SmallBooleanFunction::from_truth_table(30, 3).unwrap();
        assert_eq!(
            boolean_function
                .algebraic_normal_form()
                .get_polynomial_small()
                .unwrap(),
            30
        );

        let boolean_function = SmallBooleanFunction::from_truth_table(0, 3).unwrap();
        assert_eq!(
            boolean_function
                .algebraic_normal_form()
                .get_polynomial_small()
                .unwrap(),
            0
        );

        let boolean_function = SmallBooleanFunction::from_truth_table(110, 3).unwrap();
        assert_eq!(
            boolean_function
                .algebraic_normal_form()
                .get_polynomial_small()
                .unwrap(),
            142
        );

        let boolean_function = SmallBooleanFunction::from_truth_table(30, 3).unwrap();
        assert_eq!(
            boolean_function.algebraic_normal_form().to_string(),
            "x0*x1 + x0 + x1 + x2"
        );

        let boolean_function = SmallBooleanFunction::from_truth_table(0, 3).unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "0");

        let boolean_function = SmallBooleanFunction::from_truth_table(110, 3).unwrap();
        assert_eq!(
            boolean_function.algebraic_normal_form().to_string(),
            "x0*x1*x2 + x0*x1 + x0 + x1"
        );

        let boolean_function = SmallBooleanFunction::from_truth_table(123, 3).unwrap();
        assert_eq!(
            boolean_function.algebraic_normal_form().to_string(),
            "x0*x1 + x1*x2 + x1 + 1"
        );
    }

    #[test]
    fn test_reverse_inner() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        let reversed_boolean_function = boolean_function.reverse_inner();
        assert_eq!(reversed_boolean_function.get_truth_table_u64(), 0x55aa55aa);
        assert_eq!(reversed_boolean_function.variables_count(), 5);

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 6).unwrap();
        let reversed_boolean_function = boolean_function.reverse_inner();
        assert_eq!(
            reversed_boolean_function.get_truth_table_u64(),
            0xffffffff55aa55aa
        );
        assert_eq!(reversed_boolean_function.variables_count(), 6);

        let boolean_function = SmallBooleanFunction::from_truth_table(0, 6).unwrap();
        let reversed_boolean_function = boolean_function.reverse_inner();
        assert_eq!(
            reversed_boolean_function.get_truth_table_u64(),
            0xffffffffffffffff
        );

        let boolean_function =
            SmallBooleanFunction::from_truth_table(0xffffffffffffffff, 6).unwrap();
        let reversed_boolean_function = boolean_function.reverse_inner();
        assert_eq!(reversed_boolean_function.get_truth_table_u64(), 0);

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        let reversed_boolean_function = !boolean_function;
        assert_eq!(reversed_boolean_function.get_truth_table_u64(), 0x55aa55aa);
        assert_eq!(reversed_boolean_function.variables_count(), 5);
    }

    #[test]
    fn test_is_linear() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function = SmallBooleanFunction::from_truth_table(0x17ffe80, 5).unwrap();
        assert!(!boolean_function.is_linear());

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaaaaaaaa, 5).unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function = SmallBooleanFunction::from_truth_table(0xabce1234, 5).unwrap();
        assert!(!boolean_function.is_linear());

        let boolean_function = SmallBooleanFunction::from_truth_table(0x00000000, 5).unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function = SmallBooleanFunction::from_truth_table(0xffffffff, 5).unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function =
            SmallBooleanFunction::from_truth_table(0x0000000000000000, 6).unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function =
            SmallBooleanFunction::from_truth_table(0xffffffffffffffff, 6).unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function =
            SmallBooleanFunction::from_truth_table(0x00000000ffffffff, 6).unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function =
            SmallBooleanFunction::from_truth_table(0xabcdef0123456789, 6).unwrap();
        assert!(!boolean_function.is_linear());
    }

    #[test]
    fn test_annihilator_inner() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        let annihilator = boolean_function.annihilator_inner(1).unwrap();
        assert_eq!(annihilator.0.get_truth_table_u64(), 0x55aa55aa);
        assert_eq!(annihilator.1, 1);
        assert_eq!(annihilator.2, 1);

        let annihilator = boolean_function.annihilator_inner(6).unwrap();
        assert_eq!(annihilator.0.get_truth_table_u64(), 0x55aa55aa);
        assert_eq!(annihilator.1, 1);
        assert_eq!(annihilator.2, 16);

        let boolean_function = SmallBooleanFunction::from_truth_table(0x1e, 3).unwrap();
        let annihilator = boolean_function.annihilator_inner(1);
        assert!(annihilator.is_none());

        let annihilator = boolean_function.annihilator_inner(2).unwrap();
        assert_eq!(annihilator.0.get_truth_table_u64(), 0xe1);
        assert_eq!(annihilator.1, 2);
        assert_eq!(annihilator.2, 3);

        let boolean_function = SmallBooleanFunction::from_truth_table(0, 4).unwrap();
        let annihilator = boolean_function.annihilator_inner(0).unwrap();
        assert_eq!(annihilator.0.get_truth_table_u64(), 0xffff);
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 1);

        let annihilator = boolean_function.annihilator_inner(3).unwrap();
        assert_eq!(annihilator.0.get_truth_table_u64(), 0xffff);
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 15);

        let boolean_function =
            SmallBooleanFunction::from_truth_table(0xffffffffffffffff, 6).unwrap();
        let annihilator = boolean_function.annihilator_inner(6);
        assert!(annihilator.is_none());

        let boolean_function =
            SmallBooleanFunction::from_truth_table(0xabcdef0123456789, 6).unwrap();
        let annihilator = boolean_function.annihilator_inner(2);
        assert!(annihilator.is_none());

        let annihilator = boolean_function.annihilator_inner(3);
        assert_eq!(
            annihilator.unwrap().0.get_truth_table_u64(),
            0x1010101010101010
        );
        assert_eq!(annihilator.unwrap().1, 3);
        assert_eq!(annihilator.unwrap().2, 10);

        let annihilator = boolean_function.annihilator_inner(4);
        assert_eq!(
            annihilator.unwrap().0.get_truth_table_u64(),
            0x1010101010101010
        );
        assert_eq!(annihilator.unwrap().1, 3);
        assert_eq!(annihilator.unwrap().2, 25);

        let annihilator = boolean_function.annihilator_inner(5);
        assert_eq!(
            annihilator.unwrap().0.get_truth_table_u64(),
            0x1010101010101010
        );
        assert_eq!(annihilator.unwrap().1, 3);
        assert_eq!(annihilator.unwrap().2, 31);
    }

    #[test]
    fn test_from_walsh_hadamard_values() {
        let boolean_function =
            SmallBooleanFunction::from_walsh_hadamard_values(&[-2, 2, 2, 2]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xe);
        assert_eq!(boolean_function.variables_count(), 2);

        let boolean_function =
            SmallBooleanFunction::from_walsh_hadamard_values(&[-8, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xff);
        assert_eq!(boolean_function.variables_count(), 3);

        let boolean_function =
            SmallBooleanFunction::from_walsh_hadamard_values(&[8, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0x00);
        assert_eq!(boolean_function.variables_count(), 3);

        let boolean_function =
            SmallBooleanFunction::from_walsh_hadamard_values(&[0, 0, 0, 0, -4, 4, 4, 4]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0x1e);
        assert_eq!(boolean_function.variables_count(), 3);

        let boolean_function =
            SmallBooleanFunction::from_walsh_hadamard_values(&[0, 0, 0, 0, 4, -4, -4, -4]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xe1);
        assert_eq!(boolean_function.variables_count(), 3);

        let boolean_function =
            SmallBooleanFunction::from_walsh_hadamard_values(&[0, 0, 0, 0, 4, -4, -4]);
        assert!(boolean_function.is_err());
        assert_eq!(
            boolean_function.unwrap_err(),
            crate::BooleanFunctionError::InvalidWalshValuesCount(7)
        );

        let boolean_function = SmallBooleanFunction::from_walsh_hadamard_values(&[0]);
        assert!(boolean_function.is_err());
        assert_eq!(
            boolean_function.unwrap_err(),
            crate::BooleanFunctionError::InvalidWalshValuesCount(1)
        );
    }

    #[test]
    fn test_xor() {
        let mut boolean_function = SmallBooleanFunction::from_truth_table(0x1e, 3).unwrap();
        let boolean_function2 = SmallBooleanFunction::from_truth_table(0xab, 3).unwrap();
        let boolean_function3 = boolean_function ^ boolean_function2;
        boolean_function ^= boolean_function2;
        assert_eq!(boolean_function.get_truth_table_u64(), 0xb5);
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(boolean_function3.get_truth_table_u64(), 0xb5);
        assert_eq!(boolean_function3.variables_count(), 3);
    }

    #[test]
    fn test_biguint_truth_table() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0x1e, 3).unwrap();
        assert_eq!(boolean_function.biguint_truth_table().to_string(), "30");
    }

    #[test]
    fn test_from_walsh_fourier_values() {
        let boolean_function = SmallBooleanFunction::from_walsh_fourier_values(&[
            2, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 0, 2,
        ])
        .unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0x8001);

        let boolean_function =
            SmallBooleanFunction::from_walsh_fourier_values(&[8, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xff);

        let boolean_function =
            SmallBooleanFunction::from_walsh_fourier_values(&[0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0x00);

        let boolean_function =
            SmallBooleanFunction::from_walsh_fourier_values(&[4, -4, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xaa);

        let boolean_function = SmallBooleanFunction::from_walsh_fourier_values(&[
            64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
        ])
        .unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xffffffffffffffff);
        assert_eq!(boolean_function.variables_count(), 6);
    }

    #[test]
    fn test_iter() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0x1e, 3).unwrap();
        let mut iter = boolean_function.iter();
        assert_eq!(iter.next().unwrap(), false);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), false);
        assert_eq!(iter.next().unwrap(), false);
        assert_eq!(iter.next().unwrap(), false);
        assert!(iter.next().is_none());
    }
}
