use crate::anf_polynom::AnfPolynomial;
use crate::{BooleanFunction, BooleanFunctionError, BooleanFunctionImpl, BooleanFunctionType};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::any::Any;
use std::ops::{BitXor, BitXorAssign};
use itertools::{enumerate, Itertools};
use num_integer::binomial;
use crate::boolean_function_error::XOR_DIFFERENT_VAR_COUNT_PANIC_MSG;
use crate::utils::{fast_anf_transform_biguint, left_kernel_boolean};

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BigBooleanFunction {
    variables_count: usize,
    truth_table: BigUint,
}

impl BigBooleanFunction {
    pub fn from_truth_table(truth_table: BigUint, variables_count: usize) -> Self {
        // TODO check too big for variables_count
        BigBooleanFunction {
            variables_count,
            truth_table,
        }
    }

    pub fn derivative_inner(&self, direction: u32) -> Result<Self, BooleanFunctionError> {
        #[cfg(debug_assertions)]
        {
            let max_input_value = self.get_max_input_value();
            if direction > max_input_value {
                return Err(BooleanFunctionError::TooBigDerivativeDirection(
                    max_input_value,
                ));
            }
        }
        let mut derivative_truth_table: BigUint = BigUint::zero();
        for x in 0..=self.get_max_input_value() {
            if self.compute_cellular_automata_rule(x)
                ^ self.compute_cellular_automata_rule(x ^ direction)
            {
                derivative_truth_table |= BigUint::one() << x;
            }
        }
        Ok(BigBooleanFunction {
            variables_count: self.variables_count,
            truth_table: derivative_truth_table,
        })
    }

    pub fn reverse_inner(&self) -> Self {
        BigBooleanFunction {
            variables_count: self.variables_count,
            truth_table: self.truth_table.clone() ^ ((BigUint::one() << (1 << self.variables_count)) - BigUint::one()),
        }
    }

    pub fn annihilator_inner(&self, max_degree: usize) -> Option<(BigBooleanFunction, usize, usize)> {
        if self.truth_table == BigUint::zero() {
            let max_possible_function_tt = (BigUint::one() << (1 << self.variables_count)) - BigUint::one();
            let dim_annihilator_vec_space = (0..=max_degree).map(|i| binomial(self.variables_count as u64, i as u64)).sum::<u64>() as usize;
            return Some((Self::from_truth_table(max_possible_function_tt, self.variables_count), 0, dim_annihilator_vec_space));
        }

        let truth_table_non_zero_positions = (0u32..(1 << self.variables_count))
            .filter(|bit_pos| self.truth_table.clone() & (BigUint::one() << bit_pos) != BigUint::zero())
            .collect::<Vec<u32>>();

        let matrix_out_len = (0..=max_degree).map(|i| binomial(self.variables_count as u64, i as u64)).sum::<u64>() as usize;
        let matrix_in_len = truth_table_non_zero_positions.len();
        let mut matrix: Vec<Vec<bool>> = vec![vec![false; matrix_in_len]; matrix_out_len];

        let mut r = [AnfPolynomial::from_anf_big(&BigUint::one(), self.variables_count)].to_vec();

        for i in 1..=max_degree {
            for comb in (0..self.variables_count).combinations(i) {
                let mut bit_index = 0;
                for monomial in comb {
                    bit_index |= 1 << monomial;
                }
                let anf = BigUint::one() << bit_index;
                r.push(AnfPolynomial::from_anf_big(&anf, self.variables_count));
            }
        }

        for (i, m) in enumerate(r.iter()) {
            let truth_table = fast_anf_transform_biguint(&m.get_polynomial_big(), self.variables_count); // performance bottleneck
            for (j, v) in enumerate(truth_table_non_zero_positions.iter()) {
                matrix[i][j] = truth_table.bit(*v as u64)
            }
        }

        let left_kernel = left_kernel_boolean(&matrix);

        if left_kernel.is_empty() {
            return None;
        }

        let annihilator_anf = enumerate(r.iter())
            .filter(|(i, _)| left_kernel[0][*i])
            .map(|(_, v)| v.get_polynomial_big())
            .fold(BigUint::zero(), |mut sum, val| {sum += val; sum});

        let annihilator_function = Self::from_truth_table(
            fast_anf_transform_biguint(&annihilator_anf, self.variables_count),
            self.variables_count);

        let annihilator_degree = AnfPolynomial::from_anf_big(&annihilator_anf, self.variables_count).get_degree();

        Some((annihilator_function, annihilator_degree, left_kernel.len()))
    }

    /// Reverse walsh transform
    pub fn from_walsh_hadamard_values(walsh_values: &[i32]) -> Result<Self, BooleanFunctionError> {
        let walsh_values_count = walsh_values.len();
        if walsh_values_count < 4 || walsh_values_count.count_ones() != 1 {
            return Err(BooleanFunctionError::InvalidWalshValuesCount(walsh_values_count));
        }
        let num_variables = walsh_values_count.trailing_zeros() as usize;
        let mut truth_table = BigUint::zero();
        for i in 0..(1 << num_variables) {
            let value = walsh_values.iter().enumerate()
                .map(|(w, walsh_value)|walsh_value * ( if (w & i).count_ones() & 1 == 0 {1} else {-1}))
                .sum::<i32>() < 0;
            if value {
                truth_table.set_bit(i as u64, true);
            }
        }
        Ok(Self {
            variables_count: num_variables,
            truth_table,
        })
    }

    pub fn from_walsh_fourier_values(walsh_values: &[i32]) -> Result<Self, BooleanFunctionError> {
        let walsh_values_count = walsh_values.len();
        if walsh_values_count < 4 || walsh_values_count.count_ones() != 1 {
            return Err(BooleanFunctionError::InvalidWalshValuesCount(walsh_values_count));
        }
        let num_variables = walsh_values_count.trailing_zeros() as usize;
        let mut truth_table = BigUint::zero();
        for i in 0..(1 << num_variables) {
            let value = walsh_values.iter().enumerate()
                .map(|(w, walsh_value)|walsh_value * ( if (w & i).count_ones() & 1 == 0 {1} else {-1}))
                .sum::<i32>() != 0;
            if value {
                truth_table.set_bit(i as u64, true);
            }
        }
        Ok(Self {
            variables_count: num_variables,
            truth_table,
        })
    }
}
impl BooleanFunctionImpl for BigBooleanFunction {
    #[inline]
    fn get_num_variables(&self) -> usize {
        self.variables_count
    }

    fn get_boolean_function_type(&self) -> BooleanFunctionType {
        BooleanFunctionType::Big
    }

    fn is_balanced(&self) -> bool {
        let expected_set_number: u64 = 1 << (self.variables_count - 1);
        self.truth_table.count_ones() == expected_set_number
    }

    #[inline]
    fn compute_cellular_automata_rule(&self, input_bits: u32) -> bool {
        #[cfg(debug_assertions)] // TODO use a feature instead, like "unsafe_skip_safety_tests"
        {
            let max_input_value = self.get_max_input_value();
            if input_bits > max_input_value {
                panic!("Input bits must be less or equal than {}", max_input_value);
            }
        }
        (self.truth_table.clone() & (BigUint::one() << input_bits)) != BigUint::zero()
    }

    fn derivative(&self, direction: u32) -> Result<BooleanFunction, BooleanFunctionError> {
        Ok(Box::new(self.derivative_inner(direction)?))
    }

    fn is_linear(&self) -> bool {
        let max_input_value = self.get_max_input_value();
        [&self.truth_table, &self.reverse_inner().truth_table].iter().any(|rule| {
            let mut equivalent_xor_function = BigUint::zero();
            for i in 0..=max_input_value {
                let mut equivalent_xor_function_eval_i = false;
                for j in 0..self.variables_count {
                    if (*rule & (BigUint::one() << (1 << j))) != BigUint::zero() {
                        equivalent_xor_function_eval_i ^= (i & (1 << j)) == 0;
                    }
                }
                equivalent_xor_function |= BigUint::from(equivalent_xor_function_eval_i) << i;
            }
            **rule == equivalent_xor_function || **rule == (Self {
                variables_count: self.variables_count,
                truth_table: equivalent_xor_function,
            }).reverse_inner().truth_table
        })
    }

    fn reverse(&self) -> BooleanFunction {
        Box::new(self.reverse_inner())
    }

    fn algebraic_normal_form(&self) -> AnfPolynomial {
        let anf_form = fast_anf_transform_biguint(&self.truth_table, self.variables_count);
        AnfPolynomial::from_anf_big(&anf_form, self.variables_count)
    }

    fn annihilator(&self, max_degree: usize) -> Option<(BooleanFunction, usize, usize)> {
        let annihilator = self.annihilator_inner(max_degree)?;
        Some((Box::new(annihilator.0), annihilator.1, annihilator.2))
    }

    fn printable_hex_truth_table(&self) -> String {
        format!("{:01$x}", self.truth_table, 1 << (self.variables_count - 2))
    }

    fn biguint_truth_table(&self) -> BigUint {
        self.truth_table.clone()
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_mut_any(&mut self) -> &mut dyn Any {
        self
    }
}

impl BitXorAssign for BigBooleanFunction {
    fn bitxor_assign(&mut self, rhs: Self) {
        if self.variables_count != rhs.variables_count {  // TODO only for debug
            panic!("{}", XOR_DIFFERENT_VAR_COUNT_PANIC_MSG);
        }
        self.truth_table ^= rhs.truth_table;
    }
}

impl BitXor for BigBooleanFunction {
    type Output = Self;

    fn bitxor(mut self, rhs: Self) -> Self::Output {
        self ^= rhs;
        self
    }
}

#[cfg(test)]
mod tests {
    use crate::{BigBooleanFunction, BooleanFunctionImpl};
    use num_bigint::BigUint;
    use num_traits::{Num, One, Zero};

    #[test]
    fn test_get_num_variables() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        assert_eq!(boolean_function.get_num_variables(), 7);
    }

    #[test]
    fn test_is_balanced() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        assert_eq!(boolean_function.is_balanced(), true);

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD1", 16).unwrap(),
            7,
        );
        assert_eq!(boolean_function.is_balanced(), false);
    }

    #[test]
    fn test_compute_cellular_automata_rule() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        assert_eq!(boolean_function.compute_cellular_automata_rule(13), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(62), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(64), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(80), true);
        assert_eq!(boolean_function.compute_cellular_automata_rule(100), true);
    }

    #[test]
    fn test_derivative() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        let derivative = boolean_function.derivative(1).unwrap();
        assert_eq!(derivative.get_num_variables(), 7);
        assert_eq!(
            derivative.printable_hex_truth_table(),
            "cfffc3c00fcf0cfff003f3ccf3f0ff30"
        );

        let derivative = boolean_function.derivative(2).unwrap();
        assert_eq!(
            derivative.printable_hex_truth_table(),
            "afffa5aff0aff50f0ffaf55af5f000a0"
        );

        let derivative = boolean_function.derivative(3).unwrap();
        assert_eq!(
            derivative.printable_hex_truth_table(),
            "9000999fff90f6f0fff609690900ff60"
        );

        let derivative = boolean_function.derivative(128);
        assert!(derivative.is_err());
    }

    #[test]
    fn test_printable_hex_truth_table() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "7969817cc5893ba6ac326e47619f5ad0"
        );

        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::one(), 7);
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "00000000000000000000000000000001"
        );

        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::zero(), 7);
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "00000000000000000000000000000000"
        );
    }

    #[test]
    fn test_algebraic_normal_form() {
        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from(30u32), 3);
        assert_eq!(
            boolean_function
                .algebraic_normal_form()
                .get_polynomial_big(),
            BigUint::from(30u32)
        );

        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::zero(), 3);
        assert_eq!(
            boolean_function
                .algebraic_normal_form()
                .get_polynomial_big(),
            BigUint::zero()
        );

        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from(0xffu32), 3);
        assert_eq!(
            boolean_function
                .algebraic_normal_form()
                .get_polynomial_big(),
            BigUint::one()
        );

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2*x3 + x0*x1*x2*x4 + x0*x1*x3*x4 + x1*x2*x3*x4 + x0*x1*x2*x5 + x0*x2*x3*x5 + x1*x2*x3*x5 + x0*x2*x4*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5 + x0*x1*x2*x6 + x0*x1*x3*x6 + x0*x2*x3*x6 + x1*x2*x3*x6 + x0*x2*x4*x6 + x1*x2*x4*x6 + x0*x3*x4*x6 + x2*x3*x4*x6 + x0*x1*x5*x6 + x0*x2*x5*x6 + x1*x2*x5*x6 + x1*x3*x5*x6 + x2*x3*x5*x6 + x3*x4*x5*x6 + x0*x1*x2 + x0*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 + x0*x1*x5 + x0*x2*x5 + x1*x2*x5 + x1*x3*x5 + x2*x3*x5 + x0*x4*x5 + x2*x4*x5 + x3*x4*x5 + x0*x2*x6 + x1*x2*x6 + x3*x4*x6 + x0*x5*x6 + x2*x5*x6 + x3*x5*x6 + x0*x2 + x0*x3 + x2*x4 + x3*x5 + x0*x6 + x1*x6 + x2*x6 + x3*x6 + x5*x6 + x2 + x4 + x5");

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD1", 16).unwrap(),
            7,
        );
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2*x3*x4*x5*x6 + x0*x1*x2*x3*x4*x5 + x0*x1*x2*x3*x4*x6 + x0*x1*x2*x3*x5*x6 + x0*x1*x2*x4*x5*x6 + x0*x1*x3*x4*x5*x6 + x0*x2*x3*x4*x5*x6 + x1*x2*x3*x4*x5*x6 + x0*x1*x2*x3*x4 + x0*x1*x2*x3*x5 + x0*x1*x2*x4*x5 + x0*x1*x3*x4*x5 + x0*x2*x3*x4*x5 + x1*x2*x3*x4*x5 + x0*x1*x2*x3*x6 + x0*x1*x2*x4*x6 + x0*x1*x3*x4*x6 + x0*x2*x3*x4*x6 + x1*x2*x3*x4*x6 + x0*x1*x2*x5*x6 + x0*x1*x3*x5*x6 + x0*x2*x3*x5*x6 + x1*x2*x3*x5*x6 + x0*x1*x4*x5*x6 + x0*x2*x4*x5*x6 + x1*x2*x4*x5*x6 + x0*x3*x4*x5*x6 + x1*x3*x4*x5*x6 + x2*x3*x4*x5*x6 + x0*x2*x3*x4 + x0*x1*x3*x5 + x0*x1*x4*x5 + x0*x3*x4*x5 + x0*x1*x4*x6 + x1*x3*x4*x6 + x0*x3*x5*x6 + x0*x4*x5*x6 + x1*x4*x5*x6 + x2*x4*x5*x6 + x0*x1*x3 + x1*x2*x3 + x0*x1*x4 + x0*x2*x4 + x0*x3*x4 + x0*x3*x5 + x1*x4*x5 + x0*x1*x6 + x0*x3*x6 + x1*x3*x6 + x2*x3*x6 + x0*x4*x6 + x1*x4*x6 + x2*x4*x6 + x1*x5*x6 + x4*x5*x6 + x0*x1 + x1*x2 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x3*x4 + x0*x5 + x1*x5 + x2*x5 + x4*x5 + x4*x6 + x0 + x1 + x3 + x6 + 1");
    }

    #[test]
    fn test_reverse_inner() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        let reversed_boolean_function = boolean_function.reverse_inner();
        assert_eq!(reversed_boolean_function.printable_hex_truth_table(), "86967e833a76c45953cd91b89e60a52f");
        assert_eq!(reversed_boolean_function.get_num_variables(), 7);
    }

    #[test]
    fn test_is_linear() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        assert!(!boolean_function.is_linear());

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD1", 16).unwrap(),
            7,
        );
        assert!(!boolean_function.is_linear());

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("00000000000000000000000000000000", 16).unwrap(),
            7,
        );
        assert!(boolean_function.is_linear());

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ffffffffffffffffffffffffffffffff", 16).unwrap(),
            7,
        );
        assert!(boolean_function.is_linear());

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("0000000000000000ffffffffffffffff", 16).unwrap(),
            7,
        );
        assert!(boolean_function.is_linear());
    }

    #[test]
    fn test_annihilator_inner() {
        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("00000000000000000000000000000000", 16).unwrap(),
            7,
        );
        let annihilator = boolean_function.annihilator_inner(0).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffffffffffffffffffffffffffff");
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 1);

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ffffffffffffffffffffffffffffffff", 16).unwrap(),
            7,
        );
        let annihilator = boolean_function.annihilator_inner(7);
        assert!(annihilator.is_none());

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        let annihilator = boolean_function.annihilator_inner(2);
        assert!(annihilator.is_none());

        let annihilator = boolean_function.annihilator_inner(3).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "80921c010276c440400810a80e200425");
        assert_eq!(annihilator.1, 3);
        assert_eq!(annihilator.2, 2);

        let annihilator = boolean_function.annihilator_inner(7).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "80921c010276c440400810a80e200425");
        assert_eq!(annihilator.1, 3);
        assert_eq!(annihilator.2, 64);

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("80921c010276c440400810a80e200425", 16).unwrap(),
            7,
        );
        let annihilator = boolean_function.annihilator_inner(1);
        assert!(annihilator.is_none());

        let annihilator = boolean_function.annihilator_inner(2).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "22442244118811882244224411881188");
        assert_eq!(annihilator.1, 2);
        assert_eq!(annihilator.2, 5);

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("0000000000000000ffffffffffffffff", 16).unwrap(),
            7,
        );
        let annihilator = boolean_function.annihilator_inner(0);
        assert!(annihilator.is_none());

        let annihilator = boolean_function.annihilator_inner(1).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffffffffffff0000000000000000");
        assert_eq!(annihilator.1, 1);
        assert_eq!(annihilator.2, 1);

        let annihilator = boolean_function.annihilator_inner(4).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffffffffffff0000000000000000");
        assert_eq!(annihilator.1, 1);
        assert_eq!(annihilator.2, 42);

        let annihilator = boolean_function.annihilator_inner(7).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffffffffffff0000000000000000");
        assert_eq!(annihilator.1, 1);
        assert_eq!(annihilator.2, 64);

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("80921c010276c44224422441188118822442244118811880400810a80e200425", 16).unwrap(),
            8,
        );
        let annihilator = boolean_function.annihilator_inner(0);
        assert!(annihilator.is_none());
        let annihilator = boolean_function.annihilator_inner(1);
        assert!(annihilator.is_none());

        let annihilator = boolean_function.annihilator_inner(2).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "2244224411881188d2b4d2b4e178e178d2b4d2b4e178e1782244224411881188");
        assert_eq!(annihilator.1, 2);
        assert_eq!(annihilator.2, 2);

        let annihilator = boolean_function.annihilator_inner(5).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "2244224411881188d2b4d2b4e178e178d2b4d2b4e178e1782244224411881188");
        assert_eq!(annihilator.1, 2);
        assert_eq!(annihilator.2, 155);

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16).unwrap(),
            8,
        );
        let annihilator = boolean_function.annihilator_inner(4);
        assert!(annihilator.is_none());

        let boolean_function = BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("0000000000000000000000000000000000000000000000000000000000000000", 16).unwrap(),
            8,
        );
        let annihilator = boolean_function.annihilator_inner(4).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 163);
    }

    #[test]
    fn test_from_walsh_hadamard_values() {
        let boolean_function = BigBooleanFunction::from_walsh_hadamard_values(&[-128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "ffffffffffffffffffffffffffffffff");

        let boolean_function = BigBooleanFunction::from_walsh_hadamard_values(&[128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00000000000000000000000000000000");

        let boolean_function = BigBooleanFunction::from_walsh_hadamard_values(&[128, 0, 8, 8, 0, 0, 8, 8, -8, 8, -16, 0, 8, -8, -64, -16, 0, -16, 8, -8, 0, -16, 8, -8, 8, 8, 0, 0, -8, -8, -16, -16, -8, 8, 0, -16, -8, 8, 0, 16, 0, 0, -8, 24, 16, -48, 8, 8, 8, 8, 16, 16, 8, 8, -16, 16, 0, 16, -8, 8, -16, 32, 8, 24, 8, 8, 0, 0, -8, -8, 16, -16, 0, 16, 8, -8, 0, -16, 8, -8, -8, 8, -16, 0, 8, -8, 0, 16, -32, 0, 8, 8, 0, 0, 8, 8, 0, 0, -8, -8, -16, -16, 8, 8, 8, -8, -16, 0, 8, -8, -16, 0, 0, -16, -8, 8, -16, 0, 8, -8, -8, -8, 0, 0, -8, -8, 0, 0, 12, 4, 4, -4, -4, -12, 20, -20, 4, 12, 12, -12, 4, -20, 12, -12, -4, 4, -12, -4, 12, -12, 4, 12, -28, -4, 12, 4, 4, -4, 12, 4, 4, -4, -4, -12, -12, -20, 12, 4, 12, -12, -12, -4, 12, -12, -12, -4, 4, -20, -4, 4, -12, -4, 12, -12, -4, -12, 4, -4, -4, -12, 4, -4, -4, 4, 4, 12, -4, 4, 4, 12, -12, 12, -20, 4, 4, -4, 60, -12, -4, -12, 4, -4, -4, -12, 4, -4, 4, 12, -4, 4, -12, -4, -20, -12, -12, -20, -4, 84, -12, -20, -4, -12, -4, -28, -12, -4, 12, 52, 4, -20, 4, -20, 12, -12, 4, -20, -20, -12, -4, -12, -12, -20, -20, 4, 4, -4]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "80921c010276c44224422441188118822442244118811880400810a80e200425");

        let boolean_function = BigBooleanFunction::from_walsh_hadamard_values(&[64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -64, 0, 0, 0, 0, 0, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "22442244118811882244224411881188");

        let boolean_function = BigBooleanFunction::from_walsh_hadamard_values(&[0, 0, 0, 0, 4, -4, -4]);
        assert!(boolean_function.is_err());
        assert_eq!(boolean_function.unwrap_err(), crate::BooleanFunctionError::InvalidWalshValuesCount(7));

        let boolean_function = BigBooleanFunction::from_walsh_hadamard_values(&[0]);
        assert!(boolean_function.is_err());
        assert_eq!(boolean_function.unwrap_err(), crate::BooleanFunctionError::InvalidWalshValuesCount(1));
    }

    #[test]
    fn test_xor() {
        let mut boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("80921c010276c440400810a80e200425", 16).unwrap(), 7);
        let boolean_function2 = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("22442244118811882244224411881188", 16).unwrap(), 7);
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2;
        assert_eq!(boolean_function.printable_hex_truth_table(), "a2d63e4513fed5c8624c32ec1fa815ad");
        assert_eq!(boolean_function.get_num_variables(), 7);
        assert_eq!(boolean_function3.printable_hex_truth_table(), "a2d63e4513fed5c8624c32ec1fa815ad");
        assert_eq!(boolean_function3.get_num_variables(), 7);
    }

    #[test]
    fn test_biguint_truth_table() {
        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("80921c010276c440400810a80e200425", 16).unwrap(), 7);
        assert_eq!(boolean_function.biguint_truth_table(), BigUint::from_str_radix("80921c010276c440400810a80e200425", 16).unwrap());
    }

    #[test]
    fn test_from_walsh_fourier_values() {
        let boolean_function = BigBooleanFunction::from_walsh_fourier_values(&[2, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 0, 2]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "8001");

        let boolean_function = BigBooleanFunction::from_walsh_fourier_values(&[8, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "ff");

        let boolean_function = BigBooleanFunction::from_walsh_fourier_values(&[0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00");

        let boolean_function = BigBooleanFunction::from_walsh_fourier_values(&[4, -4, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "aa");

        let boolean_function = BigBooleanFunction::from_walsh_fourier_values(&[64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "ffffffffffffffff");
        assert_eq!(boolean_function.get_num_variables(), 6);
    }
}
