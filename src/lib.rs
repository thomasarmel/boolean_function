mod anf_polynom;
mod big_boolean_function;
mod boolean_function_error;
mod small_boolean_function;
mod utils;

use crate::anf_polynom::AnfPolynomial;
use crate::BooleanFunctionError::{
    StringHexParseError, UnexpectedError, WrongStringHexTruthTableLength,
};
pub use big_boolean_function::BigBooleanFunction;
pub use boolean_function_error::BooleanFunctionError;
use num_bigint::BigUint;
use num_traits::{Num, ToPrimitive};
pub use small_boolean_function::SmallBooleanFunction;
use std::any::Any;
use std::collections::HashMap;
use std::fmt::Debug;
use std::ops::{BitXor, BitXorAssign};
use crate::boolean_function_error::XOR_DIFFERENT_VAR_COUNT_PANIC_MSG;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanFunctionType {
    Small,
    Big,
}

pub trait BooleanFunctionImpl: Debug + Any {
    fn get_num_variables(&self) -> usize;

    fn get_boolean_function_type(&self) -> BooleanFunctionType;
    #[inline]
    fn get_max_input_value(&self) -> u32 {
        (1 << self.get_num_variables()) - 1
    }
    fn is_balanced(&self) -> bool;

    /// debug_assertions
    fn compute_cellular_automata_rule(&self, input_bits: u32) -> bool;
    fn walsh_hadamard_transform(&self, w: u32) -> i32 {
        let max_input_value = self.get_max_input_value();
        #[cfg(debug_assertions)] // TODO
        if w > max_input_value {
            panic!(
                "Too big Walsh parameter direction, must be <= {}",
                max_input_value
            );
        }
        (0..=max_input_value)
            .map(|x| {
                if (self.compute_cellular_automata_rule(x) as u32
                    + utils::fast_binary_dot_product(w, x))
                    & 1
                    == 0
                {
                    // % modulo 2
                    1
                } else {
                    -1
                }
            })
            .sum()
    }

    fn walsh_hadamard_values(&self) -> Vec<i32> {
        (0..=self.get_max_input_value())
            .map(|w| self.walsh_hadamard_transform(w))
            .collect()
    }

    fn absolute_walsh_hadamard_spectrum(&self) -> HashMap<u32, usize> {
        let mut absolute_walsh_value_count_map: HashMap<u32, usize> = HashMap::new();
        (0..=self.get_max_input_value()).for_each(|w| {
            let absolute_walsh_value = self.walsh_hadamard_transform(w).unsigned_abs();
            if !absolute_walsh_value_count_map.contains_key(&absolute_walsh_value) {
                absolute_walsh_value_count_map.insert(absolute_walsh_value, 1);
            } else {
                let count = absolute_walsh_value_count_map
                    .get_mut(&absolute_walsh_value)
                    .unwrap();
                *count += 1;
            }
        });
        absolute_walsh_value_count_map
    }

    fn walsh_fourier_transform(&self, w: u32) -> i32 {
        let max_input_value = self.get_max_input_value();
        #[cfg(debug_assertions)] // TODO
        if w > max_input_value {
            panic!(
                "Too big Walsh parameter direction, must be <= {}",
                max_input_value
            );
        }
        (0..=max_input_value)
            .map(|x| {
                if !self.compute_cellular_automata_rule(x) {
                    0
                } else {
                    if utils::fast_binary_dot_product(w, x) & 1 == 0 {
                        1
                    } else {
                        -1
                    }
                }
            })
            .sum()
    }

    fn walsh_fourier_values(&self) -> Vec<i32> {
        (0..=self.get_max_input_value())
            .map(|w| self.walsh_fourier_transform(w))
            .collect()
    }

    fn auto_correlation_transform(&self, w: u32) -> i32 {
        let max_input_value = self.get_max_input_value();
        #[cfg(debug_assertions)] // TODO
        if w > max_input_value {
            panic!(
                "Too big Walsh parameter direction, must be <= {}",
                max_input_value
            );
        }
        (0..=max_input_value)
            .map(|x| {
                if self.compute_cellular_automata_rule(x)
                    ^ self.compute_cellular_automata_rule(x ^ w)
                {
                    -1
                } else {
                    1
                }
            })
            .sum()
    }

    fn absolute_autocorrelation(&self) -> HashMap<u32, usize> {
        let mut absolute_autocorrelation_value_count_map: HashMap<u32, usize> = HashMap::new();
        (0..=self.get_max_input_value()).for_each(|w| {
            let absolute_autocorrelation_value = self.auto_correlation_transform(w).unsigned_abs();
            if !absolute_autocorrelation_value_count_map
                .contains_key(&absolute_autocorrelation_value)
            {
                absolute_autocorrelation_value_count_map.insert(absolute_autocorrelation_value, 1);
            } else {
                let count = absolute_autocorrelation_value_count_map
                    .get_mut(&absolute_autocorrelation_value)
                    .unwrap();
                *count += 1;
            }
        });
        absolute_autocorrelation_value_count_map
    }

    /// https://www.researchgate.net/publication/322383819_Distribution_of_the_absolute_indicator_of_random_Boolean_functions
    /// Max autocorrelation from w=1 to max_input_value
    fn absolute_indicator(&self) -> u32 {
        (1..=self.get_max_input_value())
            .map(|w| self.auto_correlation_transform(w).unsigned_abs())
            .max().unwrap_or(0)
    }

    fn derivative(&self, direction: u32) -> Result<BooleanFunction, BooleanFunctionError>;

    fn is_linear(&self) -> bool;

    fn reverse(&self) -> BooleanFunction;

    fn algebraic_normal_form(&self) -> AnfPolynomial;

    fn algebraic_degree(&self) -> usize {
        self.algebraic_normal_form().get_degree()
    }

    fn is_symmetric(&self) -> bool {
        let variables_count = self.get_num_variables() as u32;
        let precomputed_hamming_weights = (0..(variables_count + 1))
            .map(|i| self.compute_cellular_automata_rule(u32::MAX.checked_shr(32 - i).unwrap_or(0)))
            .collect::<Vec<bool>>();

        (0u32..(1 << variables_count)).all(|i| {
            precomputed_hamming_weights[i.count_ones() as usize] == self.compute_cellular_automata_rule(i)
        })
    }

    fn nonlinearity(&self) -> u32 {
        ((1 << self.get_num_variables()) - (0..=self.get_max_input_value()).map(|x| {
            self.walsh_hadamard_transform(x).unsigned_abs()
        }).max().unwrap_or(0)) >> 1
    }

    /// Meaning maximally nonlinear
    /// All derivative are balanced
    fn is_bent(&self) -> bool {
        if self.get_num_variables() & 1 != 0 {
            return false;
        }
        self.nonlinearity()
            == ((1 << self.get_num_variables())
                - (1 << (self.get_num_variables() >> 1)))
                >> 1
    }

    /// Function, degree and dimension of annihilator vector space
    /// Special case: annihilator of zero function is one function, by convention
    /// $ f(x).g(x) = 0 \ \ \forall x \in \mathbb{F}_2^n. $
    fn annihilator(&self, max_degree: usize) -> Option<(BooleanFunction, usize, usize)>; // TODO max degree in Option

    fn algebraic_immunity(&self) -> usize {
        match self.annihilator(self.get_num_variables()) {
            None => 0,
            Some(annihilator) => annihilator.1
        }
    }

    fn is_plateaued(&self) -> bool {
        let absolute_walsh_hadamard_spectrum = self.absolute_walsh_hadamard_spectrum();
        absolute_walsh_hadamard_spectrum.len() == 1 || (absolute_walsh_hadamard_spectrum.len() == 2 && absolute_walsh_hadamard_spectrum.contains_key(&0))
    }

    fn sum_of_square_indicator(&self) -> usize {
        (0..=self.get_max_input_value())
            .map(|w| self.auto_correlation_transform(w))
            .map(|value| (value as i64 * value as i64) as usize)
            .sum()
    }

    fn has_linear_structure(&self) -> bool {
        (1..=self.get_max_input_value())
            .any(|x| self.auto_correlation_transform(x).unsigned_abs() == 1 << self.get_num_variables())
    }

    /// Will panic if value > max_input_value
    fn is_linear_structure(&self, value: u32) -> bool {
        #[cfg(debug_assertions)] // TODO
        {
            let max_input_value = self.get_max_input_value();
            if value > max_input_value {
                panic!(
                    "Too big value parameter, must be <= {}",
                    max_input_value
                );
            }
        }
        self.auto_correlation_transform(value).unsigned_abs() == 1 << self.get_num_variables()
    }

    /// https://www.sciencedirect.com/topics/mathematics/linear-structure
    fn linear_structures(&self) -> Vec<u32> {
        (0..=self.get_max_input_value())
            .filter(|x| self.auto_correlation_transform(*x).unsigned_abs() == 1 << self.get_num_variables())
            .collect()
    }

    /// Xiao Massey theorem
    fn correlation_immunity(&self) -> usize {
        (1..=self.get_max_input_value())
            .filter(|x| self.walsh_hadamard_transform(*x) != 0)
            .map(|x| x.count_ones() as usize)
            .min().unwrap_or(self.get_num_variables() + 1) - 1
    }

    fn resiliency_order(&self) -> Option<usize> {
        if !self.is_balanced() {
            return None;
        }
        Some(self.correlation_immunity())
    }

    fn printable_hex_truth_table(&self) -> String;

    fn biguint_truth_table(&self) -> BigUint;

    fn try_u64_truth_table(&self) -> Option<u64> {
        self.biguint_truth_table().to_u64()
    }

    fn as_any(&self) -> &dyn Any;

    fn as_mut_any(&mut self) -> &mut dyn Any;

    // TODO almost bent, mul (and tt), sum, impl not, iterate on values
}

pub type BooleanFunction = Box<dyn BooleanFunctionImpl + Send + Sync>;

impl Clone for BooleanFunction {
    fn clone(&self) -> Self {
        match self.get_boolean_function_type() {
            BooleanFunctionType::Small => Box::new(self.as_any().downcast_ref::<SmallBooleanFunction>().unwrap().clone()),
            BooleanFunctionType::Big => Box::new(self.as_any().downcast_ref::<BigBooleanFunction>().unwrap().clone()),
        }
    }
}

impl PartialEq for BooleanFunction {
    fn eq(&self, other: &Self) -> bool {
        if self.get_boolean_function_type() != other.get_boolean_function_type() {
            // TODO compare type id
            return false;
        }
        match self.get_boolean_function_type() {
            BooleanFunctionType::Small => {
                let self_small_boolean_function = self
                    .as_any()
                    .downcast_ref::<SmallBooleanFunction>()
                    .unwrap();
                let other_small_boolean_function = other
                    .as_any()
                    .downcast_ref::<SmallBooleanFunction>()
                    .unwrap();
                self_small_boolean_function == other_small_boolean_function
            }
            BooleanFunctionType::Big => {
                let self_big_boolean_function =
                    self.as_any().downcast_ref::<BigBooleanFunction>().unwrap();
                let other_big_boolean_function =
                    other.as_any().downcast_ref::<BigBooleanFunction>().unwrap();
                self_big_boolean_function == other_big_boolean_function
            }
        }
    }
}

impl Eq for BooleanFunction {}

impl BitXorAssign for BooleanFunction {
    fn bitxor_assign(&mut self, rhs: Self) {
        let self_num_variables = self.get_num_variables();
        let rhs_num_variables = rhs.get_num_variables();
        if self_num_variables != rhs_num_variables {
            panic!("{}", XOR_DIFFERENT_VAR_COUNT_PANIC_MSG);
        }
        match (self.get_boolean_function_type(), rhs.get_boolean_function_type()) {
            (BooleanFunctionType::Small, BooleanFunctionType::Small) => {
                let self_small_boolean_function = self.as_mut_any().downcast_mut::<SmallBooleanFunction>().unwrap();
                let rhs_small_boolean_function = rhs.as_any().downcast_ref::<SmallBooleanFunction>().unwrap();
                *self_small_boolean_function ^= *rhs_small_boolean_function;
            },
            (BooleanFunctionType::Small, BooleanFunctionType::Big) => {
                let self_small_boolean_function = self.as_mut_any().downcast_mut::<SmallBooleanFunction>().unwrap();
                let rhs_small_boolean_function = SmallBooleanFunction::from_truth_table(
                    rhs.as_any().downcast_ref::<BigBooleanFunction>().unwrap().biguint_truth_table().to_u64().unwrap(),
                    self_num_variables
                ).unwrap();
                *self_small_boolean_function ^= rhs_small_boolean_function;
            }
            (BooleanFunctionType::Big, BooleanFunctionType::Small) => {
                let self_big_boolean_function = self.as_mut_any().downcast_mut::<BigBooleanFunction>().unwrap();
                let rhs_small_boolean_function = BigBooleanFunction::from_truth_table(
                    rhs.as_any().downcast_ref::<SmallBooleanFunction>().unwrap().biguint_truth_table(),
                    rhs_num_variables
                );
                *self_big_boolean_function ^= rhs_small_boolean_function;
            }
            (BooleanFunctionType::Big, BooleanFunctionType::Big) => {
                let self_big_boolean_function = self.as_mut_any().downcast_mut::<BigBooleanFunction>().unwrap();
                let rhs_big_boolean_function = rhs.as_any().downcast_ref::<BigBooleanFunction>().unwrap();
                *self_big_boolean_function ^= rhs_big_boolean_function.clone();
            }
        }
    }
}

impl BitXor for BooleanFunction {
    type Output = Self;

    fn bitxor(mut self, rhs: Self) -> Self::Output {
        self ^= rhs;
        self
    }
}

pub fn boolean_function_from_hex_string_truth_table(
    hex_truth_table: &str,
) -> Result<BooleanFunction, BooleanFunctionError> {
    if hex_truth_table.len().count_ones() != 1 {
        return Err(WrongStringHexTruthTableLength);
    }
    let num_variables = (hex_truth_table.len() << 2).trailing_zeros() as usize;
    if num_variables <= 6 {
        Ok(Box::new(
            SmallBooleanFunction::from_truth_table(
                u64::from_str_radix(hex_truth_table, 16).map_err(|_| StringHexParseError)?,
                num_variables,
            )
            .map_err(|_| UnexpectedError)?,
        ))
    } else {
        Ok(Box::new(BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix(hex_truth_table, 16).map_err(|_| StringHexParseError)?,
            num_variables,
        )))
    }
}

pub fn boolean_function_from_reverse_walsh_hadamard_transform(walsh_hadamard_values: &[i32]) -> Result<BooleanFunction, BooleanFunctionError> {
    const MAX_WALSH_VALUES_SMALL: usize = 64; // (2^6)
    if walsh_hadamard_values.len() > MAX_WALSH_VALUES_SMALL {
        return Ok(Box::new(BigBooleanFunction::from_walsh_hadamard_values(walsh_hadamard_values)?)); // Error is handled in BigBooleanFunction constructor
    }
    Ok(Box::new(SmallBooleanFunction::from_walsh_hadamard_values(walsh_hadamard_values)?))
}

pub fn boolean_function_from_reverse_walsh_fourier_transform(walsh_hadamard_values: &[i32]) -> Result<BooleanFunction, BooleanFunctionError> {
    const MAX_WALSH_VALUES_SMALL: usize = 64; // (2^6)
    if walsh_hadamard_values.len() > MAX_WALSH_VALUES_SMALL {
        return Ok(Box::new(BigBooleanFunction::from_walsh_fourier_values(walsh_hadamard_values)?)); // Error is handled in BigBooleanFunction constructor
    }
    Ok(Box::new(SmallBooleanFunction::from_walsh_fourier_values(walsh_hadamard_values)?))
}

// TODO from polynomial etc.

#[cfg(test)]
mod tests {
    use crate::{BooleanFunction, BooleanFunctionType};
    use std::collections::HashMap;
    use num_bigint::BigUint;
    use num_traits::Num;
    use crate::BooleanFunctionError::InvalidWalshValuesCount;

    #[test]
    fn test_boolean_function_from_hex_string_truth_table() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD");
        assert!(boolean_function.is_err());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("");
        assert!(boolean_function.is_err());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("fe1z");
        assert!(boolean_function.is_err());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.get_num_variables(), 4);

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.get_num_variables(), 7);
    }

    #[test]
    fn test_get_num_variables() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.get_num_variables(), 7);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.get_num_variables(), 4);
    }

    #[test]
    fn test_get_boolean_function_type() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6").unwrap();
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
    }

    #[test]
    fn test_get_max_input_value() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.get_max_input_value(), 127);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.get_max_input_value(), 15);
    }

    #[test]
    fn test_is_balanced() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert!(boolean_function.is_balanced());

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD1")
                .unwrap();
        assert!(!boolean_function.is_balanced());

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("aa55aa55").unwrap();
        assert!(boolean_function.is_balanced());

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("abce1234").unwrap();
        assert!(!boolean_function.is_balanced());
    }

    #[test]
    fn test_compute_cellular_automata_rule() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(boolean_function.compute_cellular_automata_rule(0), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(1), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(4), true);
        assert_eq!(boolean_function.compute_cellular_automata_rule(8), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(23), true);

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.compute_cellular_automata_rule(13), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(62), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(64), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(80), true);
        assert_eq!(boolean_function.compute_cellular_automata_rule(100), true);
    }

    #[test]
    fn test_walsh_hadamard_transform() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.walsh_hadamard_transform(0), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(1), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(7), -16);
        assert_eq!(boolean_function.walsh_hadamard_transform(15), 16);
        assert_eq!(boolean_function.walsh_hadamard_transform(126), 16);
        assert_eq!(boolean_function.walsh_hadamard_transform(127), -16);

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("aa55aa55").unwrap();
        assert_eq!(boolean_function.walsh_hadamard_transform(0), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(1), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(9), -32);
        assert_eq!(boolean_function.walsh_hadamard_transform(31), 0);

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(boolean_function.walsh_hadamard_transform(0), 2);
        assert_eq!(boolean_function.walsh_hadamard_transform(1), 6);
        assert_eq!(boolean_function.walsh_hadamard_transform(2), -2);
        assert_eq!(boolean_function.walsh_hadamard_transform(31), -6);
    }

    #[test]
    fn test_absolute_walsh_hadamard_spectrum() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.absolute_walsh_hadamard_spectrum(),
            HashMap::from([(0, 64), (16, 64)])
        );

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(
            boolean_function.absolute_walsh_hadamard_spectrum(),
            HashMap::from([(6, 10), (10, 6), (2, 16)])
        );
    }

    #[test]
    fn test_auto_correlation_transform() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.auto_correlation_transform(0), 128);
        assert_eq!(boolean_function.auto_correlation_transform(1), -24);
        assert_eq!(boolean_function.auto_correlation_transform(126), -8);
        assert_eq!(boolean_function.auto_correlation_transform(127), -32);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("03").unwrap();
        assert_eq!(boolean_function.auto_correlation_transform(0), 8);
        assert_eq!(boolean_function.auto_correlation_transform(1), 8);
        assert_eq!(boolean_function.auto_correlation_transform(2), 0);
        assert_eq!(boolean_function.auto_correlation_transform(3), 0);
        assert_eq!(boolean_function.auto_correlation_transform(4), 0);
        assert_eq!(boolean_function.auto_correlation_transform(5), 0);
        assert_eq!(boolean_function.auto_correlation_transform(6), 0);
        assert_eq!(boolean_function.auto_correlation_transform(7), 0);
    }

    #[test]
    fn test_absolute_autocorrelation_spectrum() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.absolute_autocorrelation(),
            HashMap::from([(0, 33), (8, 58), (16, 28), (24, 6), (32, 2), (128, 1)])
        );

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(
            boolean_function.absolute_autocorrelation(),
            HashMap::from([(4, 25), (12, 6), (32, 1)])
        );
    }

    #[test]
    fn test_derivative() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("aa55aa55").unwrap();
        let derivative = boolean_function.derivative(1).unwrap();
        assert_eq!(derivative.get_num_variables(), 5);
        assert_eq!(derivative.printable_hex_truth_table(), "ffffffff");

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let derivative = boolean_function.derivative(1).unwrap();
        assert_eq!(derivative.get_num_variables(), 7);
        assert_eq!(
            derivative.printable_hex_truth_table(),
            "cfffc3c00fcf0cfff003f3ccf3f0ff30"
        );
    }

    #[test]
    fn test_algebraic_normal_form() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2*x3 + x0*x1*x2*x4 + x0*x1*x3*x4 + x1*x2*x3*x4 + x0*x1*x2*x5 + x0*x2*x3*x5 + x1*x2*x3*x5 + x0*x2*x4*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5 + x0*x1*x2*x6 + x0*x1*x3*x6 + x0*x2*x3*x6 + x1*x2*x3*x6 + x0*x2*x4*x6 + x1*x2*x4*x6 + x0*x3*x4*x6 + x2*x3*x4*x6 + x0*x1*x5*x6 + x0*x2*x5*x6 + x1*x2*x5*x6 + x1*x3*x5*x6 + x2*x3*x5*x6 + x3*x4*x5*x6 + x0*x1*x2 + x0*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 + x0*x1*x5 + x0*x2*x5 + x1*x2*x5 + x1*x3*x5 + x2*x3*x5 + x0*x4*x5 + x2*x4*x5 + x3*x4*x5 + x0*x2*x6 + x1*x2*x6 + x3*x4*x6 + x0*x5*x6 + x2*x5*x6 + x3*x5*x6 + x0*x2 + x0*x3 + x2*x4 + x3*x5 + x0*x6 + x1*x6 + x2*x6 + x3*x6 + x5*x6 + x2 + x4 + x5");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD1").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2*x3*x4*x5*x6 + x0*x1*x2*x3*x4*x5 + x0*x1*x2*x3*x4*x6 + x0*x1*x2*x3*x5*x6 + x0*x1*x2*x4*x5*x6 + x0*x1*x3*x4*x5*x6 + x0*x2*x3*x4*x5*x6 + x1*x2*x3*x4*x5*x6 + x0*x1*x2*x3*x4 + x0*x1*x2*x3*x5 + x0*x1*x2*x4*x5 + x0*x1*x3*x4*x5 + x0*x2*x3*x4*x5 + x1*x2*x3*x4*x5 + x0*x1*x2*x3*x6 + x0*x1*x2*x4*x6 + x0*x1*x3*x4*x6 + x0*x2*x3*x4*x6 + x1*x2*x3*x4*x6 + x0*x1*x2*x5*x6 + x0*x1*x3*x5*x6 + x0*x2*x3*x5*x6 + x1*x2*x3*x5*x6 + x0*x1*x4*x5*x6 + x0*x2*x4*x5*x6 + x1*x2*x4*x5*x6 + x0*x3*x4*x5*x6 + x1*x3*x4*x5*x6 + x2*x3*x4*x5*x6 + x0*x2*x3*x4 + x0*x1*x3*x5 + x0*x1*x4*x5 + x0*x3*x4*x5 + x0*x1*x4*x6 + x1*x3*x4*x6 + x0*x3*x5*x6 + x0*x4*x5*x6 + x1*x4*x5*x6 + x2*x4*x5*x6 + x0*x1*x3 + x1*x2*x3 + x0*x1*x4 + x0*x2*x4 + x0*x3*x4 + x0*x3*x5 + x1*x4*x5 + x0*x1*x6 + x0*x3*x6 + x1*x3*x6 + x2*x3*x6 + x0*x4*x6 + x1*x4*x6 + x2*x4*x6 + x1*x5*x6 + x4*x5*x6 + x0*x1 + x1*x2 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x3*x4 + x0*x5 + x1*x5 + x2*x5 + x4*x5 + x4*x6 + x0 + x1 + x3 + x6 + 1");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1 + x0 + x1 + x2");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "0");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ff").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "1");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("6e").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2 + x0*x1 + x0 + x1");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7b").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1 + x1*x2 + x1 + 1");
    }

    #[test]
    fn test_algebraic_degree() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 4);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD1").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 7);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000000000000000000000000000").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffffffffffffffffffffffffffff").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 2);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0f").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 1);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ff").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("6e").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 3);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7b").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 2);
    }

    #[test]
    fn test_printable_hex_truth_table() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("0069817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "0069817cc5893ba6ac326e47619f5ad0"
        );

        let boolean_function = super::boolean_function_from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "fe12");
    }

    #[test]
    fn test_clone() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let cloned_boolean_function = boolean_function.clone();
        assert_eq!(&boolean_function, &cloned_boolean_function);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("fe12").unwrap();
        let cloned_boolean_function = boolean_function.clone();
        assert_eq!(&boolean_function, &cloned_boolean_function);
    }

    #[test]
    fn test_eq() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(&boolean_function, &boolean_function);

        let boolean_function2 =
            super::boolean_function_from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(&boolean_function2, &boolean_function2);

        assert_ne!(&boolean_function2, &boolean_function);

        let boolean_function3 =
            super::boolean_function_from_hex_string_truth_table("0000fe12").unwrap();
        assert_ne!(&boolean_function3, &boolean_function2);

        let boolean_function4 =
            super::boolean_function_from_hex_string_truth_table("0969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_ne!(&boolean_function, &boolean_function4);
    }

    #[test]
    fn test_reverse() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let reversed_boolean_function = boolean_function.reverse();
        assert_eq!(reversed_boolean_function.get_num_variables(), 7);
        assert_eq!(
            reversed_boolean_function.printable_hex_truth_table(),
            "86967e833a76c45953cd91b89e60a52f"
        );

        let boolean_function = super::boolean_function_from_hex_string_truth_table("fe12").unwrap();
        let reversed_boolean_function = boolean_function.reverse();
        assert_eq!(reversed_boolean_function.get_num_variables(), 4);
        assert_eq!(reversed_boolean_function.printable_hex_truth_table(), "01ed");
    }

    #[test]
    fn test_is_linear() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert!(!boolean_function.is_linear());

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("0000000000000000ffffffffffffffff")
                .unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("abcdef0123456789")
                .unwrap();
        assert!(!boolean_function.is_linear());

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("00000000ffffffff")
                .unwrap();
        assert!(boolean_function.is_linear());
    }

    #[test]
    fn test_is_symmetric() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("00").unwrap();
        assert!(boolean_function.is_symmetric());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ff").unwrap();
        assert!(boolean_function.is_symmetric());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("80").unwrap();
        assert!(boolean_function.is_symmetric());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert!(!boolean_function.is_symmetric());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000008").unwrap();
        assert!(!boolean_function.is_symmetric());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffffffffffffffffffffffffffff").unwrap();
        assert!(boolean_function.is_symmetric());
    }

    #[test]
    fn test_nonlinearity() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0000000a").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 2);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 28);
    }

    #[test]
    fn test_is_bent() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000").unwrap();
        assert!(!boolean_function.is_bent());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert!(boolean_function.is_bent());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000ffffffff").unwrap();
        assert!(!boolean_function.is_bent());
    }

    #[test]
    fn test_annihilator() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000").unwrap();
        let annihilator = boolean_function.annihilator(0).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffff");
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 1);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("abcdef0123456789").unwrap();
        let annihilator = boolean_function.annihilator(4).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "1010101010101010");
        assert_eq!(annihilator.1, 3);
        assert_eq!(annihilator.2, 25);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff").unwrap();
        let annihilator = boolean_function.annihilator(4);
        assert!(annihilator.is_none());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0000000000000000000000000000000000000000000000000000000000000000").unwrap();
        let annihilator = boolean_function.annihilator(4).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 163);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("80921c010276c44224422441188118822442244118811880400810a80e200425").unwrap();
        let annihilator = boolean_function.annihilator(1);
        assert!(annihilator.is_none());
        let annihilator = boolean_function.annihilator(5).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "2244224411881188d2b4d2b4e178e178d2b4d2b4e178e1782244224411881188");
        assert_eq!(annihilator.1, 2);
        assert_eq!(annihilator.2, 155);
    }

    #[test]
    fn test_algebraic_immunity() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0000000000000000000000000000000000000000000000000000000000000000").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0000").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("80921c010276c44224422441188118822442244118811880400810a80e200425").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("2244224411881188d2b4d2b4e178e178d2b4d2b4e178e1782244224411881188").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 3);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);
    }

    #[test]
    fn test_is_plateaued() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff").unwrap();
        assert!(boolean_function.is_plateaued());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("8778").unwrap();
        assert!(boolean_function.is_plateaued());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert!(boolean_function.is_plateaued());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("abcdef1234567890").unwrap();
        assert!(!boolean_function.is_plateaued());
    }

    #[test]
    fn test_sum_of_square_indicator() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 32768);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0000").unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 4096);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("abcdef1234567890abcdef1234567890").unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 84992);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 32768);
    }

    #[test]
    fn test_absolute_indicator() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 8);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 32);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0000").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 16);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("abcdef1234567890abcdef1234567890").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 128);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 32);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 32);
    }

    #[test]
    fn test_linear_structures() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0, 4]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("abcdef1234567890abcdef1234567890").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0, 64]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000ffffffff").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]);
    }

    #[test]
    fn test_has_linear_structure() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("659a").unwrap();
        assert!(boolean_function.has_linear_structure());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("dd0e").unwrap();
        assert!(!boolean_function.has_linear_structure());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000").unwrap();
        assert!(boolean_function.has_linear_structure());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffff").unwrap();
        assert!(boolean_function.has_linear_structure());

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert!(!boolean_function.has_linear_structure());
    }

    #[test]
    fn test_is_linear_structure() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("659a").unwrap();
        assert!(boolean_function.is_linear_structure(1));
        assert!(!boolean_function.is_linear_structure(7));
        assert!(boolean_function.is_linear_structure(9));
    }

    #[test]
    fn test_walsh_hadamard_values() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("dd0e").unwrap();
        assert_eq!(boolean_function.walsh_hadamard_values(), [-2, -2, 6, -2, -6, 2, 2, 2, 6, 6, -2, 6, -6, 2, 2, 2]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0000").unwrap();
        assert_eq!(boolean_function.walsh_hadamard_values(), [16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(boolean_function.walsh_hadamard_values(), [-16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_boolean_function_from_reverse_walsh_transform() {
        let boolean_function = super::boolean_function_from_reverse_walsh_hadamard_transform(&[-2, -2, 6, -2, -6, 2, 2, 2, 6, 6, -2, 6, -6, 2, 2, 2]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "dd0e");
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Small);

        let boolean_function = super::boolean_function_from_reverse_walsh_hadamard_transform(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "ffffffff00000000");
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Small);

        let boolean_function = super::boolean_function_from_reverse_walsh_hadamard_transform(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "ffffffffffffffff0000000000000000");
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Big);

        let boolean_function = super::boolean_function_from_reverse_walsh_hadamard_transform(&[16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "0000");

        let boolean_function = super::boolean_function_from_reverse_walsh_hadamard_transform(&[128, 0, 8, 8, 0, 0, 8, 8, -8, 8, -16, 0, 8, -8, -64, -16, 0, -16, 8, -8, 0, -16, 8, -8, 8, 8, 0, 0, -8, -8, -16, -16, -8, 8, 0, -16, -8, 8, 0, 16, 0, 0, -8, 24, 16, -48, 8, 8, 8, 8, 16, 16, 8, 8, -16, 16, 0, 16, -8, 8, -16, 32, 8, 24, 8, 8, 0, 0, -8, -8, 16, -16, 0, 16, 8, -8, 0, -16, 8, -8, -8, 8, -16, 0, 8, -8, 0, 16, -32, 0, 8, 8, 0, 0, 8, 8, 0, 0, -8, -8, -16, -16, 8, 8, 8, -8, -16, 0, 8, -8, -16, 0, 0, -16, -8, 8, -16, 0, 8, -8, -8, -8, 0, 0, -8, -8, 0, 0, 12, 4, 4, -4, -4, -12, 20, -20, 4, 12, 12, -12, 4, -20, 12, -12, -4, 4, -12, -4, 12, -12, 4, 12, -28, -4, 12, 4, 4, -4, 12, 4, 4, -4, -4, -12, -12, -20, 12, 4, 12, -12, -12, -4, 12, -12, -12, -4, 4, -20, -4, 4, -12, -4, 12, -12, -4, -12, 4, -4, -4, -12, 4, -4, -4, 4, 4, 12, -4, 4, 4, 12, -12, 12, -20, 4, 4, -4, 60, -12, -4, -12, 4, -4, -4, -12, 4, -4, 4, 12, -4, 4, -12, -4, -20, -12, -12, -20, -4, 84, -12, -20, -4, -12, -4, -28, -12, -4, 12, 52, 4, -20, 4, -20, 12, -12, 4, -20, -20, -12, -4, -12, -12, -20, -20, 4, 4, -4]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "80921c010276c44224422441188118822442244118811880400810a80e200425");
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Big);

        let boolean_function = super::boolean_function_from_reverse_walsh_hadamard_transform(&[128, 0, 8, 8, 0, 0, 8, 8, -8, 8, -16, 0, 8, -8, -64, -16, 0, -16, 8, -8, 0, -16, 8, -8, 8, 8, 0, 0, -8, -8, -16, -16, -8, 8, 0, -16, -8, 8, 0, 16, 0, 0, -8, 24, 16, -48, 8, 8, 8, 8, 16, 16, 8, 8, -16, 16, 0, 16, -8, 8, -16, 32, 8, 24, 8, 8, 0, 0, -8, -8, 16, -16, 0, 16, 8, -8, 0, -16, 8, -8, -8, 8, -16, 0, 8, -8, 0, 16, -32, 0, 8, 8, 0, 0, 8, 8, 0, 0, -8, -8, -16, -16, 8, 8, 8, -8, -16, 0, 8, -8, -16, 0, 0, -16, -8, 8, -16, 0, 8, -8, -8, -8, 0, 0, -8, -8, 0, 0, 12, 4, 4, -4, -4, -12, 20, -20, 4, 12, 12, -12, 4, -20, 12, -12, -4, 4, -12, -4, 12, -12, 4, 12, -28, -4, 12, 4, 4, -4, 12, 4, 4, -4, -4, -12, -12, -20, 12, 4, 12, -12, -12, -4, 12, -12, -12, -4, 4, -20, -4, 4, -12, -4, 12, -12, -4, -12, 4, -4, -4, -12, 4, -4, -4, 4, 4, 12, -4, 4, 4, 12, -12, 12, -20, 4, 4, -4, 60, -12, -4, -12, 4, -4, -4, -12, 4, -4, 4, 12, -4, 4, -12, -4, -20, -12, -12, -20, -4, 84, -12, -20, -4, -12, -4, -28, -12, -4, 12, 52, 4, -20, 4, -20, 12, -12, 4, -20, -20, -12, -4, -12, -12, -20, -20, 4, 4]);
        assert!(boolean_function.is_err());
        assert_eq!(boolean_function.unwrap_err(), InvalidWalshValuesCount(255));

        let boolean_function = super::boolean_function_from_reverse_walsh_hadamard_transform(&[128]);
        assert!(boolean_function.is_err());
        assert_eq!(boolean_function.unwrap_err(), InvalidWalshValuesCount(1));
    }

    #[test]
    fn test_correlation_immunity() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("dd0e").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("55C3AAC3").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 1);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1f").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 0);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 4);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000000000000000000000000000").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 7);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 2);
    }

    #[test]
    fn test_resiliency_order() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("dd0e").unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("55C3AAC3").unwrap();
        assert_eq!(boolean_function.resiliency_order(), Some(1));

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1f").unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00000000000000000000000000000000").unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.resiliency_order(), Some(2));
    }

    #[test]
    fn test_biguint_truth_table() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.biguint_truth_table().to_str_radix(16), "7969817cc5893ba6ac326e47619f5ad0");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.biguint_truth_table().to_str_radix(16), "1e");
    }

    #[test]
    fn test_try_u64_truth_table() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(boolean_function.try_u64_truth_table(), None);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.try_u64_truth_table(), Some(30));
    }

    #[test]
    fn test_xor() {
        let mut boolean_function = super::boolean_function_from_hex_string_truth_table("80921c010276c440400810a80e200425").unwrap();
        let boolean_function2 = super::boolean_function_from_hex_string_truth_table("22442244118811882244224411881188").unwrap();
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "a2d63e4513fed5c8624c32ec1fa815ad");
        assert_eq!(boolean_function.get_num_variables(), 7);
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Big);
        assert_eq!(boolean_function3.get_boolean_function_type(), BooleanFunctionType::Big);

        let mut boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        let boolean_function2 = super::boolean_function_from_hex_string_truth_table("ab").unwrap();
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.get_num_variables(), 3);
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Small);
        assert_eq!(boolean_function3.get_boolean_function_type(), BooleanFunctionType::Small);

        let mut boolean_function = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        let boolean_function2: BooleanFunction = Box::new(super::BigBooleanFunction::from_truth_table(BigUint::from_str_radix("ab", 16).unwrap(), 3));
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.get_num_variables(), 3);
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Small);
        assert_eq!(boolean_function2.get_boolean_function_type(), BooleanFunctionType::Big);
        assert_eq!(boolean_function3.get_boolean_function_type(), BooleanFunctionType::Small);

        let mut boolean_function: BooleanFunction = Box::new(super::BigBooleanFunction::from_truth_table(BigUint::from_str_radix("ab", 16).unwrap(), 3));
        let boolean_function2 = super::boolean_function_from_hex_string_truth_table("1e").unwrap();
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.get_num_variables(), 3);
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Big);
        assert_eq!(boolean_function2.get_boolean_function_type(), BooleanFunctionType::Small);
        assert_eq!(boolean_function3.get_boolean_function_type(), BooleanFunctionType::Big);
    }

    #[test]
    fn test_walsh_fourier_values() {
        let boolean_function = super::boolean_function_from_hex_string_truth_table("ff").unwrap();
        assert_eq!(boolean_function.walsh_fourier_values(), [8, 0, 0, 0, 0, 0, 0, 0]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("00").unwrap();
        assert_eq!(boolean_function.walsh_fourier_values(), [0, 0, 0, 0, 0, 0, 0, 0]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("0f").unwrap();
        assert_eq!(boolean_function.walsh_fourier_values(), [4, 0, 0, 0, 4, 0, 0, 0]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("55").unwrap();
        assert_eq!(boolean_function.walsh_fourier_values(), [4, 4, 0, 0, 0, 0, 0, 0]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("aa").unwrap();
        assert_eq!(boolean_function.walsh_fourier_values(), [4, -4, 0, 0, 0, 0, 0, 0]);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("8001").unwrap();
        assert_eq!(boolean_function.walsh_fourier_values(), [2, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 0, 2]);
    }

    #[test]
    fn test_boolean_function_from_reverse_walsh_fourier_transform() {
        let boolean_function = super::boolean_function_from_reverse_walsh_fourier_transform(&[8, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "ff");

        let boolean_function = super::boolean_function_from_reverse_walsh_fourier_transform(&[0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00");

        let boolean_function = super::boolean_function_from_reverse_walsh_fourier_transform(&[2, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 0, 2]).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "8001");
        assert_eq!(boolean_function.get_boolean_function_type(), BooleanFunctionType::Small);

        let boolean_function = super::boolean_function_from_hex_string_truth_table("80921c010276c44224422441188118822442244118811880400810a80e200425").unwrap();
        let walsh_fourier_values = boolean_function.walsh_fourier_values();
        let boolean_function2 = super::boolean_function_from_reverse_walsh_fourier_transform(&walsh_fourier_values).unwrap();
        assert_eq!(boolean_function2.printable_hex_truth_table(), "80921c010276c44224422441188118822442244118811880400810a80e200425");

        let boolean_function = super::boolean_function_from_hex_string_truth_table("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff").unwrap();
        let walsh_fourier_values = boolean_function.walsh_fourier_values();
        let boolean_function2 = super::boolean_function_from_reverse_walsh_fourier_transform(&walsh_fourier_values).unwrap();
        assert_eq!(boolean_function2.printable_hex_truth_table(), "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
        assert_eq!(boolean_function2.get_boolean_function_type(), BooleanFunctionType::Big);
        assert_eq!(boolean_function2.get_num_variables(), 9);
    }
}
