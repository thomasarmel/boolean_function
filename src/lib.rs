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
use num_traits::Num;
pub use small_boolean_function::SmallBooleanFunction;
use std::any::Any;
use std::collections::HashMap;
use std::fmt::Debug;

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
    fn absolute_walsh_spectrum(&self) -> HashMap<u32, usize> {
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
    fn annihilator(&self, max_degree: usize) -> Option<(BooleanFunction, usize, usize)>; // TODO max degree in Option

    fn algebraic_immunity(&self) -> usize {
        match self.annihilator(self.get_num_variables()) {
            None => 0,
            Some(annihilator) => annihilator.1
        }
    }

    fn printable_hex_truth_table(&self) -> String;

    /// Use Clone instead of this method
    fn __clone(&self) -> BooleanFunction;

    fn as_any(&self) -> &dyn Any;

    // TODO almost bent, mul (and tt), sum, impl not, iterate on values
}

pub type BooleanFunction = Box<dyn BooleanFunctionImpl>;

impl Clone for BooleanFunction {
    fn clone(&self) -> Self {
        self.__clone()
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

fn boolean_function_from_hex_string_truth_table(
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

#[cfg(test)]
mod tests {
    use crate::BooleanFunctionType;
    use std::collections::HashMap;

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
    fn test_absolute_walsh_spectrum() {
        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.absolute_walsh_spectrum(),
            HashMap::from([(0, 64), (16, 64)])
        );

        let boolean_function =
            super::boolean_function_from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(
            boolean_function.absolute_walsh_spectrum(),
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
}
