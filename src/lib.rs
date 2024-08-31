mod big_boolean_function;
mod small_boolean_function;
mod utils;
mod boolean_function_error;

use std::collections::HashMap;
use std::fmt::Debug;
use num_bigint::BigUint;
use num_traits::Num;
pub use big_boolean_function::BigBooleanFunction;
pub use small_boolean_function::SmallBooleanFunction;
pub use boolean_function_error::BooleanFunctionError;
use crate::BooleanFunctionError::{StringHexParseError, UnexpectedError, WrongStringHexTruthTableLength};

pub trait BooleanFunctionImpl: Debug {
    fn get_num_variables(&self) -> usize;
    #[inline]
    fn get_max_input_value(&self) -> u32 {
        (1 << self.get_num_variables()) - 1
    }
    fn is_balanced(&self) -> bool;

    /// debug_assertions
    fn compute_cellular_automata_rule(&self, input_bits: u32) -> bool;
    fn walsh_transform(&self, w: u32) -> i32 {
        (0..=self.get_max_input_value()).map(|x| {
            if (self.compute_cellular_automata_rule(x) as u32 + utils::fast_binary_dot_product(w, x)) & 1 == 0 { // % modulo 2
                1
            } else {
                -1
            }
        }).sum()
    }
    fn absolute_walsh_spectrum(&self) -> HashMap<u32, usize> {
        let mut absolute_walsh_value_count_map: HashMap<u32, usize> = HashMap::new();
        (0..=self.get_max_input_value())
            .for_each(|w| {
                let absolute_walsh_value = self.walsh_transform(w).unsigned_abs();
                if !absolute_walsh_value_count_map.contains_key(&absolute_walsh_value) {
                    absolute_walsh_value_count_map.insert(absolute_walsh_value, 1);
                } else {
                    let count = absolute_walsh_value_count_map.get_mut(&absolute_walsh_value).unwrap();
                    *count += 1;
                }
            });
        absolute_walsh_value_count_map
    }

    fn auto_correlation_transform(&self, w: u32) -> i32 {
        (0..=self.get_max_input_value()).map(|x| {
            if self.compute_cellular_automata_rule(x) ^ self.compute_cellular_automata_rule(x ^ w) {
                -1
            } else {
                1
            }
        }).sum()
    }

    fn absolute_autocorrelation(&self) -> HashMap<u32, usize> {
        let mut absolute_autocorrelation_value_count_map: HashMap<u32, usize> = HashMap::new();
        (0..=self.get_max_input_value())
            .for_each(|w| {
                let absolute_autocorrelation_value = self.auto_correlation_transform(w).unsigned_abs();
                if !absolute_autocorrelation_value_count_map.contains_key(&absolute_autocorrelation_value) {
                    absolute_autocorrelation_value_count_map.insert(absolute_autocorrelation_value, 1);
                } else {
                    let count = absolute_autocorrelation_value_count_map.get_mut(&absolute_autocorrelation_value).unwrap();
                    *count += 1;
                }
            });
        absolute_autocorrelation_value_count_map
    }

    fn derivative(&self, direction: u32) -> Result<BooleanFunction, BooleanFunctionError>;

    fn printable_hex_truth_table(&self) -> String;

    /// Use Clone instead of this method
    fn __clone(&self) -> BooleanFunction;
}

pub type BooleanFunction = Box<dyn BooleanFunctionImpl>;

impl Clone for BooleanFunction {
    fn clone(&self) -> Self {
        self.__clone()
    }
}

fn boolean_function_from_hex_string_truth_table(hex_truth_table: &str) -> Result<BooleanFunction, BooleanFunctionError> {
    if hex_truth_table.len().count_ones() != 1 { // TODO test + parsing
        return Err(WrongStringHexTruthTableLength);
    }
    let num_variables = (hex_truth_table.len() << 2).trailing_zeros() as usize;
    if num_variables <= 6 {
        Ok(Box::new(SmallBooleanFunction::from_truth_table(u64::from_str_radix(hex_truth_table, 16)
                                                               .map_err(|_| StringHexParseError)?, num_variables)
            .map_err(|_| UnexpectedError)?))
    } else {
        Ok(Box::new(BigBooleanFunction::from_truth_table(BigUint::from_str_radix(hex_truth_table, 16).map_err(|_| StringHexParseError)?, num_variables)))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    #[test]
    fn test_absolute_autocorrelation_spectrum() {
        let big_boolean_function = super::boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        assert_eq!(big_boolean_function.absolute_autocorrelation(), HashMap::from([(0, 33), (8, 58), (16, 28), (24, 6), (32, 2), (128, 1)]));
    }
}