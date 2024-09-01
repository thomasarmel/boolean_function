use std::any::Any;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use crate::{BooleanFunctionImpl, BooleanFunctionError, BooleanFunction, BooleanFunctionType};

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BigBooleanFunction {
    variables_count: usize,
    truth_table: BigUint,
}

impl BigBooleanFunction {
    pub fn from_truth_table(truth_table: BigUint, variables_count: usize) -> Self { // TODO check too big for variables_count
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
                return Err(BooleanFunctionError::TooBigDerivativeDirection(max_input_value));
            }
        }
        let mut derivative_truth_table: BigUint = BigUint::zero();
        for x in 0..=self.get_max_input_value() {
            if self.compute_cellular_automata_rule(x) ^ self.compute_cellular_automata_rule(x ^ direction) {
                derivative_truth_table |= BigUint::one() << x;
            }
        }
        Ok(BigBooleanFunction {
            variables_count: self.variables_count,
            truth_table: derivative_truth_table,
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

    fn printable_hex_truth_table(&self) -> String {
        format!("{:01$x}", self.truth_table, 1 << (self.variables_count - 2))
    }

    fn __clone(&self) -> BooleanFunction {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::{Num, One, Zero};
    use crate::{BigBooleanFunction, BooleanFunctionImpl};

    #[test]
    fn test_get_num_variables() {
        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(), 7);
        assert_eq!(boolean_function.get_num_variables(), 7);
    }

    #[test]
    fn test_is_balanced() {
        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(), 7);
        assert_eq!(boolean_function.is_balanced(), true);

        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD1", 16).unwrap(), 7);
        assert_eq!(boolean_function.is_balanced(), false);
    }

    #[test]
    fn test_compute_cellular_automata_rule() {
        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(), 7);
        assert_eq!(boolean_function.compute_cellular_automata_rule(13), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(62), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(64), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(80), true);
        assert_eq!(boolean_function.compute_cellular_automata_rule(100), true);
    }

    #[test]
    fn test_derivative() {
        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(), 7);
        let derivative = boolean_function.derivative(1).unwrap();
        assert_eq!(derivative.get_num_variables(), 7);
        assert_eq!(derivative.printable_hex_truth_table(), "cfffc3c00fcf0cfff003f3ccf3f0ff30");

        let derivative = boolean_function.derivative(2).unwrap();
        assert_eq!(derivative.printable_hex_truth_table(), "afffa5aff0aff50f0ffaf55af5f000a0");

        let derivative = boolean_function.derivative(3).unwrap();
        assert_eq!(derivative.printable_hex_truth_table(), "9000999fff90f6f0fff609690900ff60");

        let derivative = boolean_function.derivative(128);
        assert!(derivative.is_err());
    }

    #[test]
    fn test_printable_hex_truth_table() {
        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(), 7);
        assert_eq!(boolean_function.printable_hex_truth_table(), "7969817cc5893ba6ac326e47619f5ad0");

        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::one(), 7);
        assert_eq!(boolean_function.printable_hex_truth_table(), "00000000000000000000000000000001");

        let boolean_function = BigBooleanFunction::from_truth_table(BigUint::zero(), 7);
        assert_eq!(boolean_function.printable_hex_truth_table(), "00000000000000000000000000000000");
    }
}