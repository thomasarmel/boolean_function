use num_bigint::BigUint;
use num_traits::{One, Zero};
use crate::{BooleanFunctionImpl, BooleanFunctionError, BooleanFunction};

#[derive(Debug, Clone)]
pub struct BigBooleanFunction {
    variables_count: usize,
    truth_table: BigUint,
}

impl BigBooleanFunction {
    pub fn from_truth_table(truth_table: BigUint, variables_count: usize) -> Self {
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
}