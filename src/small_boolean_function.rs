use crate::{BooleanFunctionImpl, BooleanFunctionError, BooleanFunction};
use crate::BooleanFunctionError::TooBigVariableCount;

#[derive(Debug, Clone, Copy)]
pub struct SmallBooleanFunction {
    variables_count: usize,
    truth_table: u64,
}

impl SmallBooleanFunction {
    /// Length must be < 6
    pub fn from_truth_table(truth_table: u64, variables_count: usize) -> Result<Self, BooleanFunctionError> {
        if variables_count > 6 { // TODO skip test if needed, maybe create another function
            return Err(TooBigVariableCount(6));
        }
        Ok(SmallBooleanFunction {
            variables_count,
            truth_table,
        })
    }

    pub fn get_truth_table_u64(&self) -> u64 {
        self.truth_table
    }

    pub fn derivative_inner(&self, direction: u32) -> Result<Self, BooleanFunctionError> {
        #[cfg(debug_assertions)]
        {
            let max_input_value = self.get_max_input_value();
            if direction > max_input_value {
                return Err(BooleanFunctionError::TooBigDerivativeDirection(max_input_value));
            }
        }
        let mut derivative_truth_table: u64 = 0;
        for x in 0..=self.get_max_input_value() {
            if self.compute_cellular_automata_rule(x) ^ self.compute_cellular_automata_rule(x ^ direction) {
                derivative_truth_table |= 1 << x;
            }
        }
        Ok(SmallBooleanFunction {
            variables_count: self.variables_count,
            truth_table: derivative_truth_table,
        })
    }
}

impl BooleanFunctionImpl for SmallBooleanFunction {
    #[inline]
    fn get_num_variables(&self) -> usize {
        self.variables_count
    }

    fn is_balanced(&self) -> bool {
        let expected_set_number: u32 = 1 << (self.variables_count - 1);
        self.truth_table.count_ones() == expected_set_number
    }

    #[inline]
    fn compute_cellular_automata_rule(&self, input_bits: u32) -> bool {
        #[cfg(debug_assertions)]
        {
            let max_input_value = self.get_max_input_value();
            if input_bits > max_input_value {
                panic!("Input bits must be less or equal than {}", max_input_value);
            }
        }
        (self.truth_table & (1 << input_bits)) != 0
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

#[cfg(test)]
mod tests {
    use crate::{BooleanFunctionImpl, SmallBooleanFunction};

    #[test]
    fn test_derivative() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        let derivative = boolean_function.derivative_inner(1);
        assert_eq!(derivative.unwrap().get_truth_table_u64(), 0xffffffff);

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaaddbb55, 5).unwrap();
        let derivative = boolean_function.derivative_inner(1);
        assert_eq!(derivative.unwrap().get_truth_table_u64(), 0xff33ccff);

        let derivative = boolean_function.derivative_inner(2);
        assert_eq!(derivative.unwrap().get_truth_table_u64(), 0x00aa5500);

        let derivative = boolean_function.derivative_inner(3);
        assert_eq!(derivative.unwrap().get_truth_table_u64(), 0xff6666ff);

        let derivative = boolean_function.derivative_inner(32);
        assert!(derivative.is_err());
    }

    #[test]
    fn test_printable_hex_truth_table() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "aa55aa55");

        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55, 5).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "0000aa55");
    }
}