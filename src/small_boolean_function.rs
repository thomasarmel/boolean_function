use crate::anf_polynom::AnfPolynomial;
use crate::BooleanFunctionError::TooBigVariableCount;
use crate::{BooleanFunction, BooleanFunctionError, BooleanFunctionImpl, BooleanFunctionType};
use fast_boolean_anf_transform::fast_bool_anf_transform_unsigned;
use std::any::Any;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct SmallBooleanFunction {
    variables_count: usize,
    truth_table: u64,
}

impl SmallBooleanFunction {
    /// Length must be < 6
    pub fn from_truth_table(
        truth_table: u64,
        variables_count: usize,
    ) -> Result<Self, BooleanFunctionError> {
        if variables_count > 6 {
            // TODO skip test if needed, maybe create another function
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
}

impl BooleanFunctionImpl for SmallBooleanFunction {
    #[inline]
    fn get_num_variables(&self) -> usize {
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

    fn algebraic_normal_form(&self) -> AnfPolynomial {
        let anf_form = fast_bool_anf_transform_unsigned(self.truth_table, self.variables_count);
        AnfPolynomial::from_truth_table_small(anf_form, self.variables_count)
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
    use crate::{BooleanFunctionImpl, SmallBooleanFunction};

    #[test]
    fn test_from_truth_table() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert_eq!(boolean_function.get_truth_table_u64(), 0xaa55aa55);

        let boolean_function = SmallBooleanFunction::from_truth_table(0xffffffffffffffff, 7);
        assert!(boolean_function.is_err());
    }

    #[test]
    fn test_get_num_variables() {
        let boolean_function = SmallBooleanFunction::from_truth_table(0xaa55aa55, 5).unwrap();
        assert_eq!(boolean_function.get_num_variables(), 5);
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
        assert_eq!(derivative.get_num_variables(), 5);
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
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1 + x0 + x1 + x2");

        let boolean_function = SmallBooleanFunction::from_truth_table(0, 3).unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "0");

        let boolean_function = SmallBooleanFunction::from_truth_table(110, 3).unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2 + x0*x1 + x0 + x1");

        let boolean_function = SmallBooleanFunction::from_truth_table(123, 3).unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1 + x1*x2 + x1 + 1");
    }
}
