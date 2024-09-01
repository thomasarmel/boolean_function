use crate::anf_polynom::AnfPolynomial;
use crate::{BooleanFunction, BooleanFunctionError, BooleanFunctionImpl, BooleanFunctionType};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::any::Any;

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

    fn algebraic_normal_form(&self) -> AnfPolynomial {
        let u0 = BigUint::zero();
        let u1 = BigUint::one();

        let mut blocksize = 1usize;
        let mut anf_form = self.truth_table.clone();
        for _ in 0..self.variables_count {
            let mut source = 0usize;
            while source < (1 << self.variables_count) {
                let target = source + blocksize.clone();
                for i in 0..blocksize {
                    let f_target_i: bool = ((anf_form.clone() >> (target + i)) & u1.clone()) != u0;
                    let f_source_i: bool = ((anf_form.clone() >> (source + i)) & u1.clone()) != u0;
                    let f_target_i_xor_f_source_i = f_target_i ^ f_source_i;
                    if f_target_i_xor_f_source_i {
                        anf_form = anf_form | (u1.clone() << (target + i));
                    } else {
                        // set (target + i) bit of final_f to 0 with not operation
                        anf_form.set_bit((target + i) as u64, false);
                    }
                }
                source = source + (blocksize << 1);
            }
            blocksize = blocksize << 1;
        }
        AnfPolynomial::from_truth_table_big(anf_form, self.variables_count)
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
}
