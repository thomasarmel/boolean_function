//! Algebraic Normal Form (ANF) representation of Boolean functions.

#[cfg(not(feature = "unsafe_disable_safety_checks"))]
use crate::boolean_function_error::POLYNOMIAL_ANF_TOO_BIG_VAR_COUNT_PANIC_MSG;
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::fmt::Display;

#[derive(Debug, Clone, Eq, PartialEq)]
enum PolynomialFormat {
    Small(u64),
    Big(BigUint),
}

/// Polynomial representation of Boolean function in Algebraic Normal Form (ANF).
///
/// The ANF representation is a XOR sum of monomials, where each monomial is a AND product of variables.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct AnfPolynomial {
    polynomial: PolynomialFormat,
    num_variables: usize,
}

impl AnfPolynomial {
    pub(crate) fn from_anf_big(polynomial: &BigUint, num_variables: usize) -> Self {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if polynomial.bits() > (1 << num_variables) {
            panic!("{}", POLYNOMIAL_ANF_TOO_BIG_VAR_COUNT_PANIC_MSG);
        }
        AnfPolynomial {
            polynomial: PolynomialFormat::Big(polynomial.clone()),
            num_variables,
        }
    }

    pub(crate) fn from_anf_small(polynomial: u64, num_variables: usize) -> Self {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if polynomial > (u64::MAX >> (64 - (1 << num_variables))) {
            panic!("{}", POLYNOMIAL_ANF_TOO_BIG_VAR_COUNT_PANIC_MSG);
        }
        AnfPolynomial {
            polynomial: PolynomialFormat::Small(polynomial),
            num_variables,
        }
    }

    /// Get the polynomial in ANF representation in unsigned 64-bit integer format. The unsigned representation of ANF is explained in README.md.
    ///
    /// Returns `None` if the polynomial is too big to fit in a 64-bit integer (ie more than 6 variables).
    pub fn get_polynomial_small(&self) -> Option<u64> {
        match self.polynomial {
            PolynomialFormat::Small(p) => Some(p),
            PolynomialFormat::Big(_) => None,
        }
    }

    /// Get the polynomial in ANF representation in big unsigned integer format. The unsigned representation of ANF is explained in README.md.
    pub fn get_polynomial_big(&self) -> BigUint {
        match &self.polynomial {
            PolynomialFormat::Small(p) => BigUint::from(*p),
            PolynomialFormat::Big(p) => p.clone(),
        }
    }

    /// Get the degree of the polynomial.
    ///
    /// The degree of the polynomial is the maximum number of variables in a monomial.
    pub fn get_degree(&self) -> usize {
        let max_input_value: u32 = (1 << self.num_variables) - 1;
        match &self.polynomial {
            PolynomialFormat::Small(p) => (0..=max_input_value)
                .into_iter()
                .map(|bit_position| {
                    if p & (1 << bit_position) != 0 {
                        bit_position.count_ones() as usize
                    } else {
                        0
                    }
                })
                .max()
                .unwrap_or(0),
            PolynomialFormat::Big(p) => (0..=max_input_value)
                .into_iter()
                .map(|bit_position| {
                    if p & (BigUint::one() << bit_position) != BigUint::zero() {
                        bit_position.count_ones() as usize
                    } else {
                        0
                    }
                })
                .max()
                .unwrap_or(0),
        }
    }

    /// Returns the string representation of the polynomial.
    ///
    /// The monomials are ordered by the number of variables in the monomial in descending order, and then by the lexicographic order of the variables.
    ///
    /// Example: "x0\*x1 + x0 + x1 + x2", the '+' operator is the XOR operator, and the '\*' operator is the AND operator.
    pub fn to_string(&self) -> String {
        let mut monomials_str_list: Vec<(String, usize)> = match &self.polynomial {
            PolynomialFormat::Small(p) => (0..(1 << self.num_variables))
                .into_iter()
                .filter(|bit_position| p & (1u64 << bit_position) != 0)
                .map(|bit_position| self.bit_monomial_to_string(bit_position))
                .collect(),
            PolynomialFormat::Big(p) => (0..(1 << self.num_variables))
                .into_iter()
                .filter(|bit_position| p.bit(*bit_position))
                .map(|bit_position| self.bit_monomial_to_string(bit_position))
                .collect(),
        };
        if monomials_str_list.is_empty() {
            return String::from("0");
        }
        monomials_str_list.sort_by(|a, b| b.1.cmp(&a.1));
        monomials_str_list
            .iter()
            .map(|(monomial, _)| monomial)
            .join(" + ")
    }

    /// Monomial and number of variables in the monomial
    fn bit_monomial_to_string(&self, bit_position: u64) -> (String, usize) {
        if bit_position == 0 {
            return (String::from("1"), 0);
        }
        let monomial_coefs_list: Vec<String> = (0..self.num_variables)
            .filter(|i| bit_position & (1u64 << i) != 0)
            .map(|i| String::from(&format!("x{}", i)))
            .collect();
        (monomial_coefs_list.join("*"), monomial_coefs_list.len())
    }
}

/// Display implementation for `AnfPolynomial`.
///
/// Internally uses the [AnfPolynomial::to_string] method.
impl Display for AnfPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

#[cfg(test)]
mod tests {
    use crate::anf_polynom::AnfPolynomial;
    use num_bigint::BigUint;
    use num_traits::{Num, One, Zero};

    #[test]
    fn test_get_polynomial_small() {
        let anf_polynomial = AnfPolynomial::from_anf_small(30, 3);
        assert_eq!(anf_polynomial.get_polynomial_small(), Some(30));

        let anf_polynomial =
            AnfPolynomial::from_anf_big(&BigUint::from_str_radix("30", 16).unwrap(), 3);
        assert_eq!(anf_polynomial.get_polynomial_small(), None);
    }

    #[test]
    fn test_get_polynomial_big() {
        let anf_polynomial = AnfPolynomial::from_anf_small(30, 3);
        assert_eq!(anf_polynomial.get_polynomial_big(), BigUint::from(30u32));

        let anf_polynomial = AnfPolynomial::from_anf_big(
            &BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap(),
            7,
        );
        assert_eq!(
            anf_polynomial.get_polynomial_big(),
            BigUint::from_str_radix("7969817CC5893BA6AC326E47619F5AD0", 16).unwrap()
        );
    }

    #[test]
    fn test_get_degree() {
        let anf_polynomial = AnfPolynomial::from_anf_small(0, 3);
        assert_eq!(anf_polynomial.get_degree(), 0);

        let anf_polynomial = AnfPolynomial::from_anf_small(1, 3);
        assert_eq!(anf_polynomial.get_degree(), 0);

        let anf_polynomial = AnfPolynomial::from_anf_small(0xff, 3);
        assert_eq!(anf_polynomial.get_degree(), 3);

        let anf_polynomial = AnfPolynomial::from_anf_big(&BigUint::zero(), 3);
        assert_eq!(anf_polynomial.get_degree(), 0);

        let anf_polynomial = AnfPolynomial::from_anf_big(&BigUint::one(), 3);
        assert_eq!(anf_polynomial.get_degree(), 0);

        let anf_polynomial = AnfPolynomial::from_anf_big(
            &BigUint::from_str_radix("00000000000000000000000000000010", 16).unwrap(),
            7,
        );
        assert_eq!(anf_polynomial.get_degree(), 1);

        let anf_polynomial = AnfPolynomial::from_anf_big(
            &BigUint::from_str_radix("00000000000000000000000000001110", 16).unwrap(),
            7,
        );
        assert_eq!(anf_polynomial.get_degree(), 2);

        let anf_polynomial = AnfPolynomial::from_anf_big(
            &BigUint::from_str_radix("f0000000000000000000000000001110", 16).unwrap(),
            7,
        );
        assert_eq!(anf_polynomial.get_degree(), 7);
    }

    #[test]
    fn test_to_string() {
        let anf_polynomial = AnfPolynomial::from_anf_small(30, 3);
        assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + x1 + x2");

        let anf_polynomial = AnfPolynomial::from_anf_small(31, 3);
        assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + x1 + x2 + 1");

        let anf_polynomial = AnfPolynomial::from_anf_small(0, 3);
        assert_eq!(anf_polynomial.to_string(), "0");

        let anf_polynomial = AnfPolynomial::from_anf_small(1, 3);
        assert_eq!(anf_polynomial.to_string(), "1");

        let anf_polynomial = AnfPolynomial::from_anf_big(&BigUint::from(30u32), 3);
        assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + x1 + x2");

        let anf_polynomial = AnfPolynomial::from_anf_big(&BigUint::from(31u32), 3);
        assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + x1 + x2 + 1");

        let anf_polynomial = AnfPolynomial::from_anf_big(&BigUint::zero(), 3);
        assert_eq!(anf_polynomial.to_string(), "0");

        let anf_polynomial = AnfPolynomial::from_anf_big(&BigUint::one(), 3);
        assert_eq!(anf_polynomial.to_string(), "1");
    }
}
