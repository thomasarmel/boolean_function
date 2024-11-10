//! Algebraic Normal Form (ANF) representation of Boolean functions.

#[cfg(not(feature = "unsafe_disable_safety_checks"))]
use crate::boolean_function_error::{POLYNOMIAL_ANF_TOO_BIG_VAR_COUNT_PANIC_MSG, XOR_DIFFERENT_VAR_COUNT_PANIC_MSG, AND_DIFFERENT_VAR_COUNT_PANIC_MSG};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};
use std::fmt::Display;
use std::ops::{Add, AddAssign, BitAnd, BitAndAssign, BitXor, BitXorAssign, Mul, MulAssign};
use fast_boolean_anf_transform::fast_bool_anf_transform_unsigned;
use crate::{BigBooleanFunction, BooleanFunction, BooleanFunctionError, SmallBooleanFunction};
use crate::utils::fast_anf_transform_biguint;

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

    /// Computes the ANF polynomial from string representation
    ///
    /// Representation must be in the form "`x0*x2*x3 + x2*x3 + x1 + 1`".
    ///
    /// X's index starts at 0, meaning the maximum index is variable count - 1.
    ///
    /// # Parameters:
    /// - `anf_polynomial`: The string representation of the ANF form
    /// - `num_variables`: Variable count of the polynomial
    ///
    /// # Returns:
    /// The ANF polynomial, or an error if the input string doesn't respect the format and `unsafe_disable_safety_checks` feature is not activated.
    pub fn from_str(anf_polynomial: &str, num_variables: usize) -> Result<Self, BooleanFunctionError> {
        let anf_polynomial_filtered = anf_polynomial.chars().filter(|c| !c.is_whitespace()).collect::<String>();

        let mut anf_polynomial_obj = Self {
            polynomial: if num_variables <= 6 {
                PolynomialFormat::Small(0u64)
            } else {
                PolynomialFormat::Big(BigUint::zero())
            },
            num_variables,
        };

        if anf_polynomial.is_empty() {
            return Ok(anf_polynomial_obj);
        }

        for monomial_string in anf_polynomial_filtered.split('+') {
            if monomial_string == "0" {
                continue;
            }
            if monomial_string == "1" {
                anf_polynomial_obj.flip_bit_pos(0)?;
                continue;
            }
            let mut bit_position_to_flip = 0u64;
            for factor in monomial_string.split('*') {
                #[cfg(not(feature = "unsafe_disable_safety_checks"))]
                if factor.chars().nth(0).ok_or(BooleanFunctionError::ErrorParsingAnfString)? != 'x' {
                    return Err(BooleanFunctionError::ErrorParsingAnfString);
                }
                let factor_index = factor[1..].parse::<u64>().map_err(|_| BooleanFunctionError::ErrorParsingAnfString)?;
                #[cfg(not(feature = "unsafe_disable_safety_checks"))]
                if factor_index as usize >= num_variables {
                    return Err(BooleanFunctionError::AnfFormNVariableTooBigFactor(num_variables, factor_index as usize));
                }
                bit_position_to_flip |= 1 << factor_index;
            }
            anf_polynomial_obj.flip_bit_pos(bit_position_to_flip)?;
        }

        Ok(anf_polynomial_obj)
    }

    fn flip_bit_pos(&mut self, bit_position: u64) -> Result<(), BooleanFunctionError> {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if bit_position >= 1 << self.num_variables {
            return Err(BooleanFunctionError::UnexpectedError);
        }
        match &mut self.polynomial {
            PolynomialFormat::Small(anf) => {
                *anf ^= 1u64 << bit_position
            },
            PolynomialFormat::Big(anf) => {
                anf.set_bit(bit_position, !anf.bit(bit_position))
            }
        }
        Ok(())
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

    /// Returns Boolean function internal storage type:
    ///
    /// - [Small](crate::BooleanFunctionType::Small) for `u64` internal storage (less or equal than 6 variables Boolean function)
    /// - [Big](crate::BooleanFunctionType::Big) for `BigUInt` internal storage (more than 6 variables Boolean function)
    pub fn get_boolean_function_type(&self) -> crate::BooleanFunctionType {
        match self.polynomial {
            PolynomialFormat::Small(_) => crate::BooleanFunctionType::Small,
            PolynomialFormat::Big(_) => crate::BooleanFunctionType::Big
        }
    }

    /// Convert ANF polynomial to the corresponding Boolean Function, using fast ANF transform algorithm
    /// 
    /// # Returns
    /// A Boolean function corresponding to the polynomial
    pub fn to_boolean_function(&self) -> BooleanFunction {
        match &self.polynomial {
            PolynomialFormat::Small(polynomial) => {
                BooleanFunction::Small(
                    SmallBooleanFunction::from_truth_table(
                        fast_bool_anf_transform_unsigned(*polynomial, self.num_variables),
                        self.num_variables
                    ).unwrap()
                )
            }
            PolynomialFormat::Big(polynomial) => {
                BooleanFunction::Big(
                    BigBooleanFunction::from_truth_table(
                        fast_anf_transform_biguint(polynomial, self.num_variables),
                        self.num_variables
                    )
                )
            }
        }
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

impl Into<BooleanFunction> for AnfPolynomial {
    fn into(self) -> BooleanFunction {
        self.to_boolean_function()
    }
}

impl Into<String> for AnfPolynomial {
    fn into(self) -> String {
        self.to_string()
    }
}

/// In-place XOR operator for Boolean functions ANF polynomial
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitXorAssign for AnfPolynomial {
    fn bitxor_assign(&mut self, rhs: Self) {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if self.num_variables != rhs.num_variables {
            panic!("{}", XOR_DIFFERENT_VAR_COUNT_PANIC_MSG);
        }
        match (&mut self.polynomial, &rhs.polynomial) {
            (PolynomialFormat::Small(self_poly), PolynomialFormat::Small(rhs_poly)) => {
                *self_poly ^= rhs_poly;
            },
            (PolynomialFormat::Big(self_poly), PolynomialFormat::Small(rhs_poly)) => {
                *self_poly ^= BigUint::from_u64(*rhs_poly).unwrap();
            },
            (PolynomialFormat::Small(self_poly), PolynomialFormat::Big(rhs_poly)) => {
                *self_poly ^= rhs_poly.to_u64().unwrap();
            },
            (PolynomialFormat::Big(self_poly), PolynomialFormat::Big(rhs_poly)) => {
                *self_poly ^= rhs_poly;
            },
        }
    }
}

/// XOR operator for Boolean functions ANF polynomial
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitXor for AnfPolynomial {
    type Output = AnfPolynomial;

    fn bitxor(mut self, rhs: Self) -> Self::Output {
        self ^= rhs;
        self
    }
}

/// ADD operator for Boolean functions ANF polynomial.
///
/// It is equivalent to [crate::AnfPolynomial::bitxor] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl Add for AnfPolynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        self ^ rhs
    }
}

/// In-place ADD operator for Boolean functions ANF polynomial.
///
/// It is equivalent to [crate::AnfPolynomial::bitxor_assign] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl AddAssign for AnfPolynomial {
    fn add_assign(&mut self, rhs: Self) {
        *self ^= rhs;
    }
}

/// In-place AND operator for Boolean functions ANF polynomial
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitAndAssign for AnfPolynomial {
    fn bitand_assign(&mut self, rhs: Self) {
        *self = self.clone() & rhs;
    }
}

/// AND operator for Boolean functions ANF polynomial
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitAnd for AnfPolynomial {
    type Output = AnfPolynomial;

    fn bitand(self, rhs: Self) -> Self::Output {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if self.num_variables != rhs.num_variables {
            panic!("{}", AND_DIFFERENT_VAR_COUNT_PANIC_MSG);
        }

        fn small_anf_polynomial_multiply(left_poly: u64, right_poly: u64, num_variables: usize) -> u64 {
            let mut res = 0u64;
            for left_bit_pos in 0u64..(1 << num_variables) {
                if (left_poly >> left_bit_pos) & 1 == 0 {
                    continue;
                }
                for right_bit_pos in 0u64..(1 << num_variables) {
                    if (right_poly >> right_bit_pos) & 1 == 0 {
                        continue;
                    }
                    res ^= 1 << (left_bit_pos | right_bit_pos);
                }
            }
            res
        }

        fn big_anf_polynomial_multiply(left_poly: &BigUint, right_poly: &BigUint, num_variables: usize) -> BigUint {
            let mut res = BigUint::zero();
            for left_bit_pos in 0u64..(1 << num_variables) {
                if !left_poly.bit(left_bit_pos) {
                    continue;
                }
                for right_bit_pos in 0u64..(1 << num_variables) {
                    if !right_poly.bit(right_bit_pos) {
                        continue;
                    }
                    let pos_to_flip = left_bit_pos | right_bit_pos;
                    res.set_bit(pos_to_flip, !res.bit(pos_to_flip));
                }
            }
            res
        }

        let new_polynomial = match (self.polynomial, rhs.polynomial) {
            (PolynomialFormat::Small(self_poly), PolynomialFormat::Small(rhs_poly)) => {
                PolynomialFormat::Small(small_anf_polynomial_multiply(self_poly, rhs_poly, self.num_variables))
            },
            (PolynomialFormat::Big(self_poly), PolynomialFormat::Small(rhs_poly)) => {
                PolynomialFormat::Big(
                    BigUint::from_u64(
                        small_anf_polynomial_multiply(self_poly.to_u64().unwrap(), rhs_poly, self.num_variables)
                    ).unwrap()
                )
            },
            (PolynomialFormat::Small(self_poly), PolynomialFormat::Big(rhs_poly)) => {
                PolynomialFormat::Small(
                    small_anf_polynomial_multiply(self_poly, rhs_poly.to_u64().unwrap(), self.num_variables)
                )
            },
            (PolynomialFormat::Big(self_poly), PolynomialFormat::Big(rhs_poly)) => {
                PolynomialFormat::Big(big_anf_polynomial_multiply(&self_poly, &rhs_poly, self.num_variables))
            }
        };

        AnfPolynomial {
            polynomial: new_polynomial,
            num_variables: self.num_variables
        }
    }
}

/// MUL operator for Boolean functions ANF polynomial.
///
/// It is equivalent to [crate::AnfPolynomial::bitand] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl Mul for AnfPolynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        self & rhs
    }
}

/// In-place MUL operator for Boolean functions ANF polynomial.
///
/// It is equivalent to [crate::AnfPolynomial::bitand_assign] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl MulAssign for AnfPolynomial {
    fn mul_assign(&mut self, rhs: Self) {
        *self &= rhs;
    }
}

#[cfg(test)]
mod tests {
    use crate::anf_polynom::AnfPolynomial;
    use num_bigint::BigUint;
    use num_traits::{Num, One, Zero};
    use crate::{BooleanFunction, BooleanFunctionError, BooleanFunctionImpl};

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

    #[test]
    fn test_from_str() {
        let anf_str = "x0*x1 + x0 + x0 + x0 + 1";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 3).unwrap();
        assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + 1");
        assert_eq!(anf_polynomial.get_boolean_function_type(), crate::BooleanFunctionType::Small);

        let anf_str = "x2 + x0*x1 + x0";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 3).unwrap();
        assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + x2");
        assert_eq!(anf_polynomial.get_boolean_function_type(), crate::BooleanFunctionType::Small);

        let anf_str = "x2 + x0*x1 + x0 + 0 + 0";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 3).unwrap();
        assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + x2");
        assert_eq!(anf_polynomial.get_boolean_function_type(), crate::BooleanFunctionType::Small);

        let anf_str = "x0*x1 + x2 + x0 + 1";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 2);
        assert!(anf_polynomial.is_err());
        assert_eq!(anf_polynomial.unwrap_err(), BooleanFunctionError::AnfFormNVariableTooBigFactor(2, 2));

        let anf_str = "x0*y1 + x2 + x0 + 1";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 3);
        assert!(anf_polynomial.is_err());
        assert_eq!(anf_polynomial.unwrap_err(), BooleanFunctionError::ErrorParsingAnfString);

        let anf_str = "x0*xy1 + x2 + x0 + 1";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 3);
        assert!(anf_polynomial.is_err());
        assert_eq!(anf_polynomial.unwrap_err(), BooleanFunctionError::ErrorParsingAnfString);

        let anf_str = "x0*x1*x2*x3*x4*x5*x6 + x7 + x6 + x6 + 1";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 8).unwrap();
        assert_eq!(anf_polynomial.to_string(), "x0*x1*x2*x3*x4*x5*x6 + x7 + 1");
        assert_eq!(anf_polynomial.get_boolean_function_type(), crate::BooleanFunctionType::Big);

        let anf_str = "x0*x1*x2*x3*x4*x5*x6 + x7";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 8).unwrap();
        assert_eq!(anf_polynomial.to_string(), "x0*x1*x2*x3*x4*x5*x6 + x7");
        assert_eq!(anf_polynomial.get_boolean_function_type(), crate::BooleanFunctionType::Big);

        let anf_str = "x0x1*x2*x3*x4*x5*x6 + x7";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 8);
        assert!(anf_polynomial.is_err());
        assert_eq!(anf_polynomial.unwrap_err(), BooleanFunctionError::ErrorParsingAnfString);

        let anf_str = "";
        let polynomial = AnfPolynomial::from_str(anf_str, 3).unwrap();
        let boolean_function = polynomial.to_boolean_function();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00");

        let anf_str = "0";
        let polynomial = AnfPolynomial::from_str(anf_str, 3).unwrap();
        let boolean_function = polynomial.to_boolean_function();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00");

        let anf_str = "";
        let polynomial = AnfPolynomial::from_str(anf_str, 7).unwrap();
        let boolean_function = polynomial.to_boolean_function();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00000000000000000000000000000000");

        let anf_str = "0";
        let polynomial = AnfPolynomial::from_str(anf_str, 7).unwrap();
        let boolean_function = polynomial.to_boolean_function();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00000000000000000000000000000000");

        let anf_str = "x0*x0";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 3).unwrap();
        assert_eq!(anf_polynomial.to_string(), "x0");
        assert_eq!(anf_polynomial.get_boolean_function_type(), crate::BooleanFunctionType::Small);
    }

    #[test]
    fn test_to_boolean_function() {
        let rule_30_anf_str = "x0*x1 + x0 + x1 + x2";
        let rule_30_polynomial = AnfPolynomial::from_str(rule_30_anf_str, 3).unwrap();
        let rule_30_function = rule_30_polynomial.to_boolean_function();
        assert_eq!(rule_30_function.printable_hex_truth_table(), "1e");

        let rule_30_anf_str = "x0*x1 + x0 + x1 + x2";
        let rule_30_polynomial = AnfPolynomial::from_str(rule_30_anf_str, 3).unwrap();
        let rule_30_function: BooleanFunction = rule_30_polynomial.into();
        assert_eq!(rule_30_function.printable_hex_truth_table(), "1e");

        let anf_str = "x0*x1*x2*x3*x4*x5*x6 + x7";
        let anf_polynomial = AnfPolynomial::from_str(anf_str, 8).unwrap();
        let boolean_function = anf_polynomial.to_boolean_function();
        assert_eq!(boolean_function.printable_hex_truth_table(), "7fffffffffffffffffffffffffffffff80000000000000000000000000000000");
    }

    #[test]
    fn test_xor() {
        let anf_1 = AnfPolynomial::from_str("x0*x1*x4*x7 + x2*x3 + x1 + 0", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1*x2*x4*x7 + x2*x3 + x3 + 1", 8).unwrap();
        let anf_3 = anf_1.clone() ^ anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2*x4*x7 + x0*x1*x4*x7 + x1 + x3 + 1");
        anf_2 ^= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2*x4*x7 + x0*x1*x4*x7 + x1 + x3 + 1");

        let anf_1 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1*x2 + x0*x2 + x1*x0 + x1 + 1", 3).unwrap();
        let anf_3 = anf_1.clone() ^ anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2 + x0*x2 + x0 + x2 + 1");
        anf_2 ^= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2 + x0*x2 + x0 + x2 + 1");
    }

    #[test]
    fn test_add() {
        let anf_1 = AnfPolynomial::from_str("x0*x1*x4*x7 + x2*x3 + x1 + 0", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1*x2*x4*x7 + x2*x3 + x3 + 1", 8).unwrap();
        let anf_3 = anf_1.clone() + anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2*x4*x7 + x0*x1*x4*x7 + x1 + x3 + 1");
        anf_2 += anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2*x4*x7 + x0*x1*x4*x7 + x1 + x3 + 1");

        let anf_1 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1*x2 + x0*x2 + x1*x0 + x1 + 1", 3).unwrap();
        let anf_3 = anf_1.clone() + anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2 + x0*x2 + x0 + x2 + 1");
        anf_2 += anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2 + x0*x2 + x0 + x2 + 1");
    }

    #[test]
    fn test_and() {
        let anf_1 = AnfPolynomial::from_str("x0*x1", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x2", 3).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2");

        let anf_1 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1 + x0*x1*x2 + x1 + 1", 3).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2 + x1*x2 + x0 + x2");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2 + x1*x2 + x0 + x2");

        let anf_1 = AnfPolynomial::from_str("x0*x1 + x0*x1*x2 + x1 + 1", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2 + x1*x2 + x0 + x2");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2 + x1*x2 + x0 + x2");

        let anf_1 = AnfPolynomial::from_str("x3*x2*x1 + x0*x1 + x0*x1*x2 + x1 + 1", 4).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x2*x3 + x0*x1 + x0 + x1 + x2", 4).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2 + x0*x2*x3 + x1*x2 + x0 + x2");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2 + x0*x2*x3 + x1*x2 + x0 + x2");

        let anf_1 = AnfPolynomial::from_str("x3*x2*x1 + x0*x1 + x0*x1*x2 + x1 + 1", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x2*x3 + x0*x1 + x0 + x1 + x2", 8).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2 + x0*x2*x3 + x1*x2 + x0 + x2");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2 + x0*x2*x3 + x1*x2 + x0 + x2");

        let anf_1 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("0", 3).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "0");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "0");

        let anf_1 = AnfPolynomial::from_str("0", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "0");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "0");

        let anf_1 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("1", 3).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1 + x0 + x1 + x2");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1 + x0 + x1 + x2");

        let anf_1 = AnfPolynomial::from_str("1", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1 + x0 + x1 + x2");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1 + x0 + x1 + x2");

        let anf_1 = AnfPolynomial::from_str("x3*x7 + x3 + x4*x5 + x4*x6 + x5*x6*x7 + x6 + x7 + 1", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x3*x5*x6*x7 + x3*x6*x7 + x3*x6 + x3*x7 + x3 + x4*x5*x6 + x4*x5*x7 + x4*x7 + x4 + x6*x7 + x6 + x7 + 1", 8).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x3*x4*x5*x7 + x3*x4*x5 + x4*x5*x6 + x3*x4*x7 + x4*x5*x7 + x3*x6*x7 + x3*x4 + x3*x6 + x3*x7 + x4*x7 + x6*x7 + x3 + x4 + x6 + x7 + 1");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x3*x4*x5*x7 + x3*x4*x5 + x4*x5*x6 + x3*x4*x7 + x4*x5*x7 + x3*x6*x7 + x3*x4 + x3*x6 + x3*x7 + x4*x7 + x6*x7 + x3 + x4 + x6 + x7 + 1");

        let anf_1 = AnfPolynomial::from_str("x3*x5*x6*x7 + x3*x6*x7 + x3*x6 + x3*x7 + x3 + x4*x5*x6 + x4*x5*x7 + x4*x7 + x4 + x6*x7 + x6 + x7 + 1", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x3*x7 + x3 + x4*x5 + x4*x6 + x5*x6*x7 + x6 + x7 + 1", 8).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x3*x4*x5*x7 + x3*x4*x5 + x4*x5*x6 + x3*x4*x7 + x4*x5*x7 + x3*x6*x7 + x3*x4 + x3*x6 + x3*x7 + x4*x7 + x6*x7 + x3 + x4 + x6 + x7 + 1");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x3*x4*x5*x7 + x3*x4*x5 + x4*x5*x6 + x3*x4*x7 + x4*x5*x7 + x3*x6*x7 + x3*x4 + x3*x6 + x3*x7 + x4*x7 + x6*x7 + x3 + x4 + x6 + x7 + 1");

        let anf_1 = AnfPolynomial::from_str("0", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x3*x7 + x3 + x4*x5 + x4*x6 + x5*x6*x7 + x6 + x7 + 1", 8).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "0");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "0");

        let anf_1 = AnfPolynomial::from_str("1", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x3*x7 + x3 + x4*x5 + x4*x6 + x5*x6*x7 + x6 + x7 + 1", 8).unwrap();
        let anf_3 = anf_1.clone() & anf_2.clone();
        assert_eq!(anf_3.to_string(), "x5*x6*x7 + x4*x5 + x4*x6 + x3*x7 + x3 + x6 + x7 + 1");
        anf_2 &= anf_1;
        assert_eq!(anf_2.to_string(), "x5*x6*x7 + x4*x5 + x4*x6 + x3*x7 + x3 + x6 + x7 + 1");
    }

    #[test]
    fn test_mul() {
        let anf_1 = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x0*x1 + x0*x1*x2 + x1 + 1", 3).unwrap();
        let anf_3 = anf_1.clone() * anf_2.clone();
        assert_eq!(anf_3.to_string(), "x0*x1*x2 + x1*x2 + x0 + x2");
        anf_2 *= anf_1;
        assert_eq!(anf_2.to_string(), "x0*x1*x2 + x1*x2 + x0 + x2");

        let anf_1 = AnfPolynomial::from_str("x3*x7 + x3 + x4*x5 + x4*x6 + x5*x6*x7 + x6 + x7 + 1", 8).unwrap();
        let mut anf_2 = AnfPolynomial::from_str("x3*x5*x6*x7 + x3*x6*x7 + x3*x6 + x3*x7 + x3 + x4*x5*x6 + x4*x5*x7 + x4*x7 + x4 + x6*x7 + x6 + x7 + 1", 8).unwrap();
        let anf_3 = anf_1.clone() * anf_2.clone();
        assert_eq!(anf_3.to_string(), "x3*x4*x5*x7 + x3*x4*x5 + x4*x5*x6 + x3*x4*x7 + x4*x5*x7 + x3*x6*x7 + x3*x4 + x3*x6 + x3*x7 + x4*x7 + x6*x7 + x3 + x4 + x6 + x7 + 1");
        anf_2 *= anf_1;
        assert_eq!(anf_2.to_string(), "x3*x4*x5*x7 + x3*x4*x5 + x4*x5*x6 + x3*x4*x7 + x4*x5*x7 + x3*x6*x7 + x3*x4 + x3*x6 + x3*x7 + x4*x7 + x6*x7 + x3 + x4 + x6 + x7 + 1");
    }
}
