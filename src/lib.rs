//! # Boolean Function analysis library

#![doc = include_str!("../README.md")]
#![forbid(unsafe_code, unused_must_use)]
#![forbid(
    missing_docs,
    unreachable_pub,
    unused_import_braces,
    unused_extern_crates
)]

pub mod affine_equivalence_classes;
mod anf_polynom;
mod big_boolean_function;
mod boolean_function_error;
mod iterator;
mod small_boolean_function;
mod utils;

pub use crate::anf_polynom::AnfPolynomial;
use crate::boolean_function_error::{AND_DIFFERENT_VAR_COUNT_PANIC_MSG, XOR_DIFFERENT_VAR_COUNT_PANIC_MSG};
pub use crate::iterator::BooleanFunctionIterator;
use crate::BooleanFunctionError::{
    StringHexParseError, UnexpectedError, WrongStringHexTruthTableLength,
};
pub use big_boolean_function::BigBooleanFunction;
pub use boolean_function_error::BooleanFunctionError;
use enum_dispatch::enum_dispatch;
use gen_combinations::CombinationIterator;
use num_bigint::BigUint;
use num_traits::{Num, ToPrimitive};
pub use small_boolean_function::SmallBooleanFunction;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::ops::{Add, AddAssign, BitAnd, BitAndAssign, BitXor, BitXorAssign, Mul, MulAssign, Not};
use hackfn::hackfn;
use crate::iterator::CloseBalancedFunctionIterator;

/// Internal representation of Boolean function
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanFunctionType {
    /// Small boolean function with 6 or fewer variables, the truth table is stored in an u64
    Small,
    /// Big boolean function with more than 6 variables, the truth table is stored in a BigUint
    Big,
}

/// Trait for Boolean function implementations.
/// This trait is implemented by [SmallBooleanFunction] and [BigBooleanFunction].
///
/// You could use this trait via the [BooleanFunction] type, which encapsulates the [BooleanFunctionImpl] trait.
#[enum_dispatch]
pub trait BooleanFunctionImpl: Debug {
    /// Variable count of the Boolean function.
    fn variables_count(&self) -> usize;

    /// Internal type of the Boolean function abstraction:
    ///
    /// - [BooleanFunctionType::Big] for [BigBooleanFunction]
    ///
    /// - [BooleanFunctionType::Small] for [SmallBooleanFunction]
    fn get_boolean_function_type(&self) -> BooleanFunctionType;

    /// Maximum input value for the Boolean function, as unsigned 32-bit integer.
    ///
    /// This is equal to $2^n - 1$, where $n$ is the number of variables.
    fn get_max_input_value(&self) -> u32 {
        (1 << self.variables_count()) - 1
    }

    /// Returns `true` if the Boolean function is balanced, ie it has an equal number of 0 and 1 outputs.
    fn is_balanced(&self) -> bool;

    /// Computes the value of the Boolean function for a given input, as a 32-bit unsigned integer.
    ///
    /// This is equivalent of calling directly your Boolean function object as a function.
    ///
    /// # Parameters
    /// - `input_bits`: The input value for which to compute the Boolean function value, the least significant bit being the first variable.
    ///
    /// # Returns
    /// The value of the Boolean function for the given input bits.
    ///
    /// # Panics
    /// If the input value is greater than the maximum input value, and the `unsafe_disable_safety_checks` feature is not enabled.
    ///
    /// # Example
    /// ```rust
    /// use boolean_function::BooleanFunction;
    /// use boolean_function::BooleanFunctionImpl;
    ///
    /// let boolean_function = BooleanFunction::from_hex_string_truth_table("abce1234").unwrap();
    /// assert_eq!(boolean_function.compute_cellular_automata_rule(8), false);
    /// assert_eq!(boolean_function(8), false); // directly as a function
    /// ```
    fn compute_cellular_automata_rule(&self, input_bits: u32) -> bool;

    /// Computes the Walsh-Hadamard transform of the Boolean function for a given point.
    ///
    /// The Walsh-Hadamard transform of a Boolean function $f$, for a given point $\omega$, is defined as:
    ///
    /// $$W_f(\omega) = \sum_{x=0}^{2^n-1} (-1)^{f(x) \oplus \omega \cdot x}$$
    ///
    /// Where $\oplus$ is the XOR operation, $\cdot$ is the AND operand product, and $2^n - 1$ is the maximum input value.
    ///
    /// # Parameters
    /// - `w`: The point $\omega$ for which to compute the Walsh-Hadamard transform.
    ///
    /// # Returns
    /// The value of the Walsh-Hadamard transform for the given point.
    ///
    /// # Panics
    /// If the point is greater than the maximum input value, and the `unsafe_disable_safety_checks` feature is not enabled.
    fn walsh_hadamard_transform(&self, w: u32) -> i32 {
        let max_input_value = self.get_max_input_value();
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if w > max_input_value {
            panic!(
                "Too big Walsh parameter point, must be <= {}",
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

    /// Computes the Walsh-Hadamard values for all points.
    ///
    /// # Returns
    /// A vector containing the Walsh-Hadamard values for all points.
    fn walsh_hadamard_values(&self) -> Vec<i32> {
        (0..=self.get_max_input_value())
            .map(|w| self.walsh_hadamard_transform(w))
            .collect()
    }

    /// Computes the absolute Walsh-Hadamard spectrum of the Boolean function.
    ///
    /// The absolute Walsh-Hadamard spectrum is the number of occurrences of each absolute value of the Walsh-Hadamard transform.
    ///
    /// # Returns
    /// A hashmap containing the absolute Walsh-Hadamard values as keys, and the number of occurrences as values.
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

    /// Computes the Walsh-Fourier transform of the Boolean function for a given point.
    ///
    /// The Walsh-Fourier transform of a Boolean function $f$, for a given point $\omega$, is defined as:
    ///
    /// $$F_f(\omega) = \sum_{x=0}^{2^n-1} f(x) \cdot (-1)^{\omega \cdot x}$$
    ///
    /// Where $\cdot$ is the AND operand product, and $2^n - 1$ is the maximum input value.
    ///
    /// # Parameters
    /// - `w`: The point $\omega$ for which to compute the Walsh-Fourier transform.
    ///
    /// # Returns
    /// The value of the Walsh-Fourier transform for the given point.
    ///
    /// # Panics
    /// If the point is greater than the maximum input value, and the `unsafe_disable_safety_checks` feature is not enabled.
    fn walsh_fourier_transform(&self, w: u32) -> i32 {
        let max_input_value = self.get_max_input_value();
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if w > max_input_value {
            panic!(
                "Too big Walsh parameter point, must be <= {}",
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

    /// Computes the Walsh-Fourier values for all points.
    ///
    /// # Returns
    /// A vector containing the Walsh-Fourier values for all points.
    fn walsh_fourier_values(&self) -> Vec<i32> {
        (0..=self.get_max_input_value())
            .map(|w| self.walsh_fourier_transform(w))
            .collect()
    }

    /// Computes the autocorrelation transform of the Boolean function for a given point.
    /// The autocorrelation transform of a Boolean function $f$, for a given point $\omega$, is defined as:
    ///
    /// $$\Delta_f(\omega) = \sum_{x=0}^{2^n-1} (-1)^{f(x) \oplus f(x \oplus \omega)}$$
    ///
    /// Where $\oplus$ is the XOR operation, and $2^n - 1$ is the maximum input value.
    ///
    /// # Parameters
    /// - `w`: The point $\omega$ for which to compute the autocorrelation transform.
    ///
    /// # Returns
    /// The value of the autocorrelation transform for the given point.
    ///
    /// # Panics
    /// If the point is greater than the maximum input value, and the `unsafe_disable_safety_checks` feature is not enabled.
    fn auto_correlation_transform(&self, w: u32) -> i32 {
        let max_input_value = self.get_max_input_value();
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if w > max_input_value {
            panic!(
                "Too big auto-correlation parameter point, must be <= {}",
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

    /// Computes the absolute autocorrelation spectrum of the Boolean function.
    ///
    /// The absolute autocorrelation spectrum is the number of occurrences of each absolute value of the autocorrelation transform.
    ///
    /// # Returns
    /// A hashmap containing the absolute autocorrelation values as keys, and the number of occurrences as values.
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

    /// Returns the absolute indicator of the Boolean function.
    ///
    /// As defined [here](https://www.researchgate.net/publication/322383819_Distribution_of_the_absolute_indicator_of_random_Boolean_functions),
    /// the absolute indicator is the maximal absolute value of the auto-correlation transform for all points starting from 1 to the maximum input value $2^n - 1$,
    /// $n$ being the number of variables.
    fn absolute_indicator(&self) -> u32 {
        (1..=self.get_max_input_value())
            .map(|w| self.auto_correlation_transform(w).unsigned_abs())
            .max()
            .unwrap_or(0)
    }

    /// Computes the derivative of the Boolean function for a given direction.
    ///
    /// The derivative of a Boolean function $f$ in the direction $u$ is defined as:
    ///
    /// $$\nabla_u f(x) = f(x) \oplus f(x \oplus u)$$
    ///
    /// Where $\oplus$ is the XOR operation.
    ///
    /// # Parameters
    /// - `direction`: The direction $u$ for which to compute the derivative.
    ///
    /// # Returns
    /// The derivative of the Boolean function in the given direction, or an error if the direction is greater than the maximum input value and the `unsafe_disable_safety_checks` feature is not enabled.
    fn derivative(&self, direction: u32) -> Result<BooleanFunction, BooleanFunctionError>;

    /// Returns `true` if the Boolean function is linear, ie its algebraic degree is 0 or 1.
    fn is_linear(&self) -> bool;

    /// Reverse the Boolean function, ie swap 0 and 1 outputs.
    ///
    /// This is equivalent to the NOT operation, like `!boolean_function`.
    ///
    /// # Returns
    /// The reversed Boolean function.
    fn reverse(&self) -> BooleanFunction;

    /// Returns the Algebraic Normal Form of the Boolean function, in the form of an [AnfPolynomial].
    ///
    /// We use the Bakoev's algorithm to compute the ANF polynomial <http://www.math.bas.bg/moiuser/OCRT2017/a3.pdf>.
    ///
    /// # Returns
    /// The Algebraic Normal Form of the Boolean function.
    ///
    /// # Example
    /// ```rust
    /// // Wolfram's rule 30
    /// use boolean_function::BooleanFunctionImpl;
    /// use boolean_function::BooleanFunction;
    /// let boolean_function = BooleanFunction::from_u64_truth_table(30, 3).unwrap();
    /// let anf_polynomial = boolean_function.algebraic_normal_form();
    /// // `*` denotes the AND operation, and `+` denotes the XOR operation
    /// assert_eq!(anf_polynomial.to_string(), "x0*x1 + x0 + x1 + x2");
    /// ```
    fn algebraic_normal_form(&self) -> AnfPolynomial;

    /// Returns the algebraic degree of the Boolean function.
    ///
    /// The algebraic degree of a Boolean function is the maximum degree of the AND monomials in its Algebraic Normal Form.
    ///
    /// The degree of a monomial is the number of variables in the monomial.
    ///
    /// # Returns
    /// The algebraic degree of the Boolean function.
    fn algebraic_degree(&self) -> usize {
        self.algebraic_normal_form().get_degree()
    }

    /// Returns `true` if the Boolean function is symmetric.
    ///
    /// A Boolean function is symmetric if for all inputs $x$, the value of the function is the same as the value of the function for the bitwise reversed input.
    ///
    /// It means that the output depends only on the Hamming weight of the input.
    fn is_symmetric(&self) -> bool {
        let variables_count = self.variables_count() as u32;
        let precomputed_hamming_weights = (0..(variables_count + 1))
            .map(|i| self.compute_cellular_automata_rule(u32::MAX.checked_shr(32 - i).unwrap_or(0)))
            .collect::<Vec<bool>>();

        (0u32..(1 << variables_count)).all(|i| {
            precomputed_hamming_weights[i.count_ones() as usize]
                == self.compute_cellular_automata_rule(i)
        })
    }

    /// Returns the nonlinearity of the Boolean function.
    ///
    /// Nonlinearity is a measure of how far the Boolean function is from being linear.
    /// It's defined as the Hamming distance between the Boolean function and the closest linear function, or the minimum number of output bits that must be changed to make the function linear.
    ///
    /// # Returns
    /// The nonlinearity of the Boolean function, as an unsigned 32-bit integer.
    fn nonlinearity(&self) -> u32 {
        ((1 << self.variables_count())
            - (0..=self.get_max_input_value())
                .map(|x| self.walsh_hadamard_transform(x).unsigned_abs())
                .max()
                .unwrap_or(0))
            >> 1
    }

    /// Returns `true` if the Boolean function is bent, meaning maximally nonlinear, or perfectly nonlinear.
    ///
    /// Only functions with an even number of variables can be bent.
    ///
    /// The maximum possible nonlinearity for a Boolean function with $n$ variables is $\frac{2^n - 2^{\frac{n}{2}}}{2}$.
    ///
    /// A bent function has all its derivative balanced.
    ///
    /// The Rothaus theorem states that if an $n$-variable Boolean function $f$ is bent, then $deg(f) <= \frac{n}{2}$.
    ///
    /// A bent function cannot be balanced.
    fn is_bent(&self) -> bool {
        if self.variables_count() & 1 != 0 {
            return false;
        }
        self.nonlinearity()
            == ((1 << self.variables_count()) - (1 << (self.variables_count() >> 1))) >> 1
    }

    /// Returns `true` if the Boolean function is near-bent.
    ///
    /// A $n$-variable Boolean function is said to be *near-bent* if $n$ is odd,
    /// and its [Walsh-Hamamard](#method.walsh_hadamard_values) spectrum contains all and only the 3 values $\\{0, \pm 2^{\frac{n+1}{2}}\\}$.
    fn is_near_bent(&self) -> bool {
        if self.variables_count() & 1 == 0 {
            return false;
        }
        let absolute_walsh_allowed_value = 1 << ((self.variables_count() + 1) >> 1);
        let walsh_hadamard_spectrum = (0..=self.get_max_input_value())
            .map(|x| self.walsh_hadamard_transform(x))
            .collect::<HashSet<_>>();

        walsh_hadamard_spectrum.len() == 3 && walsh_hadamard_spectrum.iter().all(|w| {
            *w == 0 || *w == absolute_walsh_allowed_value || *w == -absolute_walsh_allowed_value
        })
    }

    /// Returns, if it exists, an annihilator function, its degree and the dimension of annihilator vector space.
    ///
    /// The annihilator of a Boolean function $f$ is a non-null Boolean function $g$ such that:
    ///
    /// $$f(x).g(x) = 0 \ \ \forall x \in \mathbb{F}_2^n$$
    ///
    /// **Special case**: annihilator of zero function is constant one function, by convention.
    ///
    /// # Parameters
    /// - `max_degree`: An optional maximum degree of the annihilator to search for. If set to `None`, the value is set to the variable count.
    ///
    /// # Returns
    /// - `None` if no annihilator is found (this is the case for constant one function).
    /// - `Some((annihilator, degree, dimension))` if an annihilator is found:
    ///    - `annihilator`: The annihilator function.
    ///    - `degree`: The degree of the returned annihilator function.
    ///    - `dimension`: The dimension of the annihilator vector space.
    fn annihilator(&self, max_degree: Option<usize>) -> Option<(BooleanFunction, usize, usize)>;

    /// Returns the algebraic immunity of the Boolean function.
    ///
    /// The algebraic immunity of a Boolean function is defined as the minimum degree of an annihilator of the function.
    ///
    /// # Returns
    /// The algebraic immunity of the Boolean function, or 0 if the function has no annihilator.
    fn algebraic_immunity(&self) -> usize {
        match self.annihilator(Some(self.variables_count())) {
            None => 0,
            Some(annihilator) => annihilator.1,
        }
    }

    /// Returns `true` if the Boolean function is plateaued.
    ///
    /// A Boolean function is plateaued if its Walsh-Hadamard spectrum has at most three values: 0, $\lambda$ and $-\lambda$, where $\lambda \in \mathbb{N}^*$.
    fn is_plateaued(&self) -> bool {
        let absolute_walsh_hadamard_spectrum = self.absolute_walsh_hadamard_spectrum();
        absolute_walsh_hadamard_spectrum.len() == 1
            || (absolute_walsh_hadamard_spectrum.len() == 2
                && absolute_walsh_hadamard_spectrum.contains_key(&0))
    }

    /// Returns the sum of the square of the indicator of the Boolean function.
    ///
    /// The sum of the square of the indicator for a $n$-variable Boolean function $f$ is defined as:
    ///
    /// $$\sum_{w=0}^{2^n-1} \Delta_f(w)^2$$
    ///
    /// Where $\Delta_f(w)$ is the [autocorrelation transform](Self::auto_correlation_transform) of the function for a given point $w$.
    ///
    /// # Returns
    /// The sum of the square of the indicator of the Boolean function.
    fn sum_of_square_indicator(&self) -> usize {
        (0..=self.get_max_input_value())
            .map(|w| self.auto_correlation_transform(w))
            .map(|value| (value as i64 * value as i64) as usize)
            .sum()
    }

    /// Returns `true` if the Boolean function has a linear structure.
    ///
    /// A $n$-variable boolean function has a linear structure if $\exists a \in \mathbb{F}_2^n$ such that
    /// $x \longmapsto f(x) \oplus f(x \oplus a)$ is a constant function.
    ///
    /// <https://www.sciencedirect.com/topics/mathematics/linear-structure>
    fn has_linear_structure(&self) -> bool {
        (1..=self.get_max_input_value()).any(|x| {
            self.auto_correlation_transform(x).unsigned_abs() == 1 << self.variables_count()
        })
    }

    /// Checks if the parameter is a linear structure of the Boolean function.
    ///
    /// A vector $a \in \mathbb{F}_2^n$ is a linear structure of a $n$-variable Boolean function $f$ if the function
    /// $x \longmapsto f(x) \oplus f(x \oplus a)$ is a constant function.
    ///
    /// # Parameters
    /// - `value`: The value to check if it is a linear structure.
    ///
    /// # Returns
    /// `true` if the value is a linear structure of the Boolean function, `false` otherwise.
    ///
    /// # Panics
    /// The function panics if the value is greater than the function maximum input value and the `unsafe_disable_safety_checks` feature is not enabled.
    fn is_linear_structure(&self, value: u32) -> bool {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        {
            let max_input_value = self.get_max_input_value();
            if value > max_input_value {
                panic!("Too big value parameter, must be <= {}", max_input_value);
            }
        }
        self.auto_correlation_transform(value).unsigned_abs() == 1 << self.variables_count()
    }

    /// List of all linear structures of the Boolean function.
    ///
    /// A vector $a \in \mathbb{F}_2^n$ is a linear structure of a $n$-variable Boolean function $f$ if the function
    /// $x \longmapsto f(x) \oplus f(x \oplus a)$ is a constant function.
    ///
    /// # Returns
    /// A vector containing all linear structures of the Boolean function.
    fn linear_structures(&self) -> Vec<u32> {
        (0..=self.get_max_input_value())
            .filter(|x| {
                self.auto_correlation_transform(*x).unsigned_abs() == 1 << self.variables_count()
            })
            .collect()
    }

    /// Returns the correlation immunity order of the Boolean function (calculated using Xiao Massey theorem).
    ///
    /// A Boolean function is said to be correlation immune of order $k$ if the output of the function is statistically independent of any subset of maximum $k$ input bits.
    ///
    /// As shown by Siegenthaler, a $n$-variable Boolean function  which is correlation immune of order $k$ has a degree $d \leq n - m$.
    /// <https://iacr.org/archive/asiacrypt2002/25010483/25010483.pdf>
    ///
    /// # Returns
    /// The correlation immunity order of the Boolean function.
    fn correlation_immunity(&self) -> usize {
        (1..=self.get_max_input_value())
            .filter(|x| self.walsh_hadamard_transform(*x) != 0)
            .map(|x| x.count_ones() as usize)
            .min()
            .unwrap_or(self.variables_count() + 1)
            - 1
    }

    /// Returns the resiliency order of the Boolean function.
    ///
    /// A boolean function is said to be resilient of order $k$ if it is balanced and correlation immune of order $k$.
    ///
    /// # Returns
    /// The resiliency order of the Boolean function, or `None` if the function is not balanced.
    fn resiliency_order(&self) -> Option<usize> {
        if !self.is_balanced() {
            return None;
        }
        Some(self.correlation_immunity())
    }

    /// Returns the support of the Boolean function.
    ///
    /// The support of a $n$-variable Boolean function is the set of all inputs $x \in \mathbb{F}_2^n$ such that $f(x) = 1$.
    ///
    /// # Returns
    /// The support of the Boolean function, as a set of unsigned 32-bit integers.
    fn support(&self) -> HashSet<u32> {
        (0..=self.get_max_input_value())
            .filter(|x| self.compute_cellular_automata_rule(*x))
            .collect()
    }

    /// Returns the maximum propagation criterion of the Boolean function.
    ///
    /// A $n$-variable Boolean function $f$ satisfies propagation criterion at order $k$ if its output changes with a probability
    /// of $\frac{1}{2}$ when changing the value of any subset of $i$ variables, $\forall i \in \llbracket 1, k \rrbracket$
    ///
    /// # Returns
    /// The maximum propagation criterion of the Boolean function.
    fn propagation_criterion(&self) -> usize {
        let num_variables = self.variables_count();
        let max_input_value = self.get_max_input_value();
        let possible_reversable_bit_position =
            (0..num_variables).into_iter().collect::<Vec<usize>>();
        (1..=num_variables)
            .into_iter()
            .take_while(|criterion_degree| {
                CombinationIterator::new(&possible_reversable_bit_position, *criterion_degree).all(
                    |combination| {
                        let mut bit_mask = 0;
                        for &bit_position in combination {
                            bit_mask |= (1 << bit_position) as u32;
                        }
                        let function_equal_mask_count = (0..=max_input_value)
                            .into_iter()
                            .filter(|&x| {
                                let x_prime = x ^ bit_mask;
                                self.compute_cellular_automata_rule(x)
                                    == self.compute_cellular_automata_rule(x_prime)
                            })
                            .count();
                        function_equal_mask_count == (1 << (num_variables - 1))
                    },
                )
            })
            .count()
    }

    /// Returns a 1-local neighbor of the Boolean function, at a specific position
    ///
    /// A 1-local neighbor of a Boolean function $f$ at position $i$ is a Boolean function $f_i$ such that:
    ///
    /// $$f_i(x) = \begin{cases}
    /// f(x) &\text{if } x \neq i \\\\
    ///    f(x) \oplus 1 &\text{if } x = i
    /// \end{cases}$$
    ///
    /// $f_i$ is said to be *connected* to $f$.
    ///
    /// # Parameters
    /// - `position`: The position $i$ at which to compute the 1-local neighbor.
    ///
    /// # Returns
    /// The 1-local neighbor of the Boolean function at the given position.
    ///
    /// # Panics
    /// If the position is greater than the maximum input value, and the `unsafe_disable_safety_checks` feature is not enabled.
    fn get_1_local_neighbor(&self, position: u32) -> BooleanFunction;

    /// Returns an iterator over the values of the Boolean function.
    ///
    /// # Returns
    /// An iterator over the values of the Boolean function.
    ///
    /// # Example
    /// ```rust
    /// // Wolfram's rule 30
    /// use boolean_function::BooleanFunctionImpl;
    /// use boolean_function::BooleanFunction;
    /// let boolean_function = BooleanFunction::from_u64_truth_table(30, 3).unwrap();
    /// let mut iterator = boolean_function.iter();
    /// assert_eq!(iterator.next(), Some(false));
    /// assert_eq!(iterator.next(), Some(true));
    /// ```
    fn iter(&self) -> BooleanFunctionIterator;

    /// Returns the truth table of the Boolean function as a hexadecimal string.
    ///
    /// # Returns
    /// The truth table of the Boolean function as a hexadecimal string.
    ///
    /// # Example
    /// ```rust
    /// // Wolfram's rule 30
    /// use boolean_function::BooleanFunctionImpl;
    /// use boolean_function::BooleanFunction;
    /// let boolean_function = BooleanFunction::from_u64_truth_table(30, 3).unwrap();
    /// assert_eq!(boolean_function.printable_hex_truth_table(), "1e");
    /// ```
    fn printable_hex_truth_table(&self) -> String;

    /// Returns the truth table of the Boolean function as a BigUint.
    ///
    /// # Returns
    /// The truth table of the Boolean function as a BigUint.
    fn biguint_truth_table(&self) -> BigUint;

    /// Returns the truth table of the Boolean function as an unsigned 64-bit integer, if it fits (meaning the Boolean function has 6 or less input variables).
    ///
    /// # Returns
    /// The truth table of the Boolean function as an unsigned 64-bit integer, or `None` if the truth table is too big to fit in an u64.
    fn try_u64_truth_table(&self) -> Option<u64> {
        self.biguint_truth_table().to_u64()
    }

    /// Returns an iterator over the closest possible balanced Boolean functions.
    ///
    /// "closest" means the minimum possible Hamming distance over the truth table.
    ///
    /// This is particularly useful if you want to extract a set of balanced functions from a bent function.
    /// So that you can generate highly nonlinear balanced functions.
    ///
    /// # Returns
    ///
    /// An iterator over close balanced Boolean functions, or an error if the function is already balanced.
    ///
    /// # Note
    /// It is assumed that the function truth table is correctly sanitized, so be careful if you generated it with `unsafe_disable_safety_checks` feature activated.
    ///
    /// # Example
    /// ```rust
    /// use boolean_function::BooleanFunction;
    ///  use crate::boolean_function::BooleanFunctionImpl;
    ///
    /// let bent_function = BooleanFunction::from_hex_string_truth_table(
    ///             "80329780469d0b85cd2ad63e1a6ba42adbd83c9a0c55e4e8c99f227b0ffc1418"
    ///         ).unwrap();
    /// let close_balanced_iterator = bent_function.close_balanced_functions_iterator();
    /// assert!(close_balanced_iterator.is_ok());
    /// let mut close_balanced_iterator = close_balanced_iterator.unwrap();
    /// let highly_nonlinear_balanced_function = close_balanced_iterator.next().unwrap();
    /// assert!(highly_nonlinear_balanced_function.is_balanced());
    /// assert_eq!(highly_nonlinear_balanced_function.nonlinearity(), 112);
    /// ```
    fn close_balanced_functions_iterator(&self) -> Result<CloseBalancedFunctionIterator, BooleanFunctionError>;
}

/// This type is used to store a boolean function with any number of variables.
///
/// It abstracts The [SmallBooleanFunction] and [BigBooleanFunction] types, by encapsulating the [BooleanFunctionImpl] trait.
///
/// Please refer to the [BooleanFunctionImpl] trait for more information.
#[enum_dispatch(BooleanFunctionImpl)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BooleanFunction {
    /// Struct representing a boolean function with a small truth table.
    Small(pub SmallBooleanFunction),
    /// Struct representing a boolean function with a big truth table.
    Big(pub BigBooleanFunction),
}

/// In-place XOR operator for Boolean functions truth tables.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitXorAssign for BooleanFunction {
    fn bitxor_assign(&mut self, rhs: Self) {
        let self_num_variables = self.variables_count();
        let rhs_num_variables = rhs.variables_count();
        if self_num_variables != rhs_num_variables {
            panic!("{}", XOR_DIFFERENT_VAR_COUNT_PANIC_MSG);
        }
        match (
            self,
            rhs,
        ) {
            (BooleanFunction::Small(self_small_boolean_function), BooleanFunction::Small(rhs_small_boolean_function)) => {
                *self_small_boolean_function ^= rhs_small_boolean_function;
            },
            (BooleanFunction::Small(self_small_boolean_function), BooleanFunction::Big(rhs_big_boolean_function)) => {
                let rhs_small_boolean_function = SmallBooleanFunction::from_truth_table(
                rhs_big_boolean_function
                        .biguint_truth_table()
                        .to_u64()
                        .unwrap(),
                    self_num_variables,
                )
                .unwrap();
                *self_small_boolean_function ^= rhs_small_boolean_function;
            },
            (BooleanFunction::Big(self_big_boolean_function), BooleanFunction::Small(rhs_small_boolean_function)) => {
                let rhs_big_boolean_function = BigBooleanFunction::from_truth_table(
                    rhs_small_boolean_function.biguint_truth_table(),
                    rhs_num_variables,
                );
                *self_big_boolean_function ^= rhs_big_boolean_function;
            },
            (BooleanFunction::Big(self_big_boolean_function), BooleanFunction::Big(rhs_big_boolean_function)) => {
                *self_big_boolean_function ^= rhs_big_boolean_function;
            }
        }
    }
}

/// XOR operator for Boolean functions truth tables.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitXor for BooleanFunction {
    type Output = Self;

    fn bitxor(mut self, rhs: Self) -> Self::Output {
        self ^= rhs;
        self
    }
}

/// ADD operator for Boolean functions truth tables.
///
/// It is equivalent to [crate::BooleanFunction::bitxor] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl Add for BooleanFunction {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self ^ rhs
    }
}

/// In-place ADD operator for Boolean functions truth tables.
///
/// It is equivalent to [crate::BooleanFunction::bitxor_assign] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl AddAssign for BooleanFunction {
    fn add_assign(&mut self, rhs: Self) {
        *self ^= rhs;
    }
}

/// In-place AND operator for Boolean functions truth tables.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitAndAssign for BooleanFunction {
    fn bitand_assign(&mut self, rhs: Self) {
        let self_num_variables = self.variables_count();
        let rhs_num_variables = rhs.variables_count();
        if self_num_variables != rhs_num_variables {
            panic!("{}", AND_DIFFERENT_VAR_COUNT_PANIC_MSG);
        }
        match (
            self,
            rhs,
        ) {
            (BooleanFunction::Small(self_small_boolean_function), BooleanFunction::Small(rhs_small_boolean_function)) => {
                *self_small_boolean_function &= rhs_small_boolean_function;
            },
            (BooleanFunction::Small(self_small_boolean_function), BooleanFunction::Big(rhs_big_boolean_function)) => {
                let rhs_small_boolean_function = SmallBooleanFunction::from_truth_table(
                    rhs_big_boolean_function
                        .biguint_truth_table()
                        .to_u64()
                        .unwrap(),
                    self_num_variables,
                )
                    .unwrap();
                *self_small_boolean_function &= rhs_small_boolean_function;
            },
            (BooleanFunction::Big(self_big_boolean_function), BooleanFunction::Small(rhs_small_boolean_function)) => {
                let rhs_big_boolean_function = BigBooleanFunction::from_truth_table(
                    rhs_small_boolean_function.biguint_truth_table(),
                    rhs_num_variables,
                );
                *self_big_boolean_function &= rhs_big_boolean_function;
            },
            (BooleanFunction::Big(self_big_boolean_function), BooleanFunction::Big(rhs_big_boolean_function)) => {
                *self_big_boolean_function &= rhs_big_boolean_function;
            }
        }
    }
}

/// AND operator for Boolean functions truth tables.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl BitAnd for BooleanFunction {
    type Output = Self;

    fn bitand(mut self, rhs: Self) -> Self::Output {
        self &= rhs;
        self
    }
}

/// MUL operator for Boolean functions truth tables.
///
/// It is equivalent to [crate::BooleanFunction::bitand] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl Mul for BooleanFunction {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        self & rhs
    }
}

/// In-place MUL operator for Boolean functions truth tables.
///
/// It is equivalent to [crate::BooleanFunction::bitand_assign] operator.
///
/// # Panics
/// If the Boolean functions have different number of variables, and the `unsafe_disable_safety_checks` feature is not enabled.
impl MulAssign for BooleanFunction {
    fn mul_assign(&mut self, rhs: Self) {
        *self &= rhs;
    }
}

/// NOT operator for Boolean functions.
///
/// This is equivalent to the [BooleanFunctionImpl::reverse] operation: it reverses each output of the Boolean function.
impl Not for BooleanFunction {
    type Output = Self;

    fn not(self) -> Self::Output {
        self.reverse()
    }
}

/// Create Boolean Function from hex string truth table representation
impl TryFrom<&str> for BooleanFunction {
    type Error = BooleanFunctionError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Self::from_hex_string_truth_table(value)
    }
}

#[hackfn]
impl BooleanFunction {
    fn call(&self, input_bits: u32) -> bool {
        self.compute_cellular_automata_rule(input_bits)
    }
}

impl BooleanFunction {
    /// Creates a new BooleanFunction from a hexadecimal string representing the truth table.
    ///
    /// The hexadecimal string must have a length of $\frac{2^n}{4}$ where $n$ is the number of variables of the Boolean function
    /// (meaning this function only accepts Boolean function with 2 or more input variables).
    ///
    /// # Parameters
    /// - `hex_truth_table`: The hexadecimal string representing the truth table of the Boolean function.
    ///
    /// # Returns
    /// The Boolean function created from the hexadecimal string, or an error if the string cannot be parsed as a hexadecimal truth table of a Boolean function.
    ///
    /// # Example
    /// ```rust
    /// use boolean_function::BooleanFunction;
    /// let boolean_function = BooleanFunction::from_hex_string_truth_table("0969817CC5893BA6")
    ///     .unwrap();
    /// ```
    pub fn from_hex_string_truth_table(
        hex_truth_table: &str,
    ) -> Result<Self, BooleanFunctionError> {
        if hex_truth_table.len().count_ones() != 1 {
            return Err(WrongStringHexTruthTableLength);
        }
        let num_variables = (hex_truth_table.len() << 2).trailing_zeros() as usize;
        if num_variables <= 6 {
            Ok(SmallBooleanFunction::from_truth_table(
                u64::from_str_radix(hex_truth_table, 16).map_err(|_| StringHexParseError)?,
                num_variables,
            )
                .map_err(|_| UnexpectedError)?
                .into())
        } else {
            Ok(BigBooleanFunction::from_truth_table(
                BigUint::from_str_radix(hex_truth_table, 16).map_err(|_| StringHexParseError)?,
                num_variables,
            )
                .into())
        }
    }

    /// Creates a new BooleanFunction from a list of integers representing the [Walsh-Hadamard transform](BooleanFunctionImpl::walsh_hadamard_values), by applying a reverse Walsh-Hadamard transform.
    ///
    /// The Walsh-Hadamard values list for an $n$-variable Boolean function must have a length of $2^n$.
    /// The function returns a [SmallBooleanFunction] if the given Walsh-Hadamard values list has a length of 64 or less, and a [BigBooleanFunction] otherwise.
    /// It won't check that the values are consistent, if not the function will return the closest Boolean function possible.
    ///
    /// # Parameters
    /// - `walsh_hadamard_values`: The list of integers representing the Walsh-Hadamard transform of the Boolean function.
    ///
    /// # Returns
    /// The Boolean function created from the Walsh-Hadamard values list, or an error if the list length is less than 4, or not a power of 2.
    pub fn from_reverse_walsh_hadamard_transform(
        walsh_hadamard_values: &[i32],
    ) -> Result<Self, BooleanFunctionError> {
        const MAX_WALSH_VALUES_SMALL: usize = 64; // (2^6)
        if walsh_hadamard_values.len() > MAX_WALSH_VALUES_SMALL {
            return Ok(BigBooleanFunction::from_walsh_hadamard_values(walsh_hadamard_values)?.into());
            // Error is handled in BigBooleanFunction constructor
        }
        Ok(SmallBooleanFunction::from_walsh_hadamard_values(walsh_hadamard_values)?.into())
    }

    /// Creates a new BooleanFunction from a list of integers representing the [Walsh-Fourier transform](BooleanFunctionImpl::walsh_fourier_values), by applying a reverse Walsh-Fourier transform.
    ///
    /// The Walsh-Fourier values list for an $n$-variable Boolean function must have a length of $2^n$.
    /// The function returns a [SmallBooleanFunction] if the given Walsh-Fourier values list has a length of 64 or less, and a [BigBooleanFunction] otherwise.
    /// It won't check that the values are consistent, if not the function will return the closest Boolean function possible.
    ///
    /// # Parameters
    /// - `walsh_fourier_values`: The list of integers representing the Walsh-Fourier transform of the Boolean function.
    ///
    /// # Returns
    /// The Boolean function created from the Walsh-Fourier values list, or an error if the list length is less than 4, or not a power of 2.
    pub fn from_reverse_walsh_fourier_transform(
        walsh_fourier_values: &[i32],
    ) -> Result<Self, BooleanFunctionError> {
        const MAX_WALSH_VALUES_SMALL: usize = 64; // (2^6)
        if walsh_fourier_values.len() > MAX_WALSH_VALUES_SMALL {
            return Ok(BigBooleanFunction::from_walsh_fourier_values(walsh_fourier_values)?.into());
            // Error is handled in BigBooleanFunction constructor
        }
        Ok(SmallBooleanFunction::from_walsh_fourier_values(walsh_fourier_values)?.into())
    }

    /// Creates a new BooleanFunction from an u64 representing the truth table (meaning the Boolean function has 6 or less input variables).
    ///
    /// The wrapped Boolean function is a [SmallBooleanFunction].
    ///
    /// # Parameters
    /// - `truth_table`: The u64 truth table of the Boolean function, where the lower bit represents the output of the Boolean function for the input 0.
    /// - `num_variables`: The number of input variables of the Boolean function.
    ///
    /// # Returns
    /// The Boolean function created from the u64 truth table.
    ///
    /// Returns an error if:
    /// - The given input variables count is greater than 6.
    /// - The given truth table is too big for the given input variables count and the `unsafe_disable_safety_checks` feature is not enabled.
    pub fn from_u64_truth_table(
        truth_table: u64,
        num_variables: usize,
    ) -> Result<Self, BooleanFunctionError> {
        Ok(SmallBooleanFunction::from_truth_table(truth_table, num_variables)?.into())
    }

    /// Creates a new BooleanFunction from a BigUint representing the truth table.
    ///
    /// The wrapped Boolean function is a [BigBooleanFunction] if the variables count is greater than 6, and a [SmallBooleanFunction] otherwise.
    ///
    /// # Parameters
    /// - `truth_table`: The BigUint truth table of the Boolean function, where the lower bit represents the output of the Boolean function for the input 0.
    /// - `num_variables`: The number of input variables of the Boolean function.
    ///
    /// # Returns
    /// The Boolean function created from the BigUint truth table.
    ///
    /// Returns an error if:
    /// - The given input variables count is greater than 31, and the `unsafe_disable_safety_checks` feature is not enabled.
    /// - The given truth table is too big for the given input variables count and the `unsafe_disable_safety_checks` feature is not enabled.
    pub fn from_biguint_truth_table(
        truth_table: &BigUint,
        num_variables: usize,
    ) -> Result<Self, BooleanFunctionError> {
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if truth_table.bits() > (1 << num_variables) {
            return Err(BooleanFunctionError::TooBigTruthTableForVarCount);
        }
        #[cfg(not(feature = "unsafe_disable_safety_checks"))]
        if num_variables > 31 {
            return Err(BooleanFunctionError::TooBigVariableCount(31));
        }
        if num_variables <= 6 {
            return Ok(SmallBooleanFunction::from_truth_table(
                truth_table.to_u64().unwrap(),
                num_variables,
            )?
                .into());
        }
        Ok(BigBooleanFunction::from_truth_table(truth_table.clone(), num_variables).into())
    }

    /// Computes Boolean Function from string ANF representation
    ///
    /// The ANF string representation must be in exact form "`x0*x2*x3 + x2*x3 + x1 + 1`".
    ///
    /// X's index starts at 0, meaning the maximum index is variable count - 1.
    ///
    /// # Parameters:
    /// - `anf_polynomial`: The string representation of the ANF form
    /// - `num_variables`: Variable count of the polynomial
    /// 
    /// # Returns
    /// The BooleanFunction corresponding to the ANF string representation, or an error if the input string doesn't respect the format and `unsafe_disable_safety_checks` feature is not activated.
    ///
    /// # Example
    /// ```rust
    /// use boolean_function::{BooleanFunction, BooleanFunctionImpl};
    ///
    /// let rule_30 = BooleanFunction::from_anf_polynomial_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
    /// assert_eq!(rule_30.printable_hex_truth_table(), "1e");
    /// ```
    pub fn from_anf_polynomial_str(anf_polynomial: &str, num_variables: usize) -> Result<Self, BooleanFunctionError> {
        Ok(AnfPolynomial::from_str(anf_polynomial, num_variables)?.to_boolean_function())
    }

    /// Computes Boolean Function from ANF polynomial
    ///
    /// # Parameters:
    /// - `anf_polynomial`: The polynomial in Algebraic Normal Form
    ///
    /// # Returns
    /// The BooleanFunction corresponding to the ANF polynomial
    pub fn from_anf_polynomial(anf_polynomial: &AnfPolynomial) -> Self {
        anf_polynomial.to_boolean_function()
    }
}

#[cfg(test)]
mod tests {
    use crate::BooleanFunctionError::InvalidWalshValuesCount;
    use crate::{AnfPolynomial, BooleanFunction, BooleanFunctionImpl, BooleanFunctionType};
    use num_bigint::BigUint;
    use num_traits::Num;
    use std::collections::HashMap;

    #[test]
    fn test_boolean_function_from_hex_string_truth_table() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD");
        assert!(boolean_function.is_err());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("");
        assert!(boolean_function.is_err());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("fe1z");
        assert!(boolean_function.is_err());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.variables_count(), 4);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.variables_count(), 7);
    }

    #[test]
    fn test_variables_count() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.variables_count(), 7);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.variables_count(), 4);
    }

    #[test]
    fn test_get_boolean_function_type() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6").unwrap();
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
    }

    #[test]
    fn test_get_max_input_value() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.get_max_input_value(), 127);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.get_max_input_value(), 15);
    }

    #[test]
    fn test_is_balanced() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert!(boolean_function.is_balanced());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD1")
                .unwrap();
        assert!(!boolean_function.is_balanced());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("aa55aa55").unwrap();
        assert!(boolean_function.is_balanced());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abce1234").unwrap();
        assert!(!boolean_function.is_balanced());
    }

    #[test]
    fn test_compute_cellular_automata_rule() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(boolean_function.compute_cellular_automata_rule(0), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(1), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(4), true);
        assert_eq!(boolean_function.compute_cellular_automata_rule(8), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(23), true);

        assert_eq!(boolean_function(0), false);
        assert_eq!(boolean_function(1), false);
        assert_eq!(boolean_function(4), true);
        assert_eq!(boolean_function(8), false);
        assert_eq!(boolean_function(23), true);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.compute_cellular_automata_rule(13), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(62), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(64), false);
        assert_eq!(boolean_function.compute_cellular_automata_rule(80), true);
        assert_eq!(boolean_function.compute_cellular_automata_rule(100), true);

        assert_eq!(boolean_function(13), false);
        assert_eq!(boolean_function(62), false);
        assert_eq!(boolean_function(64), false);
        assert_eq!(boolean_function(80), true);
        assert_eq!(boolean_function(100), true);
    }

    #[test]
    fn test_walsh_hadamard_transform() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.walsh_hadamard_transform(0), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(1), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(7), -16);
        assert_eq!(boolean_function.walsh_hadamard_transform(15), 16);
        assert_eq!(boolean_function.walsh_hadamard_transform(126), 16);
        assert_eq!(boolean_function.walsh_hadamard_transform(127), -16);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("aa55aa55").unwrap();
        assert_eq!(boolean_function.walsh_hadamard_transform(0), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(1), 0);
        assert_eq!(boolean_function.walsh_hadamard_transform(9), -32);
        assert_eq!(boolean_function.walsh_hadamard_transform(31), 0);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(boolean_function.walsh_hadamard_transform(0), 2);
        assert_eq!(boolean_function.walsh_hadamard_transform(1), 6);
        assert_eq!(boolean_function.walsh_hadamard_transform(2), -2);
        assert_eq!(boolean_function.walsh_hadamard_transform(31), -6);
    }

    #[test]
    fn test_absolute_walsh_hadamard_spectrum() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.absolute_walsh_hadamard_spectrum(),
            HashMap::from([(0, 64), (16, 64)])
        );

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(
            boolean_function.absolute_walsh_hadamard_spectrum(),
            HashMap::from([(6, 10), (10, 6), (2, 16)])
        );
    }

    #[test]
    fn test_auto_correlation_transform() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.auto_correlation_transform(0), 128);
        assert_eq!(boolean_function.auto_correlation_transform(1), -24);
        assert_eq!(boolean_function.auto_correlation_transform(126), -8);
        assert_eq!(boolean_function.auto_correlation_transform(127), -32);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("03").unwrap();
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
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.absolute_autocorrelation(),
            HashMap::from([(0, 33), (8, 58), (16, 28), (24, 6), (32, 2), (128, 1)])
        );

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abce1234").unwrap();
        assert_eq!(
            boolean_function.absolute_autocorrelation(),
            HashMap::from([(4, 25), (12, 6), (32, 1)])
        );
    }

    #[test]
    fn test_derivative() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("aa55aa55").unwrap();
        let derivative = boolean_function.derivative(1).unwrap();
        assert_eq!(derivative.variables_count(), 5);
        assert_eq!(derivative.printable_hex_truth_table(), "ffffffff");

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let derivative = boolean_function.derivative(1).unwrap();
        assert_eq!(derivative.variables_count(), 7);
        assert_eq!(
            derivative.printable_hex_truth_table(),
            "cfffc3c00fcf0cfff003f3ccf3f0ff30"
        );
    }

    #[test]
    fn test_algebraic_normal_form() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2*x3 + x0*x1*x2*x4 + x0*x1*x3*x4 + x1*x2*x3*x4 + x0*x1*x2*x5 + x0*x2*x3*x5 + x1*x2*x3*x5 + x0*x2*x4*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5 + x0*x1*x2*x6 + x0*x1*x3*x6 + x0*x2*x3*x6 + x1*x2*x3*x6 + x0*x2*x4*x6 + x1*x2*x4*x6 + x0*x3*x4*x6 + x2*x3*x4*x6 + x0*x1*x5*x6 + x0*x2*x5*x6 + x1*x2*x5*x6 + x1*x3*x5*x6 + x2*x3*x5*x6 + x3*x4*x5*x6 + x0*x1*x2 + x0*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 + x0*x1*x5 + x0*x2*x5 + x1*x2*x5 + x1*x3*x5 + x2*x3*x5 + x0*x4*x5 + x2*x4*x5 + x3*x4*x5 + x0*x2*x6 + x1*x2*x6 + x3*x4*x6 + x0*x5*x6 + x2*x5*x6 + x3*x5*x6 + x0*x2 + x0*x3 + x2*x4 + x3*x5 + x0*x6 + x1*x6 + x2*x6 + x3*x6 + x5*x6 + x2 + x4 + x5");

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD1")
                .unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "x0*x1*x2*x3*x4*x5*x6 + x0*x1*x2*x3*x4*x5 + x0*x1*x2*x3*x4*x6 + x0*x1*x2*x3*x5*x6 + x0*x1*x2*x4*x5*x6 + x0*x1*x3*x4*x5*x6 + x0*x2*x3*x4*x5*x6 + x1*x2*x3*x4*x5*x6 + x0*x1*x2*x3*x4 + x0*x1*x2*x3*x5 + x0*x1*x2*x4*x5 + x0*x1*x3*x4*x5 + x0*x2*x3*x4*x5 + x1*x2*x3*x4*x5 + x0*x1*x2*x3*x6 + x0*x1*x2*x4*x6 + x0*x1*x3*x4*x6 + x0*x2*x3*x4*x6 + x1*x2*x3*x4*x6 + x0*x1*x2*x5*x6 + x0*x1*x3*x5*x6 + x0*x2*x3*x5*x6 + x1*x2*x3*x5*x6 + x0*x1*x4*x5*x6 + x0*x2*x4*x5*x6 + x1*x2*x4*x5*x6 + x0*x3*x4*x5*x6 + x1*x3*x4*x5*x6 + x2*x3*x4*x5*x6 + x0*x2*x3*x4 + x0*x1*x3*x5 + x0*x1*x4*x5 + x0*x3*x4*x5 + x0*x1*x4*x6 + x1*x3*x4*x6 + x0*x3*x5*x6 + x0*x4*x5*x6 + x1*x4*x5*x6 + x2*x4*x5*x6 + x0*x1*x3 + x1*x2*x3 + x0*x1*x4 + x0*x2*x4 + x0*x3*x4 + x0*x3*x5 + x1*x4*x5 + x0*x1*x6 + x0*x3*x6 + x1*x3*x6 + x2*x3*x6 + x0*x4*x6 + x1*x4*x6 + x2*x4*x6 + x1*x5*x6 + x4*x5*x6 + x0*x1 + x1*x2 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x3*x4 + x0*x5 + x1*x5 + x2*x5 + x4*x5 + x4*x6 + x0 + x1 + x3 + x6 + 1");

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(
            boolean_function.algebraic_normal_form().to_string(),
            "x0*x1 + x0 + x1 + x2"
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("00").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "0");

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ff").unwrap();
        assert_eq!(boolean_function.algebraic_normal_form().to_string(), "1");

        let boolean_function = BooleanFunction::from_hex_string_truth_table("6e").unwrap();
        assert_eq!(
            boolean_function.algebraic_normal_form().to_string(),
            "x0*x1*x2 + x0*x1 + x0 + x1"
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("7b").unwrap();
        assert_eq!(
            boolean_function.algebraic_normal_form().to_string(),
            "x0*x1 + x1*x2 + x1 + 1"
        );
    }

    #[test]
    fn test_algebraic_degree() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 4);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD1")
                .unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 7);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000000000000000000000000000")
                .unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("ffffffffffffffffffffffffffffffff")
                .unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 2);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("00").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("0f").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 1);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ff").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("6e").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 3);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("7b").unwrap();
        assert_eq!(boolean_function.algebraic_degree(), 2);
    }

    #[test]
    fn test_printable_hex_truth_table() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0069817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "0069817cc5893ba6ac326e47619f5ad0"
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "fe12");
    }

    #[test]
    fn test_clone() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let cloned_boolean_function = boolean_function.clone();
        assert_eq!(&boolean_function, &cloned_boolean_function);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("fe12").unwrap();
        let cloned_boolean_function = boolean_function.clone();
        assert_eq!(&boolean_function, &cloned_boolean_function);
    }

    #[test]
    fn test_eq() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(&boolean_function, &boolean_function);

        let boolean_function2 =
            BooleanFunction::from_hex_string_truth_table("fe12").unwrap();
        assert_eq!(&boolean_function2, &boolean_function2);

        assert_ne!(&boolean_function2, &boolean_function);

        let boolean_function3 =
            BooleanFunction::from_hex_string_truth_table("0000fe12").unwrap();
        assert_ne!(&boolean_function3, &boolean_function2);

        let boolean_function4 =
            BooleanFunction::from_hex_string_truth_table("0969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_ne!(&boolean_function, &boolean_function4);
    }

    #[test]
    fn test_reverse() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let reversed_boolean_function = boolean_function.reverse();
        assert_eq!(reversed_boolean_function.variables_count(), 7);
        assert_eq!(
            reversed_boolean_function.printable_hex_truth_table(),
            "86967e833a76c45953cd91b89e60a52f"
        );

        let reversed_boolean_function = !boolean_function;
        assert_eq!(reversed_boolean_function.variables_count(), 7);
        assert_eq!(
            reversed_boolean_function.printable_hex_truth_table(),
            "86967e833a76c45953cd91b89e60a52f"
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("fe12").unwrap();
        let reversed_boolean_function = boolean_function.reverse();
        assert_eq!(reversed_boolean_function.variables_count(), 4);
        assert_eq!(
            reversed_boolean_function.printable_hex_truth_table(),
            "01ed"
        );

        let reversed_boolean_function = !boolean_function;
        assert_eq!(reversed_boolean_function.variables_count(), 4);
        assert_eq!(
            reversed_boolean_function.printable_hex_truth_table(),
            "01ed"
        );
    }

    #[test]
    fn test_is_linear() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert!(!boolean_function.is_linear());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0000000000000000ffffffffffffffff")
                .unwrap();
        assert!(boolean_function.is_linear());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abcdef0123456789").unwrap();
        assert!(!boolean_function.is_linear());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000ffffffff").unwrap();
        assert!(boolean_function.is_linear());
    }

    #[test]
    fn test_is_symmetric() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("00").unwrap();
        assert!(boolean_function.is_symmetric());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ff").unwrap();
        assert!(boolean_function.is_symmetric());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("80").unwrap();
        assert!(boolean_function.is_symmetric());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert!(!boolean_function.is_symmetric());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000008").unwrap();
        assert!(!boolean_function.is_symmetric());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("ffffffffffffffffffffffffffffffff")
                .unwrap();
        assert!(boolean_function.is_symmetric());
    }

    #[test]
    fn test_nonlinearity() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 0);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 0);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0000000a").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 2);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert_eq!(boolean_function.nonlinearity(), 28);
    }

    #[test]
    fn test_is_bent() {
        let boolean_function =
            BooleanFunction::try_from("00000000").unwrap();
        assert!(!boolean_function.is_bent());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert!(boolean_function.is_bent());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000ffffffff").unwrap();
        assert!(!boolean_function.is_bent());
    }

    #[test]
    fn test_is_near_bent() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("f9")
            .unwrap();
        assert!(boolean_function.is_near_bent());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ff")
            .unwrap();
        assert!(!boolean_function.is_near_bent());

        // even variable count -> cannot be near-bent
        let boolean_function = BooleanFunction::from_hex_string_truth_table("f9f9")
            .unwrap();
        assert!(!boolean_function.is_near_bent());
    }

    #[test]
    fn test_annihilator() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000").unwrap();
        let annihilator = boolean_function.annihilator(Some(0)).unwrap();
        assert_eq!(annihilator.0.printable_hex_truth_table(), "ffffffff");
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 1);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abcdef0123456789").unwrap();
        let annihilator = boolean_function.annihilator(Some(4)).unwrap();
        assert_eq!(
            annihilator.0.printable_hex_truth_table(),
            "1010101010101010"
        );
        assert_eq!(annihilator.1, 3);
        assert_eq!(annihilator.2, 25);

        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
        )
        .unwrap();
        let annihilator = boolean_function.annihilator(Some(4));
        assert!(annihilator.is_none());

        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();
        let annihilator = boolean_function.annihilator(Some(4)).unwrap();
        assert_eq!(
            annihilator.0.printable_hex_truth_table(),
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        );
        assert_eq!(annihilator.1, 0);
        assert_eq!(annihilator.2, 163);

        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "80921c010276c44224422441188118822442244118811880400810a80e200425",
        )
        .unwrap();
        let annihilator = boolean_function.annihilator(Some(1));
        assert!(annihilator.is_none());
        let annihilator = boolean_function.annihilator(Some(5)).unwrap();
        assert_eq!(
            annihilator.0.printable_hex_truth_table(),
            "2244224411881188d2b4d2b4e178e178d2b4d2b4e178e1782244224411881188"
        );
        assert_eq!(annihilator.1, 2);
        assert_eq!(annihilator.2, 155);
    }

    #[test]
    fn test_algebraic_immunity() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
        )
        .unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "0000000000000000000000000000000000000000000000000000000000000000",
        )
        .unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("0000").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "80921c010276c44224422441188118822442244118811880400810a80e200425",
        )
        .unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);

        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "2244224411881188d2b4d2b4e178e178d2b4d2b4e178e1782244224411881188",
        )
        .unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 3);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.algebraic_immunity(), 2);
    }

    #[test]
    fn test_is_plateaued() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
        )
        .unwrap();
        assert!(boolean_function.is_plateaued());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("8778").unwrap();
        assert!(boolean_function.is_plateaued());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert!(boolean_function.is_plateaued());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abcdef1234567890").unwrap();
        assert!(!boolean_function.is_plateaued());
    }

    #[test]
    fn test_sum_of_square_indicator() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 32768);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("0000").unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 4096);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abcdef1234567890abcdef1234567890")
                .unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 84992);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.sum_of_square_indicator(), 32768);
    }

    #[test]
    fn test_absolute_indicator() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 8);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 32);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("0000").unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 16);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abcdef1234567890abcdef1234567890")
                .unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 128);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 32);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.absolute_indicator(), 32);
    }

    #[test]
    fn test_linear_structures() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0, 4]);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("abcdef1234567890abcdef1234567890")
                .unwrap();
        assert_eq!(boolean_function.linear_structures(), [0, 64]);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert_eq!(boolean_function.linear_structures(), [0]);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000ffffffff").unwrap();
        assert_eq!(
            boolean_function.linear_structures(),
            [
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
                44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63
            ]
        );

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000").unwrap();
        assert_eq!(
            boolean_function.linear_structures(),
            [
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                23, 24, 25, 26, 27, 28, 29, 30, 31
            ]
        );

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("ffffffff").unwrap();
        assert_eq!(
            boolean_function.linear_structures(),
            [
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                23, 24, 25, 26, 27, 28, 29, 30, 31
            ]
        );
    }

    #[test]
    fn test_has_linear_structure() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("659a").unwrap();
        assert!(boolean_function.has_linear_structure());

        let boolean_function = BooleanFunction::from_hex_string_truth_table("dd0e").unwrap();
        assert!(!boolean_function.has_linear_structure());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000").unwrap();
        assert!(boolean_function.has_linear_structure());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("ffffffff").unwrap();
        assert!(boolean_function.has_linear_structure());

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("0113077C165E76A8").unwrap();
        assert!(!boolean_function.has_linear_structure());
    }

    #[test]
    fn test_is_linear_structure() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("659a").unwrap();
        assert!(boolean_function.is_linear_structure(1));
        assert!(!boolean_function.is_linear_structure(7));
        assert!(boolean_function.is_linear_structure(9));
    }

    #[test]
    fn test_walsh_hadamard_values() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("dd0e").unwrap();
        assert_eq!(
            boolean_function.walsh_hadamard_values(),
            [-2, -2, 6, -2, -6, 2, 2, 2, 6, 6, -2, 6, -6, 2, 2, 2]
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("0000").unwrap();
        assert_eq!(
            boolean_function.walsh_hadamard_values(),
            [16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(
            boolean_function.walsh_hadamard_values(),
            [-16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );
    }

    #[test]
    fn test_boolean_function_from_reverse_walsh_transform() {
        let boolean_function = BooleanFunction::from_reverse_walsh_hadamard_transform(&[
            -2, -2, 6, -2, -6, 2, 2, 2, 6, 6, -2, 6, -6, 2, 2, 2,
        ])
        .unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "dd0e");
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let boolean_function = BooleanFunction::from_reverse_walsh_hadamard_transform(&[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
        ])
        .unwrap();
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "ffffffff00000000"
        );
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let boolean_function = BooleanFunction::from_reverse_walsh_hadamard_transform(&[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ])
        .unwrap();
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "ffffffffffffffff0000000000000000"
        );
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let boolean_function = BooleanFunction::from_reverse_walsh_hadamard_transform(&[
            16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ])
        .unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "0000");

        let boolean_function = BooleanFunction::from_reverse_walsh_hadamard_transform(&[
            128, 0, 8, 8, 0, 0, 8, 8, -8, 8, -16, 0, 8, -8, -64, -16, 0, -16, 8, -8, 0, -16, 8, -8,
            8, 8, 0, 0, -8, -8, -16, -16, -8, 8, 0, -16, -8, 8, 0, 16, 0, 0, -8, 24, 16, -48, 8, 8,
            8, 8, 16, 16, 8, 8, -16, 16, 0, 16, -8, 8, -16, 32, 8, 24, 8, 8, 0, 0, -8, -8, 16, -16,
            0, 16, 8, -8, 0, -16, 8, -8, -8, 8, -16, 0, 8, -8, 0, 16, -32, 0, 8, 8, 0, 0, 8, 8, 0,
            0, -8, -8, -16, -16, 8, 8, 8, -8, -16, 0, 8, -8, -16, 0, 0, -16, -8, 8, -16, 0, 8, -8,
            -8, -8, 0, 0, -8, -8, 0, 0, 12, 4, 4, -4, -4, -12, 20, -20, 4, 12, 12, -12, 4, -20, 12,
            -12, -4, 4, -12, -4, 12, -12, 4, 12, -28, -4, 12, 4, 4, -4, 12, 4, 4, -4, -4, -12, -12,
            -20, 12, 4, 12, -12, -12, -4, 12, -12, -12, -4, 4, -20, -4, 4, -12, -4, 12, -12, -4,
            -12, 4, -4, -4, -12, 4, -4, -4, 4, 4, 12, -4, 4, 4, 12, -12, 12, -20, 4, 4, -4, 60,
            -12, -4, -12, 4, -4, -4, -12, 4, -4, 4, 12, -4, 4, -12, -4, -20, -12, -12, -20, -4, 84,
            -12, -20, -4, -12, -4, -28, -12, -4, 12, 52, 4, -20, 4, -20, 12, -12, 4, -20, -20, -12,
            -4, -12, -12, -20, -20, 4, 4, -4,
        ])
        .unwrap();
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "80921c010276c44224422441188118822442244118811880400810a80e200425"
        );
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let boolean_function = BooleanFunction::from_reverse_walsh_hadamard_transform(&[
            128, 0, 8, 8, 0, 0, 8, 8, -8, 8, -16, 0, 8, -8, -64, -16, 0, -16, 8, -8, 0, -16, 8, -8,
            8, 8, 0, 0, -8, -8, -16, -16, -8, 8, 0, -16, -8, 8, 0, 16, 0, 0, -8, 24, 16, -48, 8, 8,
            8, 8, 16, 16, 8, 8, -16, 16, 0, 16, -8, 8, -16, 32, 8, 24, 8, 8, 0, 0, -8, -8, 16, -16,
            0, 16, 8, -8, 0, -16, 8, -8, -8, 8, -16, 0, 8, -8, 0, 16, -32, 0, 8, 8, 0, 0, 8, 8, 0,
            0, -8, -8, -16, -16, 8, 8, 8, -8, -16, 0, 8, -8, -16, 0, 0, -16, -8, 8, -16, 0, 8, -8,
            -8, -8, 0, 0, -8, -8, 0, 0, 12, 4, 4, -4, -4, -12, 20, -20, 4, 12, 12, -12, 4, -20, 12,
            -12, -4, 4, -12, -4, 12, -12, 4, 12, -28, -4, 12, 4, 4, -4, 12, 4, 4, -4, -4, -12, -12,
            -20, 12, 4, 12, -12, -12, -4, 12, -12, -12, -4, 4, -20, -4, 4, -12, -4, 12, -12, -4,
            -12, 4, -4, -4, -12, 4, -4, -4, 4, 4, 12, -4, 4, 4, 12, -12, 12, -20, 4, 4, -4, 60,
            -12, -4, -12, 4, -4, -4, -12, 4, -4, 4, 12, -4, 4, -12, -4, -20, -12, -12, -20, -4, 84,
            -12, -20, -4, -12, -4, -28, -12, -4, 12, 52, 4, -20, 4, -20, 12, -12, 4, -20, -20, -12,
            -4, -12, -12, -20, -20, 4, 4,
        ]);
        assert!(boolean_function.is_err());
        assert_eq!(boolean_function.unwrap_err(), InvalidWalshValuesCount(255));

        let boolean_function =
            BooleanFunction::from_reverse_walsh_hadamard_transform(&[128]);
        assert!(boolean_function.is_err());
        assert_eq!(boolean_function.unwrap_err(), InvalidWalshValuesCount(1));
    }

    #[test]
    fn test_correlation_immunity() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("dd0e").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 0);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("55C3AAC3").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 1);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1f").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 4);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000000000000000000000000000")
                .unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 7);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.correlation_immunity(), 2);
    }

    #[test]
    fn test_resiliency_order() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("dd0e").unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("55C3AAC3").unwrap();
        assert_eq!(boolean_function.resiliency_order(), Some(1));

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1f").unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ffff").unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("00000000000000000000000000000000")
                .unwrap();
        assert_eq!(boolean_function.resiliency_order(), None);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.resiliency_order(), Some(2));
    }

    #[test]
    fn test_biguint_truth_table() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(
            boolean_function.biguint_truth_table().to_str_radix(16),
            "7969817cc5893ba6ac326e47619f5ad0"
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(
            boolean_function.biguint_truth_table().to_str_radix(16),
            "1e"
        );
    }

    #[test]
    fn test_try_u64_truth_table() {
        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        assert_eq!(boolean_function.try_u64_truth_table(), None);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.try_u64_truth_table(), Some(30));
    }

    #[test]
    fn test_xor() {
        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("80921c010276c440400810a80e200425")
                .unwrap();
        let boolean_function2 =
            BooleanFunction::from_hex_string_truth_table("22442244118811882244224411881188")
                .unwrap();
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "a2d63e4513fed5c8624c32ec1fa815ad"
        );
        assert_eq!(boolean_function.variables_count(), 7);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("ab").unwrap();
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function3 = boolean_function.clone() ^ boolean_function2.clone();
        boolean_function ^= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
    }

    #[test]
    fn test_add() {
        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("80921c010276c440400810a80e200425")
                .unwrap();
        let boolean_function2 =
            BooleanFunction::from_hex_string_truth_table("22442244118811882244224411881188")
                .unwrap();
        let boolean_function3 = boolean_function.clone() + boolean_function2.clone();
        boolean_function += boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "a2d63e4513fed5c8624c32ec1fa815ad"
        );
        assert_eq!(boolean_function.variables_count(), 7);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("ab").unwrap();
        let boolean_function3 = boolean_function.clone() + boolean_function2.clone();
        boolean_function += boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function3 = boolean_function.clone() + boolean_function2.clone();
        boolean_function += boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function3 = boolean_function.clone() + boolean_function2.clone();
        boolean_function += boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "b5");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
    }

    #[test]
    fn test_and() {
        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("4f1ead396f247a0410bdb210c006eab568ab4bfa8acb7a13b14ede67096c6eed")
                .unwrap();
        let boolean_function2 =
            BooleanFunction::from_hex_string_truth_table("c870974094ead8a96a450b2ef33486b4e61a4c5e97816f7a7bae007d4c53fc7d")
                .unwrap();
        let boolean_function3 = boolean_function.clone() & boolean_function2.clone();
        boolean_function &= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "481085000420580000050200c00482b4600a485a82816a12310e006508406c6d"
        );
        assert_eq!(boolean_function.variables_count(), 8);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("ab").unwrap();
        let boolean_function3 = boolean_function.clone() & boolean_function2.clone();
        boolean_function &= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "0a");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function3 = boolean_function.clone() & boolean_function2.clone();
        boolean_function &= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "0a");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function3 = boolean_function.clone() & boolean_function2.clone();
        boolean_function &= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "0a");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
    }

    #[test]
    fn test_mul() {
        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("4f1ead396f247a0410bdb210c006eab568ab4bfa8acb7a13b14ede67096c6eed")
                .unwrap();
        let boolean_function2 =
            BooleanFunction::from_hex_string_truth_table("c870974094ead8a96a450b2ef33486b4e61a4c5e97816f7a7bae007d4c53fc7d")
                .unwrap();
        let boolean_function3 = boolean_function.clone() & boolean_function2.clone();
        boolean_function &= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "481085000420580000050200c00482b4600a485a82816a12310e006508406c6d"
        );
        assert_eq!(boolean_function.variables_count(), 8);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("ab").unwrap();
        let boolean_function3 = boolean_function.clone() * boolean_function2.clone();
        boolean_function *= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "0a");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function =
            BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function2: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function3 = boolean_function.clone() * boolean_function2.clone();
        boolean_function *= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "0a");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let mut boolean_function: BooleanFunction = super::BigBooleanFunction::from_truth_table(
            BigUint::from_str_radix("ab", 16).unwrap(),
            3,
        ).into();
        let boolean_function2 = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let boolean_function3 = boolean_function.clone() * boolean_function2.clone();
        boolean_function *= boolean_function2.clone();
        assert_eq!(&boolean_function, &boolean_function3);
        assert_eq!(boolean_function.printable_hex_truth_table(), "0a");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Small
        );
        assert_eq!(
            boolean_function3.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
    }

    #[test]
    fn test_walsh_fourier_values() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("ff").unwrap();
        assert_eq!(
            boolean_function.walsh_fourier_values(),
            [8, 0, 0, 0, 0, 0, 0, 0]
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("00").unwrap();
        assert_eq!(
            boolean_function.walsh_fourier_values(),
            [0, 0, 0, 0, 0, 0, 0, 0]
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("0f").unwrap();
        assert_eq!(
            boolean_function.walsh_fourier_values(),
            [4, 0, 0, 0, 4, 0, 0, 0]
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("55").unwrap();
        assert_eq!(
            boolean_function.walsh_fourier_values(),
            [4, 4, 0, 0, 0, 0, 0, 0]
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("aa").unwrap();
        assert_eq!(
            boolean_function.walsh_fourier_values(),
            [4, -4, 0, 0, 0, 0, 0, 0]
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("8001").unwrap();
        assert_eq!(
            boolean_function.walsh_fourier_values(),
            [2, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 0, 2]
        );
    }

    #[test]
    fn test_boolean_function_from_reverse_walsh_fourier_transform() {
        let boolean_function =
            BooleanFunction::from_reverse_walsh_fourier_transform(&[8, 0, 0, 0, 0, 0, 0, 0])
                .unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "ff");

        let boolean_function =
            BooleanFunction::from_reverse_walsh_fourier_transform(&[0, 0, 0, 0, 0, 0, 0, 0])
                .unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "00");

        let boolean_function = BooleanFunction::from_reverse_walsh_fourier_transform(&[
            2, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 0, 2,
        ])
        .unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "8001");
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table(
            "80921c010276c44224422441188118822442244118811880400810a80e200425",
        )
        .unwrap();
        let walsh_fourier_values = boolean_function.walsh_fourier_values();
        let boolean_function2 =
            BooleanFunction::from_reverse_walsh_fourier_transform(&walsh_fourier_values)
                .unwrap();
        assert_eq!(
            boolean_function2.printable_hex_truth_table(),
            "80921c010276c44224422441188118822442244118811880400810a80e200425"
        );

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff").unwrap();
        let walsh_fourier_values = boolean_function.walsh_fourier_values();
        let boolean_function2 =
            BooleanFunction::from_reverse_walsh_fourier_transform(&walsh_fourier_values)
                .unwrap();
        assert_eq!(boolean_function2.printable_hex_truth_table(), "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff");
        assert_eq!(
            boolean_function2.get_boolean_function_type(),
            BooleanFunctionType::Big
        );
        assert_eq!(boolean_function2.variables_count(), 9);
    }

    #[test]
    fn test_support() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("00").unwrap();
        assert_eq!(boolean_function.support().len(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ff").unwrap();
        assert_eq!(boolean_function.support().len(), 8);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        assert_eq!(boolean_function.support().len(), 4);
    }

    #[test]
    fn test_propagation_criterion() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("00").unwrap();
        assert_eq!(boolean_function.propagation_criterion(), 0);

        let boolean_function = BooleanFunction::from_hex_string_truth_table("ff").unwrap();
        assert_eq!(boolean_function.propagation_criterion(), 0);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("288d1b41").unwrap();
        assert_eq!(boolean_function.propagation_criterion(), 3);
    }

    #[test]
    fn test_iter() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let mut iter = boolean_function.iter();
        assert_eq!(iter.next().unwrap(), false);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), false);
        assert_eq!(iter.next().unwrap(), false);
        assert_eq!(iter.next().unwrap(), false);
        assert!(iter.next().is_none());

        let iter = boolean_function.iter();
        assert_eq!(iter.count(), 8);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let mut iter = boolean_function.iter();
        assert_eq!(iter.next(), Some(false));
        assert_eq!(iter.next(), Some(false));
        assert_eq!(iter.next(), Some(false));
        assert_eq!(iter.next(), Some(false));
        assert_eq!(iter.next(), Some(true));

        let iter = boolean_function.iter();
        assert_eq!(iter.count(), 128);
    }

    #[test]
    fn test_boolean_function_from_u64_truth_table() {
        let boolean_function = BooleanFunction::from_u64_truth_table(30, 3).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "1e");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let boolean_function = BooleanFunction::from_u64_truth_table(30, 7);
        assert!(boolean_function.is_err());

        let boolean_function = BooleanFunction::from_u64_truth_table(300, 3);
        assert!(boolean_function.is_err());
    }

    #[test]
    fn test_boolean_function_from_biguint_truth_table() {
        let boolean_function =
            BooleanFunction::from_biguint_truth_table(&BigUint::from(30u32), 3).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "1e");
        assert_eq!(boolean_function.variables_count(), 3);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Small
        );

        let boolean_function =
            BooleanFunction::from_biguint_truth_table(&BigUint::from(30u32), 7).unwrap();
        assert_eq!(
            boolean_function.printable_hex_truth_table(),
            "0000000000000000000000000000001e"
        );
        assert_eq!(boolean_function.variables_count(), 7);
        assert_eq!(
            boolean_function.get_boolean_function_type(),
            BooleanFunctionType::Big
        );

        let boolean_function = BooleanFunction::from_u64_truth_table(300, 3);
        assert!(boolean_function.is_err());

        let boolean_function = BooleanFunction::from_u64_truth_table(300, 32);
        assert!(boolean_function.is_err());
    }

    #[test]
    fn test_get_1_local_neighbor() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let neighbor = boolean_function.get_1_local_neighbor(1);
        assert_eq!(neighbor.printable_hex_truth_table(), "1c");

        let boolean_function = BooleanFunction::from_hex_string_truth_table("80921c010276c440400810a80e200425").unwrap();
        let neighbor = boolean_function.get_1_local_neighbor(0);
        assert_eq!(neighbor.printable_hex_truth_table(), "80921c010276c440400810a80e200424");
    }

    #[test]
    fn test_close_balanced_functions_iterator() {
        let balanced_function = BooleanFunction::from_hex_string_truth_table("aaaa").unwrap();
        assert!(balanced_function.close_balanced_functions_iterator().is_err());

        let bent_function = BooleanFunction::from_hex_string_truth_table("ac90").unwrap();
        let close_balanced_iterator = bent_function.close_balanced_functions_iterator();
        assert!(close_balanced_iterator.is_ok());
        let close_balanced_iterator = close_balanced_iterator.unwrap();
        assert_eq!(close_balanced_iterator.into_iter().count(), 45); // 10 choose 2

        let mut close_balanced_iterator = bent_function.close_balanced_functions_iterator().unwrap();
        assert!(close_balanced_iterator.all(|f| f.is_balanced()));


        let balanced_function = BooleanFunction::from_hex_string_truth_table(
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        ).unwrap();
        assert!(balanced_function.close_balanced_functions_iterator().is_err());

        let bent_function = BooleanFunction::from_hex_string_truth_table(
            "80329780469d0b85cd2ad63e1a6ba42adbd83c9a0c55e4e8c99f227b0ffc1418"
        ).unwrap();
        let close_balanced_iterator = bent_function.close_balanced_functions_iterator();
        assert!(close_balanced_iterator.is_ok());
        let mut close_balanced_iterator = close_balanced_iterator.unwrap();
        for _ in 0..10 {
            assert!(close_balanced_iterator.next().unwrap().is_balanced());
        }
    }

    #[test]
    fn test_from_anf_polynomial_str() {
        let rule_30_anf_str = "x0*x1 + x0 + x1 + x2";
        let rule_30_function = BooleanFunction::from_anf_polynomial_str(rule_30_anf_str, 3).unwrap();
        assert_eq!(rule_30_function.printable_hex_truth_table(), "1e");

        let anf_str = "x0*x1*x2*x3*x4*x5*x6 + x7";
        let boolean_function = BooleanFunction::from_anf_polynomial_str(anf_str, 8).unwrap();
        assert_eq!(boolean_function.printable_hex_truth_table(), "7fffffffffffffffffffffffffffffff80000000000000000000000000000000");
    }

    #[test]
    fn test_from_anf_polynomial() {
        let rule_30_anf = AnfPolynomial::from_str("x0*x1 + x0 + x1 + x2", 3).unwrap();
        let rule_30_function = BooleanFunction::from_anf_polynomial(&rule_30_anf);
        assert_eq!(rule_30_function.printable_hex_truth_table(), "1e");

        let anf = AnfPolynomial::from_str("x0*x1*x2*x3*x4*x5*x6 + x7", 8).unwrap();
        let boolean_function = BooleanFunction::from_anf_polynomial(&anf);
        assert_eq!(boolean_function.printable_hex_truth_table(), "7fffffffffffffffffffffffffffffff80000000000000000000000000000000");
    }
}
