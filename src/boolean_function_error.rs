//! Possible errors and panic messages for the boolean function library.

use thiserror::Error;

/// Possible errors for the boolean function library.
#[derive(Error, Debug, PartialEq)]
pub enum BooleanFunctionError {
    /// Error parsing hex string: the string length must be a power of 2 to represent a truth table.
    #[error("Hex truth table length must be a power of 2")]
    WrongStringHexTruthTableLength,
    /// Error parsing hex string: probably some characters are not valid hex characters.
    #[error("Error parsing string hex number")]
    StringHexParseError,
    /// The Boolean function variable count is too big: it must be $\leq 6$ for [crate::SmallBooleanFunction] and $\leq 31$ for [crate::BooleanFunction] / [crate::BigBooleanFunction].
    /// The last limit is due to the maximum number of variables that can be represented in a 32-bit integer.
    #[error("Too big variable count, must be <= {0}")]
    TooBigVariableCount(usize),
    /// An unexpected error occurred. This should not happen as the case is supposed to be handled. Please report this issue to the crate maintainer.
    #[error("Unexpected error, this shouldn't happen. Please report this issue to the crate maintainer.")]
    UnexpectedError,
    /// The derivative direction is too big: it must be $< 2^n$, where n is the number of variables.
    #[error("Too big derivative direction, must be <= {0}")]
    TooBigDerivativeDirection(u32),
    /// The Walsh values count is invalid: it must be $2^n$, where $n \geq 2$. n is the number of variables of the Boolean function.
    #[error("Invalid number of Walsh values {0}, should be 2^n, n >= 2")]
    InvalidWalshValuesCount(usize),
    /// The truth table is too big for the number of variables: it must be $< 2^{2^n}$, where n is the number of variables.
    #[error("Truth table is too big for variables count")]
    TooBigTruthTableForVarCount,
    /// Cannot generate close balanced Boolean function iterator, as the given function is already balanced
    #[error("This Boolean function is already balanced")]
    AlreadyBalanced,
}

pub(crate) const XOR_DIFFERENT_VAR_COUNT_PANIC_MSG: &'static str =
    "XOR operation requires the same number of variables in both functions";

#[cfg(not(feature = "unsafe_disable_safety_checks"))]
pub(crate) const TRUTH_TABLE_TOO_BIG_VAR_COUNT_PANIC_MSG: &'static str =
    "Truth table is too big for variables count";
#[cfg(not(feature = "unsafe_disable_safety_checks"))]
pub(crate) const POLYNOMIAL_ANF_TOO_BIG_VAR_COUNT_PANIC_MSG: &'static str =
    "Polynomial ANF is too big for variables count";
