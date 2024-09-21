use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum BooleanFunctionError {
    #[error("Hex truth table length must be a power of 2")]
    WrongStringHexTruthTableLength,
    #[error("Error parsing string hex number")]
    StringHexParseError,
    #[error("Too big variable count, must be <= {0}")]
    TooBigVariableCount(usize),
    #[error("Unexpected error, this shouldn't happen")]
    UnexpectedError,
    #[error("Too big derivative direction, must be <= {0}")]
    TooBigDerivativeDirection(u32),
    #[error("Invalid number of Walsh values {0}, should be 2^n, n >= 2")]
    InvalidWalshValuesCount(usize),
    #[error("Truth table is too big for variables count")]
    TooBigTruthTableForVarCount,
}

pub(crate) const XOR_DIFFERENT_VAR_COUNT_PANIC_MSG: &'static str = "XOR operation requires the same number of variables in both functions";
#[cfg(not(feature = "unsafe_disable_safety_checks"))]
pub(crate) const TRUTH_TABLE_TOO_BIG_VAR_COUNT_PANIC_MSG: &'static str = "Truth table is too big for variables count";
#[cfg(not(feature = "unsafe_disable_safety_checks"))]
pub(crate) const POLYNOMIAL_ANF_TOO_BIG_VAR_COUNT_PANIC_MSG: &'static str = "Polynomial ANF is too big for variables count";