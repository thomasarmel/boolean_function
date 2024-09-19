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
}

pub(crate) const XOR_DIFFERENT_VAR_COUNT_PANIC_MSG: &'static str = "XOR operation requires the same number of variables in both functions";