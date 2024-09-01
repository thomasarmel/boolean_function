use thiserror::Error;

#[derive(Error, Debug)]
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
}
