//! Iterators for Boolean functions.

use std::slice::Iter;
use itertools::{Combinations, Itertools};
use num_bigint::BigUint;
use ouroboros::self_referencing;
use crate::{BigBooleanFunction, BooleanFunction, BooleanFunctionImpl, SmallBooleanFunction};

/// Iterator for the successive values of a Boolean function.
///
/// Example:
/// ```rust
/// use boolean_function::BooleanFunctionImpl;
/// use boolean_function::BooleanFunction;
/// use boolean_function::BooleanFunctionIterator;
///
/// let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
/// let mut iterator = boolean_function.iter();
/// assert_eq!(iterator.next(), Some(false));
/// assert_eq!(iterator.next(), Some(true));
/// assert_eq!(iterator.next(), Some(true));
/// assert_eq!(iterator.next(), Some(true));
/// assert_eq!(iterator.next(), Some(true));
/// assert_eq!(iterator.next(), Some(false));
/// assert_eq!(iterator.next(), Some(false));
/// assert_eq!(iterator.next(), Some(false));
/// assert_eq!(iterator.next(), None);
/// ```
pub struct BooleanFunctionIterator {
    current_index: u32,
    max_index: u32,
    inner_bool_func: crate::BooleanFunction,
}

impl BooleanFunctionIterator {
    pub(crate) fn new(bool_func: crate::BooleanFunction) -> Self {
        BooleanFunctionIterator {
            current_index: 0,
            max_index: bool_func.get_max_input_value(),
            inner_bool_func: bool_func,
        }
    }
}

impl Iterator for BooleanFunctionIterator {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_index > self.max_index {
            // limit of 31 variables
            return None;
        }
        let result = self
            .inner_bool_func
            .compute_cellular_automata_rule(self.current_index);
        self.current_index += 1;
        Some(result)
    }
}

#[self_referencing]
pub struct SmallCloseBalancedFunctionIterator {
    original_truth_table: u64,
    var_count: usize,
    flippable_bit_positions: Vec<usize>,
    flips_count: usize,
    #[borrows(flippable_bit_positions)]
    #[not_covariant]
    combination_generator: Combinations<Iter<'this, usize>>,
}

impl SmallCloseBalancedFunctionIterator {
    pub(crate) fn create(original_function: &SmallBooleanFunction, flippable_bit_positions: Vec<usize>, flips_count: usize) -> Self {
        SmallCloseBalancedFunctionIteratorBuilder {
            original_truth_table: original_function.get_truth_table_u64(),
            var_count: original_function.variables_count(),
            flippable_bit_positions,
            flips_count,
            combination_generator_builder: |flippable_pos| flippable_pos.iter().combinations(flips_count),
        }.build()
    }
}

impl Iterator for SmallCloseBalancedFunctionIterator {
    type Item = SmallBooleanFunction;

    fn next(&mut self) -> Option<Self::Item> {
        let mut new_tt = self.borrow_original_truth_table().clone();
        let var_count = self.borrow_var_count().clone();

        let bits_to_flip = self.with_mut(|fields| {
            fields.combination_generator.next()
        })?;
        for bit_flip_pos in bits_to_flip {
            new_tt ^= 1u64 << bit_flip_pos;
        }

        Some(SmallBooleanFunction::from_truth_table_unchecked(new_tt, var_count))
    }
}

#[self_referencing]
pub struct BigCloseBalancedFunctionIterator {
    original_truth_table: BigUint,
    var_count: usize,
    flippable_bit_positions: Vec<usize>,
    flips_count: usize,
    #[borrows(flippable_bit_positions)]
    #[not_covariant]
    combination_generator: Combinations<Iter<'this, usize>>,
}

impl BigCloseBalancedFunctionIterator {
    pub(crate) fn create(original_function: &BigBooleanFunction, flippable_bit_positions: Vec<usize>, flips_count: usize) -> Self {
        BigCloseBalancedFunctionIteratorBuilder {
            original_truth_table: original_function.biguint_truth_table(),
            var_count: original_function.variables_count(),
            flippable_bit_positions,
            flips_count,
            combination_generator_builder: |flippable_pos| flippable_pos.iter().combinations(flips_count),
        }.build()
    }
}

impl Iterator for BigCloseBalancedFunctionIterator {
    type Item = BigBooleanFunction;

    fn next(&mut self) -> Option<Self::Item> {
        let mut new_tt = self.borrow_original_truth_table().clone();
        let var_count = self.borrow_var_count().clone();

        let bits_to_flip = self.with_mut(|fields| {
            fields.combination_generator.next()
        })?;
        for bit_flip_pos in bits_to_flip {
            let bit_flip_pos = *bit_flip_pos as u64;
            new_tt.set_bit(bit_flip_pos, !new_tt.bit(bit_flip_pos))
        }

        Some(BigBooleanFunction::from_truth_table(new_tt, var_count))
    }
}

pub enum CloseBalancedFunctionIterator {
    Small(SmallCloseBalancedFunctionIterator),
    Big(BigCloseBalancedFunctionIterator),
}

impl Iterator for CloseBalancedFunctionIterator {
    type Item = BooleanFunction;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            CloseBalancedFunctionIterator::Small(it) => Some(it.next()?.into()),
            CloseBalancedFunctionIterator::Big(it) => Some(it.next()?.into())
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::BooleanFunction;
    use crate::iterator::BooleanFunctionIterator;
    #[test]
    fn test_boolean_function_iterator() {
        let boolean_function = BooleanFunction::from_hex_string_truth_table("1e").unwrap();
        let mut iterator = BooleanFunctionIterator::new(boolean_function);
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(true));
        assert_eq!(iterator.next(), Some(true));
        assert_eq!(iterator.next(), Some(true));
        assert_eq!(iterator.next(), Some(true));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), None);

        let boolean_function =
            BooleanFunction::from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0")
                .unwrap();
        let mut iterator = BooleanFunctionIterator::new(boolean_function);
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(true));
    }
}
