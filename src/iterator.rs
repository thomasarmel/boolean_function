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
        if self.current_index > self.max_index { // limit of 31 variables
            return None;
        }
        let result = self.inner_bool_func.compute_cellular_automata_rule(self.current_index);
        self.current_index += 1;
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use crate::boolean_function_from_hex_string_truth_table;
    use crate::iterator::BooleanFunctionIterator;

    #[test]
    fn test_boolean_function_iterator() {
        let boolean_function = boolean_function_from_hex_string_truth_table("1e").unwrap();
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

        let boolean_function = boolean_function_from_hex_string_truth_table("7969817CC5893BA6AC326E47619F5AD0").unwrap();
        let mut iterator = BooleanFunctionIterator::new(boolean_function);
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(false));
        assert_eq!(iterator.next(), Some(true));
    }
}