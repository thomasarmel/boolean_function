use boolean_function::BooleanFunctionImpl;
use boolean_function::BooleanFunction;
use num_bigint::BigUint;
use num_traits::{Num, Zero};
use rayon::iter::{ParallelBridge, ParallelIterator};

fn main() {
    num_iter::range_inclusive(BigUint::zero(), BigUint::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap())
        .into_iter().par_bridge()
        .map(|i| BooleanFunction::from_biguint_truth_table(&i, 8).unwrap())
        .filter(|f| f.is_balanced())
        .filter(|f| f.nonlinearity() >= 112)
        .for_each(|f| println!("{}", f.printable_hex_truth_table()));
}
