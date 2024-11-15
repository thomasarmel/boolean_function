use std::fs::File;
use std::io::Read;
use num_bigint::BigUint;
use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use boolean_function::{BigBooleanFunction, BooleanFunctionImpl};

// Thanks https://langevin.univ-tln.fr/project/genbent/genbent.html
const FILENAME: &'static str = "data/bent-1.zip";

fn main() {
    let mut bent_file = File::open(FILENAME).unwrap();
    let mut bent_function_buffer = [0u8; 32];

    loop {
        if bent_file.read(&mut bent_function_buffer).unwrap() == 0 {
            println!("EOF");
            break;
        }
        let bent_function_tt = BigUint::from_bytes_le(&bent_function_buffer);
        let bent_function = BigBooleanFunction::from_truth_table(bent_function_tt.clone(), 8);

        let close_balanced_iterator = bent_function.close_balanced_functions_iterator().unwrap();
        close_balanced_iterator.par_bridge()
            .filter(|balanced_function| balanced_function.nonlinearity() > 116)
            .for_each(|balanced_function| {
                println!("{}: {}", balanced_function.printable_hex_truth_table(), balanced_function.nonlinearity());
            });
    }
}
