# Boolean function analysis library

*Mathematical analysis of Boolean functions.*

---

## Introduction

Boolean functions play a dominant role in computer science and cryptography.
They are used in many applications, such as block ciphers, hash functions, error-correcting codes, cellular automata, boolean networks, and more.

These functions take as input a vector of bits and return a single bit. An $n$-variable Boolean function $f$ is noted $f: \mathbb{F}_2^n \to \mathbb{F}_2$.

This library provides tools to analyze and manipulate boolean functions.
It is inspired by the excellent [Sage Math implementation](https://doc.sagemath.org/html/en/reference/cryptography/sage/crypto/boolean_function.html).
This Rust implementation aims to be faster, and easier to parallelize, which could be interesting for some exhaustive searches.

## Usage

In order to use this library, add the crate as a dependency in your `Cargo.toml` file:

```bash
cargo add boolean_function
```

Here is a simple example to illustrate how to use this library:

```rust
use boolean_function::*;
use num_bigint::BigUint;
use num_traits::Num;

fn main() {
    // Create a Boolean function from its string truth table
    const TRUTH_TABLE_STR: &'static str = "0113077C165E76A8";
    let f: BooleanFunction = BooleanFunction::from_hex_string_truth_table(TRUTH_TABLE_STR)
        .unwrap();
    
    // How many variables does the function have?
    assert_eq!(f.variables_count(), 6);

    // Check if the function is bent
    assert!(f.is_bent());
    
    // If you already know that the function has a small number of variables,
    // you can use SmallBooleanFunction
    // It will be faster runtime and use less memory
    let f_small: SmallBooleanFunction = SmallBooleanFunction::from_truth_table(0xac90, 4)
        .unwrap();
    assert!(f_small.is_bent());
    assert!(!f_small.is_balanced()); // A bent function is not balanced
    // ANF of the function, `*` is AND and `+` is XOR
    assert_eq!(
        f_small.algebraic_normal_form().to_string(),
        "x0*x2 + x1*x2 + x1*x3 + x2*x3 + x2");
    
    // If your function has more than 6 variables, you can use BigBooleanFunction
    // So that you will avoid runtime calls to the V-table
    const AES_COMP_TT: &'static str = "4f1ead396f247a0410bdb210c006eab568ab4bfa8acb7a13b14ede67096c6eed";
    let f_big: BigBooleanFunction = BigBooleanFunction::from_truth_table(
        BigUint::from_str_radix(AES_COMP_TT, 16).unwrap(),
        8);
    assert_eq!(f_big.nonlinearity(), 112); // AES S-Box has a nonlinearity of 112
}
```

As you can see, this library provides 3 types of Boolean functions:
- `BooleanFunction`: for Boolean functions with an arbitrary number of variables (up to 31). This is a wrapper for `Box<dyn BooleanFunctionImpl>`.
- `SmallBooleanFunction`: for Boolean functions with up to 6 variables. Internally, the truth table is stored as a `u64`.
- `BigBooleanFunction`: for Boolean functions with more than 6 variables. Internally, the truth table is stored as a `BigUint`.

All these types have the methods of the `BooleanFunctionImpl` trait, which provides methods to analyze and manipulate Boolean functions.

See [examples](examples/) directory for more examples.

### Limitations

- Some constructor methods, like `BooleanFunction::from_hex_string_truth_table`, don't work for Boolean functions with 0 or 1 input variables (because there cannot be expressed as hex string). As these functions are trivial, it should not be a problem.
- This crate doesn't work for Boolean functions with more than 31 input variables. This is because internally some calculations are done with `u32`, which would overflow for more than 31 variables. But your computer probably doesn't have enough memory to store the truth table anyway, neither the calculation power to iterate over any subset of the $3.103 \cdot 10^{1292913986}$ 32-variable Boolean functions.

## Performance

This library aims to be as fast as possible. Here are some tricks you could use to improve the performance of your program:

### Parallelize your code

This crate is designed to make Boolean function exhaustive search easy to parallelize.

#### Example: Check that there are $\binom{16}{8} = 12870$ balanced 4-variable Boolean functions

```rust
use rayon::prelude::*;
use boolean_function::*;

fn main() {
    // Parallel exhaustive search on all 4-variable Boolean functions
    let count = (0..(1 << 16)).into_par_iter()
        .filter(|truth_table| {
            let f: BooleanFunction = BooleanFunction::from_u64_truth_table(*truth_table, 4)
                .unwrap();
            f.is_balanced()
        }).count();
    assert_eq!(count, 12870);
}
```

### Use `SmallBooleanFunction` when possible

This crate can handle Boolean functions with arbitrary numbers of variables (up to 31).
For an $n$-variable Boolean function, the truth table has $2^n$ Boolean entries.

The crate is optimized to use a `u64` to store the truth table when the number of variables is less than or equal to 6, and `BigUint` otherwise (as today most computers have 64-bit CPUs).

When using `BooleanFunction` type, the crate will detect runtime whether to use `SmallBooleanFunction` (truth table stored as `u64`) or `BigBooleanFunction` (truth table stored as `BigUint`) depending on the number of variables.

If you already know that your function will have less than 6 variables, you can use the `SmallBooleanFunction` type directly.
So that you can avoid using polymorphism, and the runtime check to the V-table.

### Disable safety checks

This crate is designed to be safe, meaning no undefined behavior should occur.
So there are some safety checks that are performed at runtime, for example to ensure that the truth table has the correct size, or the number of variables is not too large.

However, once you are sure that your code is correct, you can disable the safety checks:
simply add the `unsafe_disable_safety_checks` feature to the `boolean_function` dependency in your `Cargo.toml` file.

**NB:** Activating this feature could lead to undefined behaviour if the code is incorrect, so please test it with the safety checks enabled first.

### Native CPU target

If you already know that your program will run on a CPU that supports the same instruction set extension as the one you are compiling for, you can specify native compilation.

In `.cargo/config.toml`:

```toml
[build]
rustflags = ["-C", "target-cpu=native"]
```

Here is a simple example to illustrate how it changes the produced assembly.
```rust
pub fn boolean_function_balanced(truth_table: u64, variables_count: usize) -> bool {
    let expected_set_number: u32 = 1 << (variables_count - 1);
    truth_table.count_ones() == expected_set_number
}
```
Let's compile using rustc 1.81.0 for x86_64 target and O3 optimization level.

Without native target:

```assembly
boolean_function_balanced:
        lea     ecx, [rsi - 1]
        mov     eax, 1
        shl     eax, cl
        mov     rcx, rdi
        shr     rcx
        movabs  rdx, 6148914691236517205
        and     rdx, rcx
        sub     rdi, rdx
        movabs  rcx, 3689348814741910323
        mov     rdx, rdi
        and     rdx, rcx
        shr     rdi, 2
        and     rdi, rcx
        add     rdi, rdx
        mov     rcx, rdi
        shr     rcx, 4
        add     rcx, rdi
        movabs  rdx, 1085102592571150095
        and     rdx, rcx
        movabs  rcx, 72340172838076673
        imul    rcx, rdx
        shr     rcx, 56
        cmp     eax, ecx
        sete    al
        ret
```

With native target:

```assembly
boolean_function_balanced:
        dec     sil
        mov     eax, 1
        shlx    eax, eax, esi
        popcnt  rcx, rdi
        cmp     eax, ecx
        sete    al
        ret
```

## Contributing

This library is open-source, and contributions are welcome. Please do not hesitate to open a pull request or an issue!

---

# Some theory

## Algebraic Normal Form

ANF transform is a method to convert a boolean function **from its truth table** representation **to its Algebraic Normal Form** (ANF) representation.
ANF is a representation of a boolean function as **a XOR sum of AND monomials**.

For example, let's consider the following 2-variables boolean function truth table:

| x1 | x0 | f(x1, x2) |
|----|----|-----------|
| 0  | 0  | 1         |
| 0  | 1  | 1         |
| 1  | 0  | 0         |
| 1  | 1  | 0         |

The ANF representation of this function is: `f(x1, x0) = 1 XOR x1`.

### Boolean function degree

The degree of a boolean function is the maximum degree of its monomials in its ANF representation.

For example, the degree of the function `f(x1, x0) = x0 XOR x0.x1` is 2.

### Unsigned integer representation

Truth tables are often represented as unsigned integers.
For example, the truth table of the following function:

| x1 | x0 | f(x1, x2) |
|----|----|-----------|
| 0  | 0  | 1         |
| 0  | 1  | 1         |
| 1  | 0  | 0         |
| 1  | 1  | 0         |

has the unsigned integer representation: `0b0011 = 3`.

In the same way, the ANF representation of boolean functions can be represented as unsigned integers.
For example, the ANF representation of the previous function `f(x1, x0) = 1 XOR x1`.
This can be written as `0b0101 = 5`.
To find the function from the ANF representation in binary, separate each bit and indicate its position starting **from the least significant bit**.

Position **0 corresponds to constant 1**, for the other positions convert them to binary.
This binary number gives the monomial, starting with x0 for the least significant bit.
All you have to do is multiply this monomial by the corresponding bit of the ANF representation.

|                               |               |            |            |           |
|-------------------------------|---------------|------------|------------|-----------|
| **ANF binary representation** | 0             | 1          | 0          | 1         |
| **bit position**              | 3             | 2          | 1          | 0         |
| **bit position in binary**    | 0b11          | 0b10       | 0b01       | /         |
| **monomial**                  | x1.x0         | x1         | x0         | 1         |
| **function**                  | **0**.(x1.x0) | **1**.(x1) | **0**.(x0) | **1**.(1) |

The final function is: `f(x1, x0) = x1 XOR 1`. Its degree is 1.