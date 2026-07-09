# POLSYS-RS

This crate provides a Rust wrapper around the `POLSYS_PLP` Fortran library
by Wise, Sommese, and Watson in their [ACM paper](https://doi.org/10.1145/347837.347885).
The library can solve a system of N polynomial equations with complex
coefficients with N unknowns using a globally convergence homotopy method.

This is currently **very experimental** since doing Fortran bindings is hard!

## Example

The following snippet solves a system of 4 polynomial equations in 4 complex
unknowns:

```math
\begin{cases}
x_1^2 + x_2 - 1 = 0, \\
x_2^2 + x_3 - 1 = 0, \\
x_3^2 + x_4 - 1 = 0, \\
x_4^2 + x_1 - 1 = 0,
\end{cases}
```

using the homotopy method and prints all the computed roots.

```rust
use num::complex::Complex64;

use polsys::{Polynomial, SolveState, make_homogeneous_partition};

fn main() {
    // Each equation is given as a list of (degree_tuple, coefficient) pairs,
    // where the degree tuple holds the exponents of (x1, x2, x3, x4).
    let mut poly = Polynomial::<4>::new(vec![
        vec![
            ([2, 0, 0, 0], Complex64::new(1.0, 0.0)),
            ([0, 1, 0, 0], Complex64::new(1.0, 0.0)),
            ([0, 0, 0, 0], Complex64::new(-1.0, 0.0)),
        ],
        vec![
            ([0, 2, 0, 0], Complex64::new(1.0, 0.0)),
            ([0, 0, 1, 0], Complex64::new(1.0, 0.0)),
            ([0, 0, 0, 0], Complex64::new(-1.0, 0.0)),
        ],
        vec![
            ([0, 0, 2, 0], Complex64::new(1.0, 0.0)),
            ([0, 0, 0, 1], Complex64::new(1.0, 0.0)),
            ([0, 0, 0, 0], Complex64::new(-1.0, 0.0)),
        ],
        vec![
            ([0, 0, 0, 2], Complex64::new(1.0, 0.0)),
            ([1, 0, 0, 0], Complex64::new(1.0, 0.0)),
            ([0, 0, 0, 0], Complex64::new(-1.0, 0.0)),
        ],
    ])
    .unwrap();

    // Use a standard (1-homogeneous) partition of the variables.
    let mut part = make_homogeneous_partition(poly.len()).unwrap();

    // Solve with a globally convergent homotopy.
    let result = SolveState::solve(
        &mut poly, &mut part,
        1.0e-8, 1.0e-8, 1.0e-8, 1, false, false,
    )
    .unwrap();

    println!("Found {} roots:", result.roots.len());
    for (i, root) in result.roots.iter().enumerate() {
        println!("{i:2}: {:?}", &root[..poly.len()]);
    }
}
```

A runnable version of this example lives in `examples/solve_4d.rs`.

## Development

### Links

- [Code](https://github.com/alexfikl/polsys-rs).

Other known libraries that may do a better job if they had Rust bindings:

- [pypolsys](https://github.com/nennigb/pypolsys) (Python). This has been a
  major inspiration for this wrapper, since we have similar goals.
- [PHCpack](https://github.com/janverschelde/PHCpack) (Ada, C, Julia, Maple, Python).
- [Bertini 2.0](https://github.com/bertiniteam/b2) (C++, Python).
- [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) (Julia).

## License

The Rust and Fortran code in this library is MIT licensed. However, the Fortran code
for the `POLSYS_PLP` library is under the [ACM Software License Agreement](https://www.acm.org/publications/policies/software-copyright-notice).
If you want to use this library for anything serious, make sure to enquire with the
original authors.
