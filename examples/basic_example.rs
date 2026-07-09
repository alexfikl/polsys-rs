// SPDX-FileCopyrightText: 2026 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use polsys::{Polynomial, SolveState, make_homogeneous_partition, term};

fn main() {
    // A system of 4 polynomial equations in 4 complex unknowns.
    //
    //     f1(x) = x1^2 + x2      - 1
    //     f2(x) = x2^2 + x3      - 1
    //     f3(x) = x3^2 + x4      - 1
    //     f4(x) = x4^2 + x1      - 1
    //
    // Each equation is given as a list of (degree_tuple, coefficient) pairs,
    // where the degree tuple holds the exponents of (x1, x2, x3, x4).
    let mut poly = Polynomial::<4>::new(vec![
        vec![
            term([2, 0, 0, 0], 1.0),
            term([0, 1, 0, 0], 1.0),
            term([0, 0, 0, 0], -1.0),
        ],
        vec![
            term([0, 2, 0, 0], 1.0),
            term([0, 0, 1, 0], 1.0),
            term([0, 0, 0, 0], -1.0),
        ],
        vec![
            term([0, 0, 2, 0], 1.0),
            term([0, 0, 0, 1], 1.0),
            term([0, 0, 0, 0], -1.0),
        ],
        vec![
            term([0, 0, 0, 2], 1.0),
            term([1, 0, 0, 0], 1.0),
            term([0, 0, 0, 0], -1.0),
        ],
    ])
    .unwrap();

    // Use a standard (1-homogeneous) partition of the variables.
    let mut part = make_homogeneous_partition(poly.len()).unwrap();

    // Solve the system using a globally convergent homotopy.
    let result = SolveState::solve(
        &mut poly, &mut part, 1.0e-8, 1.0e-8, 1.0e-8, 1, false, false,
    )
    .unwrap();

    println!("Found {} roots:", result.roots.len());
    for (i, root) in result.roots.iter().enumerate() {
        println!("{i:2}: {:?}", &root[..poly.len()]);
    }
}
