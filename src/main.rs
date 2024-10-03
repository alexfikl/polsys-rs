// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

mod bindings;
mod plp;

use std::collections::HashMap;

use num::complex::c64;

fn main() {
    let mut ierr: i32 = 0;
    let n = 2;
    let n_coeffs_per_eq = [3, 3];
    let coefficients = [
        c64(3.0, 1.3),
        c64(1.0, 0.1),
        c64(-1.0, 0.2),
        c64(2.0, 0.5),
        c64(1.0, 0.5),
        c64(-3.0, 1.0),
    ];
    let degrees = [[2, 0], [0, 2], [0, 0], [1, 1], [0, 2], [0, 0]].as_flattened();

    let mut poly = plp::Polynomial::new(vec![
        HashMap::from([
            ([2, 0], c64(3.0, 1.3)),
            ([0, 2], c64(1.0, 0.1)),
            ([0, 0], c64(-1.0, 0.2)),
        ]),
        HashMap::from([
            ([1, 1], c64(2.0, 0.5)),
            ([0, 2], c64(1.0, 0.5)),
            ([0, 0], c64(-3.0, 1.0)),
        ]),
    ]);
    poly.init().unwrap();
    println!("System size: {}", poly.len());
    println!("System degrees: {:?}", poly.degrees());
    println!("Total degree: {}", poly.total_degree());

    unsafe {
        bindings::init_polynomial(
            n,
            coefficients.len() as i32,
            n_coeffs_per_eq.as_ptr(),
            coefficients.as_ptr(),
            degrees.as_ptr(),
            &mut ierr,
        );
    }

    println!("{:?}", plp::PathTrackingResult::from(2));

    println!("ierr: {}", ierr);

    let mut part = plp::make_homogeneous_partition(3);
    part.init().unwrap();
    println!("n_sets_per_partition: {:?}", part.n_sets_per_partition);

    let mut part = plp::make_m_homogeneous_partition(3, vec![vec![1, 2], vec![3]]);
    part.init().unwrap();
    println!("n_sets_per_partition: {:?}", part.n_sets_per_partition);

    let n_sets_per_partition = [2, 2, 2, 2, 2, 2, 2, 2];
    let n_indices_per_set = [
        [4, 4],
        [4, 4],
        [4, 4],
        [4, 4],
        [4, 4],
        [4, 4],
        [4, 4],
        [4, 4],
    ]
    .as_flattened();
    let indices = [
        [[1, 2, 5, 6], [3, 4, 7, 8]],
        [[1, 2, 5, 6], [3, 4, 7, 8]],
        [[1, 2, 5, 6], [3, 4, 7, 8]],
        [[1, 2, 5, 6], [3, 4, 7, 8]],
        [[1, 2, 5, 6], [3, 4, 7, 8]],
        [[1, 2, 5, 6], [3, 4, 7, 8]],
        [[1, 2, 5, 6], [3, 4, 7, 8]],
        [[1, 2, 5, 6], [3, 4, 7, 8]],
    ]
    .as_flattened()
    .as_flattened();

    unsafe {
        bindings::init_partition(
            n_sets_per_partition.len() as i32,
            2,
            4,
            n_sets_per_partition.as_ptr(),
            n_indices_per_set.as_ptr(),
            indices.as_ptr(),
            &mut ierr,
        );
    }
    println!("ierr: {}", ierr);
}
