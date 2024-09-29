// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

mod bindings;

use bindings::c64;

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

    unsafe {
        bindings::init_polynomial(
            n,
            coefficients.len() as i32,
            n_coeffs_per_eq.as_ptr(),
            coefficients.as_ptr(),
            degrees.as_ptr(),
            &mut ierr,
        );

        bindings::init_polynomial(
            n,
            coefficients.len() as i32,
            n_coeffs_per_eq.as_ptr(),
            coefficients.as_ptr(),
            degrees.as_ptr(),
            &mut ierr,
        );
    }

    println!("ierr: {}", ierr);

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
