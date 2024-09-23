// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

mod bindings;

use bindings::c64;

fn main() {
    let mut ierr: i32 = 0;
    let n = 2;
    let n_coeffs_per_eq = [3, 3];
    let coefficients = [
        c64(3.0, 0.0),
        c64(1.0, 0.0),
        c64(-1.0, 0.0),
        c64(2.0, 0.0),
        c64(1.0, 0.0),
        c64(-3.0, 0.0),
    ];
    let degrees = [[2, 0], [0, 2], [0, 0], [1, 1], [0, 2], [0, 0]];

    unsafe {
        bindings::init_polynomial(
            &n,
            n_coeffs_per_eq.as_ptr(),
            coefficients.as_ptr(),
            degrees.as_flattened().as_ptr(),
            &mut ierr,
        );
    }
    println!("ierr: {}", ierr);
    // unsafe {
    //     bindings::init_partition(&mut ierr);
    // }
    // println!("ierr: {}", ierr);
}
