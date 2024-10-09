// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use num::complex::Complex64;

extern "C" {
    pub fn init_polynomial(
        n: i32,
        m: i32,
        n_coeffs_per_eq: *const i32,
        coefficients: *const Complex64,
        degrees: *const i32,
        ierr: *mut i32,
    );

    pub fn is_polynomial_allocated(flag: *mut i32);

    pub fn deallocate_polynomial(ierr: *mut i32);

    pub fn init_partition(
        n: i32,
        m: i32,
        p: i32,
        n_sets_per_partition: *const i32,
        n_indices_per_set: *const i32,
        indices: *const i32,
        ierr: *mut i32,
    );

    pub fn is_partition_allocated(flag: *mut i32);

    pub fn deallocate_partition(ierr: *mut i32);

    pub fn bezout_plp_wrap(n: i32, maxt: i32, tol: f64, bplp: *mut i32, ierr: *mut i32);

    pub fn polsys_plp_wrap(
        n: i32,
        tracktol: f64,
        finaltol: f64,
        singtol: f64,
        sspar: *mut f64,
        bplp: i32,
        iflag1: *mut i32,
        iflag2: *mut i32,
        arclen: *mut f64,
        lambda: *mut f64,
        roots: *mut Complex64,
        nfe: *mut i32,
        scale_factors: *mut f64,
        numrr: i32,
        recall: i32,
        no_scaling: i32,
        user_f_df: i32,
    );
}
