// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use num::complex::Complex64;

use std::ffi::{c_double, c_int};

extern "C" {
    pub fn init_polynomial(
        n: c_int,
        m: c_int,
        n_coeffs_per_eq: *const c_int,
        coefficients: *const Complex64,
        degrees: *const c_int,
        ierr: *mut c_int,
    );

    pub fn init_partition(
        n: c_int,
        m: c_int,
        p: c_int,
        n_sets_per_partition: *const c_int,
        n_indices_per_set: *const c_int,
        indices: *const c_int,
        ierr: *mut c_int,
    );

    pub fn bezout_plp_wrapper(
        n: c_int,
        maxt: c_int,
        tol: c_double,
        bplp: *mut c_int,
        ierr: *mut c_int,
    );

    // pub fn polsys_plp_(
    //     n: *const c_int,
    //     tracktol: *const c_double,
    //     finaltol: *const c_double,
    //     singtol: *mut c_double,
    //     sspar: *mut c_double,
    //     bplp: *mut c_int,
    //     iflag1: *mut c_int,
    //     iflag2: *mut c_int,
    //     arclen: *mut c_double,
    //     lambda: *mut c_double,
    //     roots: *mut c_double_complex,
    //     nfe: *mut c_int,
    //     scale_factors: *mut c_double,
    //     numrr: *const c_int,
    //     recall: *const c_int,
    //     no_scaling: *const c_int,
    //     user_f_df: *const c_int,
    // );
}
