// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use std::ffi::c_int;

#[repr(C)]
pub struct __BindgenComplex<T> {
    pub re: T,
    pub im: T,
}

#[allow(non_camel_case_types)]
pub type c_double_complex = __BindgenComplex<f64>;

pub fn c64<T: Into<f64>>(re: T, im: T) -> c_double_complex {
    __BindgenComplex {
        re: re.into(),
        im: im.into(),
    }
}

extern "C" {
    pub fn init_polynomial(
        n: c_int,
        m: c_int,
        n_coeffs_per_eq: *const c_int,
        coefficients: *const c_double_complex,
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

    // pub fn bezout_plp_(
    //     n: *const c_int,
    //     maxt: *const c_int,
    //     tol: *const c_double,
    //     bplp: *mut c_int,
    // );

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
