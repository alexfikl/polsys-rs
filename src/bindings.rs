// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

use std::ffi::{c_double, c_int};

#[repr(C)]
pub struct __BindgenComplex<T> {
    pub re: T,
    pub im: T,
}

pub type c_double_complex = __BindgenComplex<f64>;

pub type c_float_complex = __BindgenComplex<f32>;

extern "C" {
    pub fn bezout_plp_(n: *const c_int, maxt: *const c_int, tol: *const c_double, bplp: *mut c_int);

    pub fn polsys_plp_(
        n: *const c_int,
        tracktol: *const c_double,
        finaltol: *const c_double,
        singtol: *mut c_double,
        sspar: *mut c_double,
        bplp: *mut c_int,
        iflag1: *mut c_int,
        iflag2: *mut c_int,
        arclen: *mut c_double,
        lambda: *mut c_double,
        roots: *mut c_double_complex,
        nfe: *mut c_int,
        scale_factors: *mut c_double,
        numrr: *const c_int,
        recall: *const c_int,
        no_scaling: *const c_int,
        user_f_df: *const c_int,
    );
}
