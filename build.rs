// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

fn main() {
    cc::Build::new()
        .files(["polsys-plp/polsys_plp.f90", "polsys-plp/lapack_plp.f"])
        .compiler("gfortran")
        .flag("-std=legacy")
        .flag("-Wno-maybe-uninitialized") // suppress the maybe-uninitialized warnings
        .flag("-O3") // optimize level 3
        .static_flag(true)
        .compile("polsys_plp");

    println!("cargo:rustc-link-search=native=polsys_plp");
    println!("cargo:rustc-link-lib=static=polsys_plp");
}
