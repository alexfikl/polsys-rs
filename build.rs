// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

fn main() {
    cc::Build::new()
        .files([
            "polsys-plp/polsys_plp.f90",
            "polsys-plp/lapack_plp.f",
            "polsys-plp/polsys_plp_wrapper.f90",
        ])
        .compiler("gfortran")
        // .flag("-std=legacy")
        .flag("-Wno-compare-reals")
        .flag("-Wno-do-subscript")
        .flag("-Wno-function-elimination")
        .flag("-Wno-maybe-uninitialized")
        .flag("-Wno-unused-dummy-argument")
        .flag("-Wno-unused-label")
        .flag("-Wno-tabs")
        .flag("-O3")
        .static_flag(true)
        .cargo_metadata(true)
        .compile("polsys_plp");

    println!("cargo::rustc-link-lib=gfortran");
    println!("cargo::rerun-if-changed=polsys-plp/polsys_plp_wrapper.f90");
    println!("cargo::rerun-if-changed=polsys-plp/polsys_plp.f90");
    println!("cargo::rerun-if-changed=polsys-plp/lapack_plp.f");
}
