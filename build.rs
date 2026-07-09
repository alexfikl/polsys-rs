// SPDX-FileCopyrightText: 2024 Alexandru Fikl <alexfikl@gmail.com>
// SPDX-License-Identifier: MIT

fn main() {
    let mut files = vec![
        "polsys-plp/polsys_plp.f90",
        "polsys-plp/polsys_plp_wrapper.f90",
    ];

    let use_system_lapack = std::env::var("CARGO_FEATURE_SYSTEM_LAPACK").as_deref() == Ok("1");

    if use_system_lapack {
        if let Ok(libs) = std::env::var("POLSYS_LAPACK_LIBS") {
            for lib in libs.split(',') {
                println!("cargo::rustc-link-lib={lib}");
            }
        } else {
            pkg_config::probe_library("lapack").unwrap_or_else(|e| {
                panic!(
                    "pkg-config could not find lapack ({e}). \
                     Set POLSYS_LAPACK_LIBS=openblas (or lapack,blas, etc.)"
                )
            });
        }
    } else {
        files.push("polsys-plp/lapack_plp.f");
    }

    cc::Build::new()
        .files(&files)
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
        .cargo_metadata(true)
        .compile("polsys_plp");

    println!("cargo::rustc-link-lib=gfortran");
    println!("cargo::rerun-if-changed=polsys-plp/polsys_plp_wrapper.f90");
    println!("cargo::rerun-if-changed=polsys-plp/polsys_plp.f90");
    println!("cargo::rerun-if-changed=polsys-plp/lapack_plp.f");
}
