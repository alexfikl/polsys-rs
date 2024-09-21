fn main(){
    cc::Build::new()
        .files(["polsys-plp/polsys_plp.f90", "polsys-plp/lapack_plp.f"])
        .compiler("gfortran")
        .flag("-std=legacy")
        .flag("-fdefault-real-8") // use 8 bytes for all floats
        .flag("-Wno-maybe-uninitialized") // suppress the maybe-unitialized warnings
        .flag("-O2") // optimize level 3
        .compile("polsys_plp");

    println!("cargo:rustc-link-search=native=polsys_plp");
    println!("cargo:rustc-link-lib=static=polsys_plp");
}
