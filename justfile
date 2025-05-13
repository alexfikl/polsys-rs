PYTHON := "python -X dev"

_default:
    @just --list

# {{{ formatting

[doc("Reformat all source code")]
format: rustfmt fprettify justfmt

[doc("Run rustfmt over all Rust source code")]
rustfmt:
    cargo fmt -- src/*.rs
    cargo fmt -- build.rs
    @echo -e "\e[1;32mrustfmt clean!\e[0m"

[doc("Run fprettify over all Fortran source code")]
fprettify:
    fprettify -i 4 -l 88 -w 3 \
        polsys-plp/polsys_plp_wrapper.f90
    @echo -e "\e[1;32mfprettify clean!\e[0m"

[doc("Run just --fmt over the justfiles")]
justfmt:
    just --unstable --fmt
    @echo -e "\e[1;32mjust --fmt clean!\e[0m"

# }}}
# {{{ linting

[doc("Run all linting checks over the source code")]
lint: typos reuse fortitude clippy

[doc("Run typos over the source code and documentation")]
typos format="brief":
    typos --sort --format={{ format }}
    @echo -e "\e[1;32mtypos clean!\e[0m"

[doc("Check REUSE license compliance")]
reuse:
    {{ PYTHON }} -m reuse lint
    @echo -e "\e[1;32mREUSE compliant!\e[0m"

[doc("Run fortitude lint checks (Fortran)")]
fortitude:
    fortitude check polsys-plp/polsys_plp_wrapper.f90
    @echo -e "\e[1;32m[fort] fortitude clean!\e[0m"

[doc("Run clippy lint checks (Rust)")]
clippy:
    cargo clippy --all-targets --all-features
    @echo -e "\e[1;32m[rust] clippy clean!\e[0m"

# }}}
# {{{ develop

[doc("Update Cargo.lock")]
update:
    cargo update --verbose

[doc("Build project in debug mode")]
debug:
    cargo build --locked --all-features --verbose

[doc("Build project in release mode")]
release:
    cargo build --locked --all-features --release

[doc("Run rust tests")]
test $RUST_BACKTRACE="1":
    cargo test --all-features -- --test-threads 1

[doc("Remove various generated files")]
clean:
    rm -rf *.so
    rm -rf *.mod

[doc("Remove all generated files and caches")]
purge: clean
    rm -rf target

# }}}
