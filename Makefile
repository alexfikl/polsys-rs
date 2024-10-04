PYTHON?=python -X dev

all: help

help: 			## Show this help
	@echo -e "Specify a command. The choices are:\n"
	@grep -E '^[0-9a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[0;36m%-12s\033[m %s\n", $$1, $$2}'
	@echo ""
.PHONY: help

# {{{ formatting

format: rustfmt fprettify		## Run all formatting scripts
.PHONY: format

rustfmt:						## Run rustfmt
	cargo fmt -- --config 'max_width=88' src/*.rs
	@echo -e "\e[1;32mrustfmt clean!\e[0m"
.PHONY: rustfmt

fprettify: 						## Run fprettify
	fprettify -i 4 -l 88 -w 3 \
		polsys-plp/polsys_plp_wrapper.f90
	@echo -e "\e[1;32mfprettify clean!\e[0m"
.PHONY: fprettify

# }}}

# {{{ linting

lint: typos reuse fortitude clippy		## Run linting checks
.PHONY: lint

typos:			## Run typos over the source code and documentation
	typos --sort
	@echo -e "\e[1;32mtypos clean!\e[0m"
.PHONY: typos

reuse:			## Check REUSE license compliance
	$(PYTHON) -m reuse lint
	@echo -e "\e[1;32mREUSE compliant!\e[0m"
.PHONY: reuse

clippy:			## Run clippy lint checks (Rust)
	cargo clippy --all-targets --all-features
	@echo -e "\e[1;32mclippy clean!\e[0m"
.PHONY: clippy

fortitude:		## Run fortitude link checks (Fortran)
	fortitude check --line-length 88 \
		--ignore T032 \
		polsys-plp/polsys_plp_wrapper.f90
	@echo -e "\e[1;32mclippy clean!\e[0m"
.PHONY: fortitude

# }}}

# {{{ building

update:			## Update lock file
	cargo update
.PHONY: update

build:			## Build the project in debug mode
	cargo build --locked --all-features --verbose
.PHONY: build

release:		## Build project in release mode
	cargo build --locked --all-features --release
.PHONY: release

test:			## Run tests
	cargo test
.PHONY: test

# }}}

clean:			## Remove temporary files
	rm -rf *.mod
.PHONY: clean

purge: clean	## Remove all generated files
	rm -rf target
.PHONY: purge
