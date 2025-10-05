.PHONY: build test venv build_py build_rs install_py

BIN=venv/bin/


test:
	cargo test -p rs-nucflag --release


test_remake_images:
	eval $$(cd core/test && rm -f results/*.done && snakemake -p -s regenerate_plots.smk -c 12 > /dev/null);


venv:
	python -m venv venv
	$(BIN)pip install maturin pyo3-stubgen


build_py:
	$(BIN)maturin build --release -m py/Cargo.toml


# This will delete the original! Doesn't include classes or types.
# TODO: Merge files. https://stackoverflow.com/a/9123512
build_py_stubs:
	$(MAKE) install_py
	$(BIN)pyo3-stubgen py_nucflag py


install_py:
	$(BIN)pip install --force-reinstall $(shell find target/wheels -name "*.whl" | sort -r | head -1)


build_rs:
	cargo build --release --manifest-path core/Cargo.toml
