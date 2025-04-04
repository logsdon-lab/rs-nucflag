.PHONY: build test venv build_py build_rs install_py

BIN=venv/bin/


test:
	cargo test


venv:
	python -m venv venv
	$(BIN)pip install maturin


build_py:
	$(BIN)maturin build --release -m py/Cargo.toml


install_py:
	$(BIN)pip install --force-reinstall $(shell find target/wheels -name "*.whl" | sort -r | head -1)


build_rs:
	cargo build --release --manifest-path core/Cargo.toml
