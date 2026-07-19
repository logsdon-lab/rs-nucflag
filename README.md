# `rs_nucflag`
[![CI](https://github.com/koisland/rs-nucflag/actions/workflows/ci.yaml/badge.svg)](https://github.com/koisland/rs-nucflag/actions/workflows/ci.yaml)

Library to call misassemblies in genome assemblies from long-read alignments.

## Getting started
To install the `py_nucflag` library.
```bash
pip install py-nucflag
```

To install the `rs_nucflag` library.
```bash
cargo add rs-nucflag
```

> [!NOTE]
> Git LFS required to run tests.

## Project Overview
Three directories:
1. `core`
    * Rust `nucflag` library as `rs-nucflag` crate.
2. `py`
    * Python bindings for `core` as `py-nucflag` library.
3. `examples`
    * Example command-line programs using `rs-nucflag` and `py-nucflag` library.

## Build
Ensure [`pixi`](https://pixi.prefix.dev/latest/#getting-started) is installed for building python bindings.

```bash
# Build py_nucflag.
make build_py
# Build rs-nucflag
make build_rs
```

## Tests
```bash
# Test py_nucflag
make test_py
# Test rs-nucflag
make test_rs
```
