# `rs_nucflag`
[![CI](https://github.com/koisland/rs-nucflag/actions/workflows/ci.yaml/badge.svg)](https://github.com/koisland/rs-nucflag/actions/workflows/ci.yaml)

Library to call misassemblies in genome assemblies from long-read alignments.

## Getting started
To install the `py_nucflag` library.
```bash
pip install git+https://github.com/koisland/rs-nucflag.git#subdirectory=py
```

To install the `rs_nucflag` library.
```bash
cargo add --git https://github.com/koisland/rs-nucflag.git --package nucflag
```

> [!NOTE]
> Git LFS required to run tests.

## Project Overview
Three directories:
1. `core`
    * Rust `nucflag` library.
2. `py`
    * Python bindings for `core` as `py_nucflag` library.
3. `examples`
    * Example command-line program using `nucflag` library.

## Build
```bash
make venv && make py_nucflag
```

```bash
make rs_nucflag
```

## Tests
```bash
make test
```
