# py_nucflag
Python library for `nucflag`.

## Structure
```
py/
├── Cargo.toml
├── LICENSE
├── py_nucflag.pyi  # Python stub files for type annotations
├── pyproject.toml
├── README.md
└── src
    ├── lib.rs      # Maturin functions
    └── utils.rs    # Utility functions
```

## Commands
Make venv.
```bash
make venv
```

Build package and install it into venv.
```bash
make build_py
make install_py
```

Rebuild stubs:
```bash
# This will delete the original! Doesn't include classes or types.
make build_py_stubs
```
