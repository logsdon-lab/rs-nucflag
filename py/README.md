# py_nucflag
Python library for `nucflag`.

## Structure
```
py/
├── Cargo.toml
├── LICENSE
├── pixi.lock
├── py_nucflag.pyi  # Python stub files for type annotations
├── pyproject.toml
├── README.md
├── tests
│   ├── cfg.toml    # Test config file
│   └── test_run.py # Tests
└── src
    ├── lib.rs      # Maturin functions
    └── utils.rs    # Utility functions
```

## Commands
Ensure [`pixi`](https://pixi.prefix.dev/latest/#getting-started) is installed.
```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

Build package and install.
```bash
make build_py
# Check that library is installed
pixi run -m py pip show py_nucflag
```

Test bindings.
```bash
make test_py
```

Rebuild stubs:
```bash
# This will delete the original! Doesn't include classes or types.
make build_stubs
```
