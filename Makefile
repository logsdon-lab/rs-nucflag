.PHONY: build_rs test_rs test_rs_remake_images test_py build_py build_py_stubs

test_rs:
	cargo test -p rs-nucflag --release

build_rs:
	cargo build --release --manifest-path core/Cargo.toml

test_rs_remake_images:
	eval $$(cd core/test && rm -f results/*.done && snakemake -p -s regenerate_plots.smk -c 12 > /dev/null);

test_py:
	pixi run -m py/ test

build_py:
	pixi run -m py/ build

# This will delete the original! Doesn't include classes or types.
# TODO: Merge files. https://stackoverflow.com/a/9123512
build_py_stubs:
	pixi run -m py/ test
