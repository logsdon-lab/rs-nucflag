use std::{fs::File, path::Path};

use polars::prelude::*;

pub fn write_tsv(df: &mut DataFrame, path: impl AsRef<Path>) -> eyre::Result<()> {
    let mut file = File::create(path)?;
    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(df)?;
    Ok(())
}
