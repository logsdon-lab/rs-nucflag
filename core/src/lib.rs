#![doc=include_str!("../../README.md")]

pub(crate) mod binning;
pub(crate) mod classify;
pub(crate) mod config;
pub(crate) mod intervals;
pub(crate) mod io;
pub(crate) mod misassembly;
pub(crate) mod nucflag;
pub(crate) mod peak;
pub(crate) mod pileup;
pub(crate) mod postprocess;
pub(crate) mod preset;
pub(crate) mod repeats;

pub use classify::NucFlagResult;
pub use config::*;
pub use misassembly::*;
pub use nucflag::nucflag;
pub use preset::*;
pub use repeats::Repeat;
