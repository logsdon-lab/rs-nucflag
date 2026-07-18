#![doc=include_str!("../../README.md")]

pub(crate) mod binning;
pub(crate) mod intervals;
pub(crate) mod misassembly;
pub(crate) mod nucflag;
pub(crate) mod peak;
pub(crate) mod postprocess;
pub(crate) mod preset;

// Need to be public for py_nucflag
#[doc(hidden)]
pub mod classify;
#[doc(hidden)]
pub mod config;
#[doc(hidden)]
pub mod io;
#[doc(hidden)]
pub mod pileup;
#[doc(hidden)]
pub mod repeats;

pub use classify::NucFlagResult;
#[doc(inline)]
pub use config::*;
#[doc(inline)]
pub use misassembly::*;
pub use nucflag::nucflag;
#[doc(inline)]
pub use preset::*;
pub use repeats::Repeat;
