use std::str::FromStr;

use eyre::bail;
use serde::Deserialize;

use crate::config::{Config, CoverageConfig, IndelConfig, MinimumSizeConfig, MismatchConfig};


/// Sequencing data preset.
#[derive(Deserialize, Debug, Default, Clone)]
pub enum Preset {
    /// PacBio Hifi. Default option.
    #[default]
    PacBioHiFi,
    /// ONT R9. Removes mismatch as signal due to error rate.
    OntR9,
}

impl FromStr for Preset {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "pacbio" | "hifi" | "pacbiohifi" | "pacbio_hifi" => Ok(Preset::PacBioHiFi),
            "ont" | "ontr9" | "ont_r9" | "r9" => Ok(Preset::OntR9),
            _ => bail!("Invalid preset. {s}")
        }
    }
}

impl From<Preset> for Config {
    fn from(value: Preset) -> Self {
        match value {
            Preset::PacBioHiFi => {
                Config::default()
            },
            Preset::OntR9 => {
                Config {
                    mismatch: MismatchConfig {
                        rolling_mean_window: Some(51),
                        ..Default::default()
                    },
                    cov: CoverageConfig {
                        n_zscores_low: 5.0,
                        rolling_mean_window: Some(51),
                        ..Default::default()
                    },
                    indel: IndelConfig {
                        rolling_mean_window: Some(51),
                        ..Default::default()
                    },
                    minimum_size: Some(MinimumSizeConfig {
                        false_dupe: usize::MIN,
                        ..Default::default()
                    }),
                    ..Default::default()
                }
            },
        }
    }
}