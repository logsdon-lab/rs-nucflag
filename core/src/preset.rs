use std::str::FromStr;

use eyre::bail;
use serde::Deserialize;

use crate::config::{Config, CoverageConfig, IndelConfig, MinimumSizeConfig, MismatchConfig};

/// Sequencing data preset.
#[derive(Deserialize, Debug, Default, Clone)]
pub enum Preset {
    /// PacBio Hifi. Default option. Accuracy ~99.9%.
    #[default]
    PacBioHiFi,
    /// ONT R9. Smooths mismatch as signal due to error rate. Accuracy of ~95%.
    OntR9,
}

impl FromStr for Preset {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "pacbio" | "hifi" | "pacbiohifi" | "pacbio_hifi" => Ok(Preset::PacBioHiFi),
            "ont" | "ontr9" | "ont_r9" | "r9" => Ok(Preset::OntR9),
            _ => bail!("Invalid preset. {s}"),
        }
    }
}

impl From<Preset> for Config {
    fn from(value: Preset) -> Self {
        match value {
            Preset::PacBioHiFi => Config::default(),
            Preset::OntR9 => Config {
                mismatch: MismatchConfig {
                    rolling_mean_window: Some(31),
                    ..Default::default()
                },
                cov: CoverageConfig {
                    n_zscores_high: 3.0,
                    n_zscores_low: 3.0,
                    ratio_misjoin: 0.2,
                    ratio_collapse: 1.5,
                    rolling_mean_window: Some(31),
                    ..Default::default()
                },
                indel: IndelConfig {
                    rolling_mean_window: None,
                    min_ins_size: 20,
                    min_del_size: 20,
                    ..Default::default()
                },
                minimum_size: Some(MinimumSizeConfig {
                    false_dupe: 2,
                    ..Default::default()
                }),
                ..Default::default()
            },
        }
    }
}
