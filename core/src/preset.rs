use std::str::FromStr;

use eyre::bail;
use serde::Deserialize;

use crate::config::{
    Config, CoverageConfig, GeneralConfig, IndelConfig, MinimumSizeConfig, MismatchConfig,
};

/// Sequencing data preset.
#[derive(Deserialize, Debug, Default, Clone)]
pub enum Preset {
    /// PacBio Hifi. Default option. Accuracy ~99.9%.
    #[default]
    PacBioHiFi,
    /// ONT R9. Smooths mismatch signal due to error rate. Removes small indels (<20 bp). Accuracy of ~95%.
    OntR9,
    /// ONT R10. Removes small indels (<10 bp). Accuracy of ~99%.
    OntR10,
}

impl FromStr for Preset {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "pacbio" | "hifi" | "pacbiohifi" | "pacbio_hifi" => Ok(Preset::PacBioHiFi),
            "ont" | "ontr9" | "ont_r9" | "r9" => Ok(Preset::OntR9),
            "ontr10" | "ont_r10" | "r10" => Ok(Preset::OntR10),
            _ => bail!("Invalid preset. {s}"),
        }
    }
}

impl From<Preset> for Config {
    fn from(value: Preset) -> Self {
        match value {
            Preset::PacBioHiFi => Config::default(),
            Preset::OntR9 => Config {
                general: GeneralConfig {
                    // Remove small reads that will mismap due to low base accuracy.
                    bp_min_aln_length: 30_000,
                    ..Default::default()
                },
                mismatch: MismatchConfig {
                    rolling_mean_window: Some(25),
                    ..Default::default()
                },
                cov: CoverageConfig {
                    // Dips tend to be less prominent due to avg length or read.
                    n_zscores_high: 4.0,
                    n_zscores_low: 1.0,
                    ratio_collapse: 1.5,
                    rolling_mean_window: Some(11),
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
            Preset::OntR10 => Config {
                cov: CoverageConfig {
                    n_zscores_high: 4.0,
                    n_zscores_low: 1.0,
                    ratio_collapse: 1.5,
                    rolling_mean_window: Some(15),
                    ..Default::default()
                },
                indel: IndelConfig {
                    rolling_mean_window: None,
                    min_ins_size: 10,
                    min_del_size: 10,
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
