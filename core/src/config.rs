use std::collections::{HashMap, HashSet};

use serde::Deserialize;

use crate::misassembly::MisassemblyType;

#[derive(Deserialize, Debug, Clone)]
pub struct Config {
    /// General config.
    pub general: GeneralConfig,
    /// Coverage config.
    pub cov: CoverageConfig,
    /// Mismatch base config.
    pub mismatch: MismatchConfig,
    /// Indel base config.
    pub indel: IndelConfig,
    /// Softclip base config.
    pub softclip: SoftClipConfig,
    /// Repeat detection config.
    pub repeat: Option<RepeatConfig>,
    /// Bin pileup based on self-identity. Requires fasta in input.
    pub group_by_ani: Option<GroupByANIConfig>,
    /// Minimum size of misassemblies.
    pub minimum_size: Option<MinimumSizeConfig>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            general: GeneralConfig::default(),
            cov: CoverageConfig::default(),
            mismatch: MismatchConfig::default(),
            indel: IndelConfig::default(),
            softclip: SoftClipConfig::default(),
            group_by_ani: Some(GroupByANIConfig::default()),
            minimum_size: Some(MinimumSizeConfig::default()),
            repeat: Some(RepeatConfig::default()),
        }
    }
}

impl Config {
    /// Merge two config structs take self as base. Only used for optional config sections.
    pub(crate) fn merge(self, other: Config) -> Self {
        Self {
            general: self.general,
            cov: CoverageConfig {
                n_zscores_high: other.cov.n_zscores_high,
                n_zscores_low: other.cov.n_zscores_low,
                ratio_misjoin: other.cov.ratio_misjoin,
                rolling_mean_window: other.cov.rolling_mean_window,
                ..self.cov
            },
            mismatch: MismatchConfig {
                rolling_mean_window: other.mismatch.rolling_mean_window,
                ..self.mismatch
            },
            indel: IndelConfig {
                rolling_mean_window: other.indel.rolling_mean_window,
                ..self.indel
            },
            softclip: self.softclip,
            group_by_ani: self.group_by_ani,
            minimum_size: other.minimum_size,
            repeat: self.repeat,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
pub struct MinimumSizeConfig {
    pub collapse: usize,
    pub misjoin: usize,
    pub low_quality: usize,
    pub false_dupe: usize,
    pub softclip: usize,
    pub indel: usize,
    pub homopolymer: usize,
    pub simple_repeat: usize,
    pub scaffold: usize,
}

impl TryFrom<&MinimumSizeConfig> for HashMap<&str, u64> {
    type Error = eyre::Error;

    fn try_from(cfg: &MinimumSizeConfig) -> Result<Self, Self::Error> {
        Ok(HashMap::from_iter([
            ("good", 1),
            ("collapse", cfg.collapse.try_into()?),
            ("false_dupe", cfg.false_dupe.try_into()?),
            ("indel", cfg.indel.try_into()?),
            ("low_quality", cfg.low_quality.try_into()?),
            ("misjoin", cfg.misjoin.try_into()?),
            ("softclip", cfg.softclip.try_into()?),
            ("homopolymer", cfg.homopolymer.try_into()?),
            ("simple_repeat", cfg.simple_repeat.try_into()?),
            ("scaffold", cfg.scaffold.try_into()?),
        ]))
    }
}

impl Default for MinimumSizeConfig {
    fn default() -> Self {
        Self {
            collapse: 100,
            misjoin: 1,
            low_quality: 1,
            false_dupe: 20_000,
            softclip: 1,
            indel: 1,
            homopolymer: 3,
            simple_repeat: 4,
            scaffold: 1,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
pub struct GroupByANIConfig {
    /// Size of window to calculate self-identity.
    pub window_size: usize,
    /// Minimum group size.
    /// * Smaller sizes may result in more false-positives due to coverage changes in transition regions.
    pub min_grp_size: usize,
    /// Minimum identity of group.
    pub min_ident: f32,
}

impl Default for GroupByANIConfig {
    fn default() -> Self {
        Self {
            window_size: 5_000,
            min_grp_size: 100_000,
            min_ident: 80.0,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Config for generated plots.
pub struct GeneralConfig {
    /// Display verbose logging.
    pub verbose: bool,
    /// Number of bases to merge misassembly intervals.
    pub bp_merge: usize,
    /// Whole genome window size in base pairs. Only used if no BED file is provided.
    pub bp_wg_window: usize,
}

impl Default for GeneralConfig {
    fn default() -> Self {
        Self {
            verbose: true,
            bp_merge: 50_000,
            bp_wg_window: 10_000_000,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for coverage signal.
pub struct CoverageConfig {
    /// Number of z-scores above the median to be considered a misassembly.
    pub n_zscores_high: f32,
    /// Number of z-scores below the median to be considered a misassembly.
    pub n_zscores_low: f32,
    /// Minimum coverage ratio required for a misjoin.
    pub ratio_misjoin: f32,
    /// Minimum coverage ratio required for a collapse.
    pub ratio_collapse: f32,
    /// Minimum coverage ratio required for a false dupe.
    pub ratio_false_dupe: f32,
    /// Baseline coverage used for false-duplication classification. Defaults to average coverage of region.
    pub baseline: Option<u32>,
    /// Window to apply rolling mean filter. Reduces noise.
    pub rolling_mean_window: Option<usize>,
}

impl Default for CoverageConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 3.5,
            n_zscores_low: 3.5,
            ratio_misjoin: 0.01,
            ratio_collapse: 1.5,
            ratio_false_dupe: 0.5,
            rolling_mean_window: None,
            baseline: None,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the mismatch signal.
pub struct MismatchConfig {
    /// Number of z-scores above the median to flag.
    pub n_zscores_high: f32,
    /// Ratio used to split hets from small collapses.
    pub ratio_het: f32,
    /// Window to apply rolling mean filter. Reduces noise.
    pub rolling_mean_window: Option<usize>,
}

impl Default for MismatchConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 3.5,
            ratio_het: 0.2,
            rolling_mean_window: None,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the base indel coverage signal.
pub struct IndelConfig {
    /// Number of z-scores above the median to flag.
    pub n_zscores_high: f32,
    /// Ratio used to call indel peaks.
    pub ratio_indel: f32,
    /// Window to apply rolling mean filter. Reduces noise.
    pub rolling_mean_window: Option<usize>,
}

impl Default for IndelConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 10.0,
            ratio_indel: 0.5,
            rolling_mean_window: Some(11),
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the base softclip coverage signal.
pub struct SoftClipConfig {
    /// Number of z-scores above the median to flag.
    pub n_zscores_high: f32,
    /// Ratio used to call softclipped peaks.
    pub ratio_softclip: f32,
}

impl Default for SoftClipConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 3.5,
            ratio_softclip: 0.5,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for repeat detection from misassemblies. Requires providing fasta.
pub struct RepeatConfig {
    /// Which misassembles to check for repeats.
    /// * Usually types associated with drops in coverage.
    /// * Defaults to [`MisassemblyType::Misjoin`] and [`MisassemblyType::Indel`].
    pub check_types: HashSet<MisassemblyType>,
    /// Ratio required of checked region to call as repeat.
    /// * Defaults to a majority.
    pub ratio_repeat: f32,
    /// Extend region checked by n bases on both ends.
    /// By default is the misassembled regions length.
    /// * Sometimes this is not enough is only 1 position long.
    /// * Defaults to 2 bp.  
    pub bp_extend: usize,
}

impl Default for RepeatConfig {
    fn default() -> Self {
        Self {
            check_types: HashSet::from_iter([MisassemblyType::Misjoin, MisassemblyType::Indel]),
            ratio_repeat: 0.5,
            bp_extend: 2,
        }
    }
}
