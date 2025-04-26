use serde::Deserialize;

#[derive(Deserialize, Debug, Default, Clone)]
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
    /// Bin pileup based on self-identity. Requires fasta in input.
    pub group_by_ani: Option<GroupByANIConfig>,
    /// Minimum size of misassemblies.
    pub minimum_size: Option<MinimumSizeConfig>,
}

impl Config {
    /// Merge two config structs take self as base. Only used for optional config sections.
    pub(crate) fn merge(self, other: Config) -> Self {
        Self {
            general: self.general,
            cov: CoverageConfig {
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
}

impl Default for MinimumSizeConfig {
    fn default() -> Self {
        Self {
            collapse: usize::MIN,
            misjoin: usize::MIN,
            low_quality: usize::MIN,
            false_dupe: 20_000,
            softclip: usize::MIN,
            indel: usize::MIN,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
pub struct GroupByANIConfig {
    /// Size of window to calculate self-identity.
    pub window_size: usize,
    /// Minimum group size.
    pub min_grp_size: usize,
}

impl Default for GroupByANIConfig {
    fn default() -> Self {
        Self {
            window_size: 2000,
            min_grp_size: 50_000,
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
    /// Merge across misassembly type.
    pub merge_across_type: bool,
    /// Store pileup data. Toggle off to reduce memory usage.
    pub store_pileup: bool,
}

impl Default for GeneralConfig {
    fn default() -> Self {
        Self {
            verbose: true,
            bp_merge: 5000,
            bp_wg_window: 10_000_000,
            merge_across_type: false,
            store_pileup: true,
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
    /// Baseline coverage used for false-duplication classification. Defaults to average coverage of region.
    pub baseline: Option<u64>,
    /// Window to apply rolling mean filter. Reduces noise.
    pub rolling_mean_window: Option<usize>,
}

impl Default for CoverageConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 5.0,
            n_zscores_low: 3.0,
            rolling_mean_window: Some(11),
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
            n_zscores_high: 3.4,
            ratio_het: 0.34,
            rolling_mean_window: Some(11),
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
            n_zscores_high: 4.0,
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
            n_zscores_high: 3.4,
            ratio_softclip: 0.5,
        }
    }
}
