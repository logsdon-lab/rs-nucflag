use serde::Deserialize;

#[derive(Deserialize, Debug, Default, Clone)]
pub struct Config {
    /// General config.
    pub general: GeneralConfig,
    /// First most-common base config.
    pub first: FirstConfig,
    /// Second most-common base config.
    pub second: SecondConfig,
}

#[derive(Deserialize, Debug, Clone)]
/// Config for generated plots.
pub struct GeneralConfig {
    /// Number of bases to merge intervals.
    pub bp_merge: usize,
    /// Minimum misassembly size.
    pub bp_min: usize,
    /// Window in base pairs. Only used if no BED file is provided.
    pub bp_window: usize,
    /// Baseline coverage used for false-duplication classification. Defaults to average coverage of region.
    pub cov: Option<u64>,
    /// Store coverage data. Toggle off to reduce memory usage.
    pub store_coverage: bool,
    /// Merge across misassembly type.
    pub merge_across_type: bool,
}

impl Default for GeneralConfig {
    fn default() -> Self {
        Self {
            bp_merge: 5000,
            bp_min: 1,
            bp_window: 10_000_000,
            merge_across_type: false,
            cov: None,
            store_coverage: true,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the first-common base coverage signal.
pub struct FirstConfig {
    /// Number of z-scores above the mean to be considered a misassembly.
    pub n_zscores_high: f32,
    /// Number of z-scores below the mean to be considered a misassembly.
    pub n_zscores_low: f32,
    /// Number of z-scores below the mean to be considered a false-dupe.
    pub n_zscores_false_dupe: f32,
}

impl Default for FirstConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 4.0,
            n_zscores_low: 4.0,
            n_zscores_false_dupe: 2.0,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the second-common base coverage signal.
pub struct SecondConfig {
    /// Number of z-scores above the mean to flag.
    pub n_zscores_high: f32,
    /// Ratio used to split hets from small collapses.
    pub ratio_het: f32,
}

impl Default for SecondConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 3.4,
            ratio_het: 0.2,
        }
    }
}
