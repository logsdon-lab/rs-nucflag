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
    pub min_bp: usize,
    /// Merge across misassembly type.
    pub merge_across_type: bool,
    /// Window used if no bedfile provided.
    pub window: usize,
    /// Baseline coverage. Defaults to average coverage of region.
    pub baseline_cov: Option<u64>,
}

impl Default for GeneralConfig {
    fn default() -> Self {
        Self {
            bp_merge: 5000,
            min_bp: 1,
            merge_across_type: false,
            window: 10_000_000,
            baseline_cov: None
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
}

impl Default for FirstConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 3.4,
            n_zscores_low: 3.4,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the second-common base coverage signal.
pub struct SecondConfig {
    /// Number of z-scores above the mean to flag.
    pub n_zscores_high: f32,
    /// Threshold to remove background noise signal will have u=0.
    pub min_perc: f32,
    /// Ratio used to split hets from small collapses.
    pub thr_het_ratio: f32,
}

impl Default for SecondConfig {
    fn default() -> Self {
        Self {
            n_zscores_high: 3.4,
            min_perc: 0.25,
            thr_het_ratio: 0.2,
        }
    }
}
