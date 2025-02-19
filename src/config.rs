use serde::Deserialize;

#[derive(Deserialize, Debug, Default, Clone)]
pub struct Config {
    /// Plot config.
    pub plot: PlotConfig,
    /// First most-common base config.
    pub first: FirstConfig,
    /// Second most-common base config.
    pub second: SecondConfig,
}

#[derive(Deserialize, Debug, Clone)]
/// Config for generated plots.
pub struct PlotConfig {
    /// Dimensions of individual plot.
    pub dim: (u32, u32),
}

impl Default for PlotConfig {
    fn default() -> Self {
        Self { dim: (4000, 750) }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the first-common base coverage signal.
pub struct FirstConfig {
    /// Number of z-scores above the mean to be considered a collapse.
    pub n_zscores_collapse: f32,
    /// Number of z-scores below the mean to be considered a misjoin.
    pub n_zscores_misjoin: f32,
    /// Window size to calculate mean and stdev.
    pub window_size: usize,
    /// Number of bases to merge peaks/valleys.
    pub bp_merge: usize,
    /// Minimum peak/valley size.
    pub min_bp_peak: usize,
}

impl Default for FirstConfig {
    fn default() -> Self {
        Self {
            n_zscores_collapse: 3.4,
            n_zscores_misjoin: 3.4,
            window_size: 500_000,
            bp_merge: 5000,
            min_bp_peak: 1,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the second-common base coverage signal.
pub struct SecondConfig {
    /// Number of z-scores above the mean to flag.
    pub n_zscores: f32,
    /// Window size.
    pub window_size: usize,
    /// Number of bases to merge peaks.
    pub bp_merge: usize,
    /// Threshold to remove background noise signal will have u=0.
    pub min_perc: f32,
    /// Minimum size of peak.
    pub min_bp_peak: usize,
    /// Ratio used to split hets from small collapses.
    pub thr_het_ratio: f32,
}

impl Default for SecondConfig {
    fn default() -> Self {
        Self {
            n_zscores: 3.4,
            window_size: 500_000,
            bp_merge: 5000,
            min_bp_peak: 1,
            min_perc: 0.25,
            thr_het_ratio: 0.2,
        }
    }
}
