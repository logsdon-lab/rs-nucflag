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
    /// Dimensions of individual plot.
    pub plot_dim: (u32, u32),
    /// Number of bases to merge intervals.
    pub bp_merge: usize,
    /// Minimum misassembly size.
    pub min_bp: usize,
}

impl Default for GeneralConfig {
    fn default() -> Self {
        Self {
            plot_dim: (4000, 750),
            bp_merge: 5000,
            min_bp: 1,
        }
    }
}

#[derive(Deserialize, Debug, Clone)]
/// Configuration for the first-common base coverage signal.
pub struct FirstConfig {
    /// Number of z-scores above the mean to be considered a misassembly.
    /// Absolute value. For both collapses and misjoins
    pub n_zscores: f32,
    /// Window size to calculate mean and stdev.
    pub window_size: usize,
}

impl Default for FirstConfig {
    fn default() -> Self {
        Self {
            n_zscores: 3.4,
            window_size: 500_000,
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
    /// Threshold to remove background noise signal will have u=0.
    pub min_perc: f32,
    /// Ratio used to split hets from small collapses.
    pub thr_het_ratio: f32,
}

impl Default for SecondConfig {
    fn default() -> Self {
        Self {
            n_zscores: 3.4,
            window_size: 500_000,
            min_perc: 0.25,
            thr_het_ratio: 0.2,
        }
    }
}
