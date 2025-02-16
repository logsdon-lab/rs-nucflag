use coitrees::{GenericInterval, Interval};

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Peak {
    Low,
    High,
}

pub struct PeaksDetector {
    threshold: f64,
    influence: f64,
    window: Vec<f64>,
}

impl PeaksDetector {
    pub fn new(lag: usize, threshold: f64, influence: f64) -> PeaksDetector {
        PeaksDetector {
            threshold,
            influence,
            window: Vec::with_capacity(lag),
        }
    }

    pub fn signal(&mut self, value: f64) -> Option<Peak> {
        if self.window.len() < self.window.capacity() {
            self.window.push(value);
            None
        } else if let (Some((mean, stddev)), Some(&window_last)) =
            (self.stats(), self.window.last())
        {
            self.window.remove(0);
            if (value - mean).abs() > (self.threshold * stddev) {
                let next_value = (value * self.influence) + ((1. - self.influence) * window_last);
                self.window.push(next_value);
                Some(if value > mean { Peak::High } else { Peak::Low })
            } else {
                self.window.push(value);
                None
            }
        } else {
            None
        }
    }

    pub fn stats(&self) -> Option<(f64, f64)> {
        if self.window.is_empty() {
            None
        } else {
            let window_len = self.window.len() as f64;
            let mean = self.window.iter().fold(0., |a, v| a + v) / window_len;
            let sq_sum = self
                .window
                .iter()
                .fold(0., |a, v| a + ((v - mean) * (v - mean)));
            let stddev = (sq_sum / window_len).sqrt();
            Some((mean, stddev))
        }
    }
}

pub struct PeaksIter<I, F> {
    source: I,
    signal: F,
    detector: PeaksDetector,
}

pub trait PeaksFilter<I>
where
    I: Iterator,
{
    fn peaks<F>(self, detector: PeaksDetector, signal: F) -> PeaksIter<I, F>
    where
        F: FnMut(&I::Item) -> f64;
}

impl<I> PeaksFilter<I> for I
where
    I: Iterator,
{
    fn peaks<F>(self, detector: PeaksDetector, signal: F) -> PeaksIter<I, F>
    where
        F: FnMut(&I::Item) -> f64,
    {
        PeaksIter {
            source: self,
            signal,
            detector,
        }
    }
}

impl<I, F> Iterator for PeaksIter<I, F>
where
    I: Iterator,
    F: FnMut(&I::Item) -> f64,
{
    type Item = (I::Item, Peak);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(item) = self.source.next() {
            let value = (self.signal)(&item);
            if let Some(peak) = self.detector.signal(value) {
                return Some((item, peak));
            }
        }
        None
    }
}

pub fn merge_peaks(
    peaks: impl Iterator<Item = (usize, Peak)>,
    bp_merge: Option<usize>,
    thr_bp_peak: Option<usize>,
) -> eyre::Result<Vec<Interval<Peak>>> {
    let mut final_peaks: Vec<Interval<Peak>> = vec![];
    let mut curr_peak: Vec<(usize, Peak)> = vec![];
    let bp_merge = bp_merge.unwrap_or(5000);
    let thr_bp_peak = thr_bp_peak.unwrap_or(1).try_into()?;

    fn build_interval(curr_peak: &Vec<(usize, Peak)>) -> eyre::Result<Interval<Peak>> {
        let (curr_peak_first, curr_peak_last) =
            (curr_peak.first().unwrap(), curr_peak.last().unwrap());
        Ok(Interval::new(
            curr_peak_first.0.try_into()?,
            curr_peak_last.0.try_into()?,
            curr_peak_first.1,
        ))
    };

    for pk in peaks {
        if let Some(prev_pk) = curr_peak.last().filter(|prev_pk| prev_pk.1 == pk.1) {
            let dst_diff = pk.0 - prev_pk.0;
            // Greater than merge distance, add current peak elements. Then store new peak.
            if dst_diff > bp_merge {
                final_peaks.push(build_interval(&curr_peak)?);
                curr_peak.clear();
            }
        }
        curr_peak.push(pk);
    }
    if !curr_peak.is_empty() {
        final_peaks.push(build_interval(&curr_peak)?);
        curr_peak.clear();
    }
    final_peaks.retain(|pk| pk.len() > thr_bp_peak);
    Ok(final_peaks)
}
