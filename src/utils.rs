use std::num::NonZeroU64;

#[derive(Debug, Clone)]
pub struct Interval<T: Clone + std::fmt::Debug> {
    pub st: NonZeroU64,
    pub end: NonZeroU64,
    pub metadata: T,
}

impl<T: Clone + std::fmt::Debug> Interval<T> {
    pub fn new(st: u64, end: u64, metadata: T) -> eyre::Result<Self> {
        assert!(end >= st, "Invalid interval coordinates.");
        Ok(Self {
            st: st.try_into()?,
            end: end.try_into()?,
            metadata,
        })
    }

    pub fn length(&self) -> u64 {
        self.end.get() - self.st.get()
    }
}
