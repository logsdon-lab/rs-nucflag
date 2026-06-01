use coitrees::Interval;
use itertools::Itertools;

/// Merge intervals. Includes book-ended intervals.
///
/// # Arguments
/// * `intervals`: Intervals to merge. Elements are cloned.
/// * `dst`: Distance to merge over.
/// * `fn_cmp`: Function to enforce additional check before merging.
/// * `fn_reducer`: Function to reduce metadata.
///
/// # Returns
/// * Merged overlapping intervals.
pub fn merge_intervals<I, T>(
    intervals: I,
    dst: i32,
    fn_cmp: impl Fn(&Interval<T>, &Interval<T>) -> bool,
    fn_reducer: impl Fn(&Interval<T>, &Interval<T>) -> T,
) -> Vec<Interval<T>>
where
    I: Iterator<Item = Interval<T>>,
    T: Clone,
{
    let mut iter_intervals = intervals.sorted_by(|a, b| a.first.cmp(&b.first));
    let Some(itv_first) = iter_intervals.next() else {
        return vec![];
    };
    let mut merged: Vec<Interval<T>> = vec![itv_first];

    for itv in iter_intervals {
        if let Some(prev) = merged.last_mut().filter(|prev| {
            let dst_between = itv.first - prev.last;
            let added_check = fn_cmp(prev, &itv);
            (dst_between <= dst) & added_check
        }) {
            prev.last = itv.last;
            prev.metadata = fn_reducer(prev, &itv)
        } else {
            merged.push(itv);
        }
    }
    merged
}

// Adapted from https://github.com/chaimleib/intervaltree/blob/1bc406e1f441577c4e421fc51aba2ab67fbd97fb/intervaltree/interval.py#L56
pub fn overlap_length(a_first: i32, a_last: i32, b_first: i32, b_last: i32) -> i32 {
    // No overlap
    if !(a_first < b_last && a_last > b_first) {
        return 0;
    }
    // a  |---|
    // b |---|
    let max_st = std::cmp::max(a_first, b_first);
    let min_end = std::cmp::min(a_last, b_last);
    min_end - max_st
}

#[cfg(test)]
mod tests {
    use super::merge_intervals;
    use coitrees::Interval;
    use std::fmt::Debug;

    fn reduce_to_a<'a>(a: &Interval<usize>, _b: &Interval<usize>) -> usize {
        a.metadata
    }

    fn no_added_check<'a>(_a: &Interval<usize>, _b: &Interval<usize>) -> bool {
        true
    }

    fn assert_itvs_equal<T: Clone + PartialEq + Debug>(
        itvs_1: &[Interval<T>],
        itvs_2: &[Interval<T>],
    ) {
        itertools::assert_equal(
            itvs_1
                .iter()
                .map(|itv| (itv.first, itv.last, itv.metadata.clone())),
            itvs_2
                .iter()
                .map(|itv| (itv.first, itv.last, itv.metadata.clone())),
        );
    }

    #[test]
    fn test_no_merge_intervals() {
        let itvs = vec![
            Interval::new(1, 2, 1),
            Interval::new(3, 5, 2),
            Interval::new(6, 9, 3),
        ];
        let merged_itvs = merge_intervals(itvs.clone().into_iter(), 0, no_added_check, reduce_to_a);
        assert_itvs_equal(&itvs, &merged_itvs);
    }

    #[test]
    fn test_merge_intervals_single() {
        let itvs = vec![
            Interval::new(1, 3, 1),
            Interval::new(3, 5, 2),
            Interval::new(6, 9, 3),
        ];
        let merged_itvs = merge_intervals(itvs.into_iter(), 0, no_added_check, reduce_to_a);
        let exp_itvs = vec![Interval::new(1, 5, 1), Interval::new(6, 9, 3)];

        assert_itvs_equal(&exp_itvs, &merged_itvs);
    }

    #[test]
    fn test_merge_intervals_multiple() {
        let itvs = vec![
            Interval::new(1, 3, 1),
            Interval::new(6, 9, 3),
            Interval::new(3, 6, 2),
        ];
        let merged_itvs = merge_intervals(itvs.into_iter(), 0, no_added_check, reduce_to_a);
        let exp_itvs = vec![Interval::new(1, 9, 1)];
        assert_itvs_equal(&exp_itvs, &merged_itvs);
    }

    #[test]
    fn test_merge_condition() {
        let itvs = vec![
            Interval::new(1, 2, 2),
            Interval::new(3, 4, 2),
            Interval::new(5, 6, 3),
        ];
        let exp_itvs = vec![Interval::new(1, 4, 2), Interval::new(5, 6, 3)];

        let merged_itvs = merge_intervals(
            itvs.clone().into_iter(),
            1,
            |a, b| (a.metadata % 2 == 0) & (b.metadata % 2 == 0),
            reduce_to_a,
        );
        assert_itvs_equal(&merged_itvs, &exp_itvs);
    }
}
