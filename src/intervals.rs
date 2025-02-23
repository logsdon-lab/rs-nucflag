use coitrees::{COITree, Interval, IntervalTree};
use itertools::Itertools;
use std::collections::VecDeque;

#[allow(unused)]
/// Get overlapping intervals in a [`COITree`]. Simply wraps [`COITree::query`].
///
/// # Arguments
/// * `itree`: Interval tree to check for overlaps.
/// * `start`: Start position.
/// * `stop`: Stop position.
///
/// # Returns
/// * [`Vec`] of intervals as a tuple.
///     * [`COITree`]'s `IntervalNode` doesn't impl [`Debug`].
pub fn get_overlapping_intervals<T>(
    itree: &COITree<T, usize>,
    start: i32,
    stop: i32,
) -> Vec<Interval<T>>
where
    T: Clone,
{
    let mut overlapping_itvs = vec![];
    itree.query(start, stop, |n| {
        overlapping_itvs.push(Interval::new(n.first, n.last, n.metadata.clone()));
    });
    overlapping_itvs
}

/// Merge intervals in a [`COITree`]. Includes book-ended intervals.
///
/// # Arguments
/// * `intervals`: Intervals to merge. Elements are cloned.
/// * `data_reducer`: Function to reduce metadata.
/// * `data_finalizer`: Function to apply some final operation on intervals.
///
/// # Returns
/// * Merged overlapping intervals.
pub fn merge_overlapping_intervals<I, T>(
    intervals: I,
    data_reducer: impl Fn(&Interval<T>, &Interval<T>) -> T,
    data_finalizer: impl Fn(Interval<T>) -> Interval<T>,
) -> Vec<Interval<T>>
where
    I: Iterator<Item = Interval<T>>,
    T: Clone,
{
    // let mut merged: Vec<Interval<T>> = Vec::with_capacity(intervals.len());
    let mut merged: Vec<Interval<T>> = Vec::new();
    let mut intervals: VecDeque<Interval<T>> = intervals
        .into_iter()
        .sorted_by(|a, b| a.first.cmp(&b.first))
        .collect();
    while !intervals.is_empty() {
        let Some(itv_1) = intervals.pop_front() else {
            unreachable!()
        };
        let Some(itv_2) = intervals.pop_front() else {
            merged.push(itv_1);
            break;
        };
        // (if) First case:
        // 1-2
        //     3-4
        // (else) Second case:
        // 1-2
        //   2-3
        // (else) Third case:
        // 1-2
        // 1-2
        if itv_1.last < itv_2.first {
            merged.push(itv_1);
            intervals.push_front(itv_2);
        } else {
            let new_data = data_reducer(&itv_1, &itv_2);
            let merged_interval = Interval::new(itv_1.first, itv_2.last, new_data);
            intervals.push_front(merged_interval);
        }
    }
    // Apply finalizer function
    merged.into_iter().map(data_finalizer).collect_vec()
}

/// Subtract interval from st and end position
pub fn subtract_interval<T: Clone>(itvs: &[Interval<T>], st: i32, end: i32) -> Vec<Interval<T>> {
    let mut sub_itvs = vec![];

    let itree: COITree<T, usize> = COITree::new(itvs);
    let overlapping_itvs  = get_overlapping_intervals(&itree, st, end);
    let mut sort_itvs = overlapping_itvs.into_iter().sorted_by(|a, b| a.first.cmp(&b.first)).peekable();

    if let Some(first_itv) = sort_itvs.peek().filter(|itv| itv.first > st) {
        sub_itvs.push(Interval::new(st, first_itv.first, first_itv.metadata.clone()));
    }

    while let Some(itv) = sort_itvs.next() {
        let Some(next_itv) = sort_itvs.peek() else {
            if itv.last < end {
                sub_itvs.push(Interval::new(itv.last, end, itv.metadata.clone()));
            }
            break;
        };
        sub_itvs.push(Interval::new(itv.last, next_itv.first, itv.metadata));
    }
    sub_itvs
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use coitrees::Interval;

    use crate::intervals::subtract_interval;

    use super::merge_overlapping_intervals;

    fn reduce_to_a<'a>(a: &Interval<usize>, _b: &Interval<usize>) -> usize {
        a.metadata
    }

    fn noop(a: Interval<usize>) -> Interval<usize> {
        a
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
        let merged_itvs = merge_overlapping_intervals(itvs.clone().into_iter(), reduce_to_a, noop);
        assert_itvs_equal(&itvs, &merged_itvs);
    }

    #[test]
    fn test_merge_intervals_single() {
        let itvs = vec![
            Interval::new(1, 3, 1),
            Interval::new(3, 5, 2),
            Interval::new(6, 9, 3),
        ];
        let merged_itvs = merge_overlapping_intervals(itvs.into_iter(), reduce_to_a, noop);
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
        let merged_itvs = merge_overlapping_intervals(itvs.into_iter(), reduce_to_a, noop);
        let exp_itvs = vec![Interval::new(1, 9, 1)];
        assert_itvs_equal(&exp_itvs, &merged_itvs);
    }

    #[test]
    fn test_subtract_itv() {
        let itvs = [
            Interval::new(1, 3, ()),
            Interval::new(4, 5, ()),
            Interval::new(8, 12, ()),
            Interval::new(12, 20, ()),
        ];

        let sub_itvs = subtract_interval(&itvs, 1, 20);
        assert_itvs_equal(
            &sub_itvs,
            &[
                Interval::new(3, 4, ()),
                Interval::new(5, 8, ()),
                Interval::new(12, 12, ()),
            ]
        )
    }

    #[test]
    fn test_subtract_itv_last_remainder() {
        let itvs = [
            Interval::new(1, 3, ()),
            Interval::new(4, 5, ()),
            Interval::new(8, 12, ()),
        ];

        let sub_itvs = subtract_interval(&itvs, 1, 20);
        assert_itvs_equal(
            &sub_itvs,
            &[
                Interval::new(3, 4, ()),
                Interval::new(5, 8, ()),
                Interval::new(12, 20, ()),
            ]
        )
    }
    #[test]
    fn test_subtract_itv_first_remainder() {
        let itvs = [
            Interval::new(1, 3, ()),
            Interval::new(4, 5, ()),
            Interval::new(8, 12, ()),
        ];

        let sub_itvs = subtract_interval(&itvs, 10, 20);
        assert_itvs_equal(
            &sub_itvs,
            &[
                Interval::new(12, 20, ()),
            ]
        )
    }
}
