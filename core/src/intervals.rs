use coitrees::{IntWithMax, Interval, IntervalNode};
use itertools::Itertools;
use std::collections::VecDeque;

/// Merge intervals in a [`COITree`]. Includes book-ended intervals.
///
/// # Arguments
/// * `intervals`: Intervals to merge. Elements are cloned.
/// * `dst`: Distance to merge over.
/// * `fn_cmp`: Function to enforce additional check before merging.
/// * `fn_reducer`: Function to reduce metadata.
/// * `fn_finalizer`: Function to apply some final operation on intervals.
///
/// # Returns
/// * Merged overlapping intervals.
pub fn merge_intervals<I, T>(
    intervals: I,
    dst: i32,
    fn_cmp: impl Fn(&Interval<T>, &Interval<T>) -> bool,
    fn_reducer: impl Fn(&Interval<T>, &Interval<T>) -> T,
    fn_finalizer: impl Fn(Interval<T>) -> Interval<T>,
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
        let dst_between = itv_2.first - itv_1.last;
        let added_check = fn_cmp(&itv_1, &itv_2);
        if (dst_between <= dst) & added_check {
            let new_data = fn_reducer(&itv_1, &itv_2);
            let merged_interval = Interval::new(itv_1.first, itv_2.last, new_data);
            intervals.push_front(merged_interval);
        } else {
            merged.push(itv_1);
            intervals.push_front(itv_2);
        }
    }
    // Apply finalizer function
    merged.into_iter().map(fn_finalizer).collect_vec()
}

pub(crate) fn trim_coords<T: Clone, I: IntWithMax>(
    st: &mut i32,
    end: &mut i32,
    ovl_itv: &IntervalNode<T, I>,
    status: &mut String,
    split_itv_1: &mut Option<(i32, i32)>,
    split_itv_2: &mut Option<(i32, i32)>,
) {
    // Cases:
    //    |---|
    // * |---|
    // Current interval overlaps right
    // Trim interval end to ignore interval start.
    //
    //   |---|
    // *  |---|
    // Current interval overlaps left
    // Trim interval start to ignore interval end.
    //
    //   |---|
    // * |---|
    // Current interval contained in ignored interval.
    // Set status of interval to good.
    //
    //    |-|
    // * |---|
    // Ignore interval contained in current interval.
    // Create two intervals.
    // * Interval start to ignore interval start w/status
    // * Ignore interval end to interval end w/status
    if *end > ovl_itv.first && *end < ovl_itv.last {
        // Trim to overlap.
        *end = ovl_itv.first
    } else if *st < ovl_itv.last && *st > ovl_itv.first {
        *st = ovl_itv.last
    } else if *st >= ovl_itv.first && *end <= ovl_itv.last {
        // Remove splits since entire interval overlapped.
        split_itv_1.take();
        split_itv_2.take();

        status.clear();
        status.push_str("good");
    } else if ovl_itv.first > *st && ovl_itv.last < *end {
        // Take min or max coordinate trimming to bounds of contained, overlapping interval.
        if let Some(old_split_itv_1) = split_itv_1 {
            *old_split_itv_1 = (*st, std::cmp::min(ovl_itv.first, old_split_itv_1.1))
        } else {
            *split_itv_1 = Some((*st, ovl_itv.first))
        }
        if let Some(old_split_itv_2) = split_itv_2 {
            *old_split_itv_2 = (std::cmp::max(ovl_itv.last, old_split_itv_2.0), *end)
        } else {
            *split_itv_2 = Some((ovl_itv.last, *end))
        }

        // Ignore center.
        *st = ovl_itv.first;
        *end = ovl_itv.last;
        status.clear();
        status.push_str("good");
    }
}

#[cfg(test)]
mod tests {
    use super::{merge_intervals, trim_coords};
    use coitrees::{Interval, IntervalNode};
    use std::fmt::Debug;

    const ST: i32 = 4;
    const END: i32 = 8;

    fn reduce_to_a<'a>(a: &Interval<usize>, _b: &Interval<usize>) -> usize {
        a.metadata
    }

    fn noop(a: Interval<usize>) -> Interval<usize> {
        a
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
        let merged_itvs = merge_intervals(
            itvs.clone().into_iter(),
            0,
            no_added_check,
            reduce_to_a,
            noop,
        );
        assert_itvs_equal(&itvs, &merged_itvs);
    }

    #[test]
    fn test_merge_intervals_single() {
        let itvs = vec![
            Interval::new(1, 3, 1),
            Interval::new(3, 5, 2),
            Interval::new(6, 9, 3),
        ];
        let merged_itvs = merge_intervals(itvs.into_iter(), 0, no_added_check, reduce_to_a, noop);
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
        let merged_itvs = merge_intervals(itvs.into_iter(), 0, no_added_check, reduce_to_a, noop);
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
            0,
            |a, b| (a.metadata % 2 == 0) & (b.metadata % 2 == 0),
            reduce_to_a,
            noop,
        );
        assert_itvs_equal(&merged_itvs, &exp_itvs);
    }

    #[test]
    fn test_trim_ovl_itv_none() {
        // 12345678
        // |-|
        //    |---|
        let (mut st, mut end) = (ST, END);
        let mut status = String::from("bad");
        let ovl_itv: IntervalNode<String, usize> = IntervalNode::new(1, 3, String::from("left"));
        let (mut split_itv_1, mut split_itv_2) = (None, None);

        trim_coords(
            &mut st,
            &mut end,
            &ovl_itv,
            &mut status,
            &mut split_itv_1,
            &mut split_itv_2,
        );
        assert!(
            st == 4
                && end == 8
                && &status == "bad"
                && split_itv_1.is_none()
                && split_itv_2.is_none()
        );
    }

    #[test]
    fn test_trim_ovl_itv_left() {
        // 12345678
        // |---|
        //    |---|
        let (mut st, mut end) = (ST, END);
        let mut status = String::from("bad");
        let ovl_itv: IntervalNode<String, usize> = IntervalNode::new(1, 5, String::from("left"));
        let (mut split_itv_1, mut split_itv_2) = (None, None);

        trim_coords(
            &mut st,
            &mut end,
            &ovl_itv,
            &mut status,
            &mut split_itv_1,
            &mut split_itv_2,
        );
        assert!(
            st == 5
                && end == 8
                && &status == "bad"
                && split_itv_1.is_none()
                && split_itv_2.is_none()
        );
    }

    #[test]
    fn test_trim_ovl_itv_right() {
        // 123456789X
        //      |---|
        //    |---|
        let (mut st, mut end) = (ST, END);
        let mut status = String::from("bad");
        let ovl_itv: IntervalNode<String, usize> = IntervalNode::new(6, 10, String::from("right"));
        let (mut split_itv_1, mut split_itv_2) = (None, None);

        trim_coords(
            &mut st,
            &mut end,
            &ovl_itv,
            &mut status,
            &mut split_itv_1,
            &mut split_itv_2,
        );
        assert!(
            st == 4
                && end == 6
                && &status == "bad"
                && split_itv_1.is_none()
                && split_itv_2.is_none()
        );
    }

    #[test]
    fn test_trim_ovl_itv_center() {
        // 123456789X
        //     |-|
        //    |---|
        let (mut st, mut end) = (ST, END);
        let mut status = String::from("bad");
        let ovl_itv: IntervalNode<String, usize> = IntervalNode::new(5, 7, String::from("center"));
        let (mut split_itv_1, mut split_itv_2) = (None, None);

        trim_coords(
            &mut st,
            &mut end,
            &ovl_itv,
            &mut status,
            &mut split_itv_1,
            &mut split_itv_2,
        );
        assert!(
            st == 5
                && end == 7
                && &status == "good"
                && split_itv_1 == Some((4, 5))
                && split_itv_2 == Some((7, 8))
        );
        // dbg!(st, end, status, split_itv_1, split_itv_2);
    }

    #[test]
    fn test_trim_ovl_itv_contains() {
        // 123456789X
        //    |---|
        //    |---|
        let (mut st, mut end) = (ST, END);
        let mut status = String::from("bad");
        let ovl_itv: IntervalNode<String, usize> =
            IntervalNode::new(ST, END, String::from("contains"));
        let (mut split_itv_1, mut split_itv_2) = (None, None);

        trim_coords(
            &mut st,
            &mut end,
            &ovl_itv,
            &mut status,
            &mut split_itv_1,
            &mut split_itv_2,
        );
        assert!(
            st == ST
                && end == END
                && &status == "good"
                && split_itv_1.is_none()
                && split_itv_2.is_none()
        );

        // Also check that removes any pre-existing split intervals.
        let (mut split_itv_1, mut split_itv_2) = (Some((0, 1)), Some((1, 2)));
        trim_coords(
            &mut st,
            &mut end,
            &ovl_itv,
            &mut status,
            &mut split_itv_1,
            &mut split_itv_2,
        );
        assert!(
            st == ST
                && end == END
                && &status == "good"
                && split_itv_1.is_none()
                && split_itv_2.is_none()
        );
    }
}
