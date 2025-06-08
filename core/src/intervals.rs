use coitrees::Interval;
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

pub fn overlap_length(a_first: i32, a_last: i32, b_first: i32, b_last: i32) -> i32 {
    if a_first >= b_first && a_last >= b_last {
        // a  |---|
        // b |---|
        b_last - a_first
    } else if a_first <= b_first && a_last <= b_last {
        // a |---|
        // b  |---|
        a_last - b_first
    } else if a_first <= b_first && a_last >= b_last {
        // a |-----|
        // b  |---|
        b_last - b_first
    } else if a_first >= b_first && a_last <= b_last {
        // a  |-|
        // b |---|
        a_last - a_first
    } else {
        0
    }
}

/// Subtract interval by a list of non-overlapping intervals.
/// * See [`merge_intervals`].
pub fn subtract_intervals<T: Clone>(
    itv: Interval<T>,
    other: impl Iterator<Item = Interval<T>>,
) -> Vec<Interval<T>> {
    let mut split_intervals = Vec::new();
    let mut st = itv.first;
    let mut last = itv.last;
    for ovl_itv in other.into_iter().sorted_by(|a, b| a.first.cmp(&b.first)) {
        if last >= ovl_itv.first && last <= ovl_itv.last {
            //    |---|
            // * |---|
            last = ovl_itv.first;
        } else if st <= ovl_itv.last && st >= ovl_itv.first {
            //   |---|
            // *  |---|
            st = ovl_itv.last;
        } else if st >= ovl_itv.first && last <= ovl_itv.last {
            //   |---|
            // * |---|
            break;
        } else if ovl_itv.first > st && ovl_itv.last < last {
            //    |-|
            // * |---|
            split_intervals.push(Interval::new(st, ovl_itv.first, itv.metadata.clone()));
            st = ovl_itv.last;
        }
    }
    // Add remainder.
    if st != last {
        split_intervals.push(Interval::new(st, last, itv.metadata.clone()));
    }
    split_intervals
}

#[cfg(test)]
mod tests {
    use super::{merge_intervals, subtract_intervals};
    use coitrees::Interval;
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
            1,
            |a, b| (a.metadata % 2 == 0) & (b.metadata % 2 == 0),
            reduce_to_a,
            noop,
        );
        assert_itvs_equal(&merged_itvs, &exp_itvs);
    }

    #[test]
    fn test_subtract_no_ovl() {
        // 12345678
        // |-|
        //    |---|
        let itv = Interval::new(ST, END, "");
        let itvs = [Interval::new(1, 3, "")];
        let res = subtract_intervals(itv, itvs.into_iter());
        assert_itvs_equal(&[itv], &res);
    }

    #[test]
    fn test_subtract_left_ovl() {
        // 12345678
        // |---|
        //    |---|
        let itv = Interval::new(ST, END, "");
        let itvs = [Interval::new(1, 5, "")];
        let res = subtract_intervals(itv, itvs.into_iter());
        assert_itvs_equal(&[Interval::new(5, 8, "")], &res);
    }

    #[test]
    fn test_subtract_right_ovl() {
        // 123456789X
        //      |---|
        //    |---|
        let itv = Interval::new(ST, END, "");
        let itvs = [Interval::new(6, 10, "")];
        let res = subtract_intervals(itv, itvs.into_iter());
        assert_itvs_equal(&[Interval::new(4, 6, "")], &res);
    }

    #[test]
    fn test_subtract_center_ovl() {
        // 123456789X
        //     |-|
        //    |---|
        let itv = Interval::new(ST, END, "");
        let itvs = [Interval::new(5, 7, "")];
        let res = subtract_intervals(itv, itvs.into_iter());
        assert_itvs_equal(&[Interval::new(4, 5, ""), Interval::new(7, 8, "")], &res);
    }

    #[test]
    fn test_subtract_contained_ovl() {
        // 123456789X
        //    |---|
        //    |---|
        let itv = Interval::new(ST, END, "");
        let itvs = [Interval::new(ST, END, "")];
        let res = subtract_intervals(itv, itvs.into_iter());
        assert!(res.is_empty());
    }

    #[test]
    fn test_subtract_multiple_ovl() {
        // 123456789X
        //  |-||-| ||
        // |--------|
        let itv = Interval::new(1, 10, "");
        let itvs = [
            Interval::new(2, 4, ""),
            Interval::new(5, 7, ""),
            Interval::new(9, 10, ""),
        ];
        let res = subtract_intervals(itv, itvs.into_iter());
        assert_itvs_equal(
            &[
                Interval::new(1, 2, ""),
                Interval::new(4, 5, ""),
                Interval::new(7, 9, ""),
            ],
            &res,
        );
    }
}
