// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A pair Hidden Markov Model for calculating the probability that two sequences are related to
//! each other. Depending on the used parameters, this can, e.g., be used to calculate the
//! probability that a certain sequencing read comes from a given position in a reference genome.

use std::cmp;
use std::mem;
use std::usize;

pub use crate::stats::pairhmm::{
    EmissionParameters, GapParameters, StartEndGapParameters, XYEmission,
};
use crate::stats::LogProb;

/// Fast approximation of sum over the three given proabilities. If the largest is sufficiently
/// large compared to the others, we just return that instead of computing the full (expensive)
/// sum.
#[inline]
fn ln_sum3_exp_approx(mut p0: LogProb, mut p1: LogProb, mut p2: LogProb) -> LogProb {
    if p1 < p2 {
        mem::swap(&mut p1, &mut p2);
    }
    if p1 > p0 {
        mem::swap(&mut p1, &mut p0);
    }
    if *(p0 - p1) > 10.0 {
        // if p0 is strong enough compared to second, just return the maximum
        p0
    } else {
        // calculate accurate sum (slower)
        LogProb::ln_sum_exp(&[p0, p1, p2])
    }
}

/// A pair Hidden Markov Model for comparing sequences x and y as described by
/// Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
/// Current Topics in Genome Analysis 2008. http://doi.org/10.1017/CBO9780511790492.
#[derive(Debug, Clone)]
pub struct PairHMM {
    fm: [Vec<LogProb>; 2],
    fx: [Vec<LogProb>; 2],
    fy: [Vec<LogProb>; 2],
    min_edit_dist: [Vec<usize>; 2],
    prob_cols: Vec<LogProb>,
}

impl Default for PairHMM {
    fn default() -> Self {
        PairHMM {
            fm: [Vec::new(), Vec::new()],
            fx: [Vec::new(), Vec::new()],
            fy: [Vec::new(), Vec::new()],
            min_edit_dist: [Vec::new(), Vec::new()],
            prob_cols: Vec::new(),
        }
    }
}

impl PairHMM {
    pub fn new() -> Self {
        Default::default()
    }

    /// Calculate the probability of sequence x being related to y via any alignment.
    ///
    /// # Arguments
    ///
    /// * `gap_params` - parameters for opening or extending gaps
    /// * `emission_params` - parameters for emission
    /// * `max_edit_dist` - maximum edit distance to consider; if not `None`, perform banded alignment
    pub fn prob_related<G, E>(
        &mut self,
        gap_params: &G,
        emission_params: &E,
        max_edit_dist: Option<usize>,
    ) -> LogProb
    where
        G: GapParameters + StartEndGapParameters,
        E: EmissionParameters,
    {
        for k in 0..2 {
            self.fm[k].clear();
            self.fx[k].clear();
            self.fy[k].clear();
            self.min_edit_dist[k].clear();
            self.prob_cols.clear();

            self.fm[k].resize(emission_params.len_y() + 1, LogProb::ln_zero());
            self.fx[k].resize(emission_params.len_y() + 1, LogProb::ln_zero());
            self.fy[k].resize(emission_params.len_y() + 1, LogProb::ln_zero());
            self.min_edit_dist[k].resize(emission_params.len_y() + 1, usize::MAX);

            if gap_params.free_end_gap_x() {
                let c = (emission_params.len_x() * 3).saturating_sub(self.prob_cols.capacity());
                self.prob_cols.reserve_exact(c);
            }
        }

        // cache probs
        let prob_no_gap = gap_params
            .prob_gap_x()
            .ln_add_exp(gap_params.prob_gap_y())
            .ln_one_minus_exp();
        let prob_no_gap_x_extend = gap_params.prob_gap_x_extend().ln_one_minus_exp();
        let prob_no_gap_y_extend = gap_params.prob_gap_y_extend().ln_one_minus_exp();
        let prob_gap_x = gap_params.prob_gap_x();
        let prob_gap_y = gap_params.prob_gap_y();
        let prob_gap_x_extend = gap_params.prob_gap_x_extend();
        let prob_gap_y_extend = gap_params.prob_gap_y_extend();
        let do_gap_y_extend = prob_gap_y_extend != LogProb::ln_zero();
        let do_gap_x_extend = prob_gap_x_extend != LogProb::ln_zero();

        let mut prev = 0;
        let mut curr = 1;
        self.fm[prev][0] = LogProb::ln_one();

        // iterate over x
        for i in 0..emission_params.len_x() {
            // allow alignment to start from offset in x (if prob_start_gap_x is set accordingly)
            self.fm[prev][0] = self.fm[prev][0].ln_add_exp(gap_params.prob_start_gap_x(i));
            if gap_params.free_start_gap_x() {
                self.min_edit_dist[prev][0] = 0;
            }

            let prob_emit_x = emission_params.prob_emit_x(i);

            // TODO: in the case of no gap extensions, we can reduce the number of columns of y that need to be looked at (by cone).
            let (j_min, j_max) = (0, emission_params.len_y());

            // iterate over y
            for j in j_min..j_max {
                let j_ = j + 1;
                let j_minus_one = j_ - 1;

                let min_edit_dist_topleft = self.min_edit_dist[prev][j_minus_one];
                let min_edit_dist_top = self.min_edit_dist[curr][j_minus_one];
                let min_edit_dist_left = self.min_edit_dist[prev][j_];

                if let Some(max_edit_dist) = max_edit_dist {
                    if cmp::min(
                        min_edit_dist_topleft,
                        cmp::min(min_edit_dist_top, min_edit_dist_left),
                    ) > max_edit_dist
                    {
                        // skip this cell if best edit dist is already larger than given maximum
                        continue;
                    }
                }

                let (prob_match_mismatch, prob_gap_x, prob_gap_y, min_edit_dist) = {
                    let fm_curr = &self.fm[curr];
                    let fm_prev = &self.fm[prev];
                    let fx_prev = &self.fx[prev];
                    let fy_curr = &self.fy[curr];
                    let fy_prev = &self.fy[prev];

                    // match or mismatch
                    let emit_xy = emission_params.prob_emit_xy(i, j);
                    let prob_match_mismatch = emit_xy.prob()
                        + ln_sum3_exp_approx(
                            prob_no_gap + fm_prev[j_minus_one],
                            // coming from state X
                            prob_no_gap_x_extend + fx_prev[j_minus_one],
                            // coming from state Y
                            prob_no_gap_y_extend + fy_prev[j_minus_one],
                        );

                    // gap in y
                    let mut prob_gap_y = prob_emit_x
                        + (
                            // open gap
                            prob_gap_y + fm_prev[j_]
                        );
                    if do_gap_y_extend {
                        prob_gap_y = prob_gap_y.ln_add_exp(
                            // extend gap
                            prob_gap_y_extend + fx_prev[j_],
                        );
                    }

                    // gap in x
                    let mut prob_gap_x = emission_params.prob_emit_y(j)
                        + (
                            // open gap
                            prob_gap_x + fm_curr[j_minus_one]
                        );
                    if do_gap_x_extend {
                        prob_gap_x = prob_gap_x.ln_add_exp(
                            // extend gap
                            prob_gap_x_extend + fy_curr[j_minus_one],
                        );
                    }

                    // calculate minimal number of mismatches
                    let min_edit_dist = if max_edit_dist.is_some() {
                        cmp::min(
                            if emit_xy.is_match() {
                                // a match, so nothing changes
                                min_edit_dist_topleft
                            } else {
                                // one new mismatch
                                min_edit_dist_topleft.saturating_add(1)
                            },
                            cmp::min(
                                // gap in y (no new mismatch)
                                min_edit_dist_left.saturating_add(1),
                                // gap in x (no new mismatch)
                                min_edit_dist_top.saturating_add(1),
                            ),
                        )
                    } else {
                        0
                    };

                    (prob_match_mismatch, prob_gap_x, prob_gap_y, min_edit_dist)
                };

                self.fm[curr][j_] = prob_match_mismatch;
                self.fx[curr][j_] = prob_gap_y;
                self.fy[curr][j_] = prob_gap_x;
                if max_edit_dist.is_some() {
                    self.min_edit_dist[curr][j_] = min_edit_dist;
                }
            }

            if gap_params.free_end_gap_x() {
                // Cache column probabilities or simply record the last probability.
                // We can put all of them in one array since we simply have to sum in the end.
                // This is also good for numerical stability.
                self.prob_cols.push(self.fm[curr].last().unwrap().clone());
                self.prob_cols.push(self.fx[curr].last().unwrap().clone());
                // TODO check removing this (we don't want open gaps in x):
                self.prob_cols.push(self.fy[curr].last().unwrap().clone());
            }

            // next column
            mem::swap(&mut curr, &mut prev);
            // reset next column to zeros
            for v in &mut self.fm[curr] {
                *v = LogProb::ln_zero();
            }
        }

        let p = if gap_params.free_end_gap_x() {
            LogProb::ln_sum_exp(&self.prob_cols)
        } else {
            LogProb::ln_sum_exp(&[
                *self.fm[prev].last().unwrap(),
                *self.fx[prev].last().unwrap(),
                *self.fy[prev].last().unwrap(),
            ])
        };
        // take the minimum with 1.0, because sum of paths can exceed probability 1.0
        // especially in case of repeats
        assert!(!p.is_nan());
        if p > LogProb::ln_one() {
            LogProb::ln_one()
        } else {
            p
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::stats::{LogProb, Prob};

    use super::*;

    // Single base insertion and deletion rates for R1 according to Schirmer et al.
    // BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
    static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
    static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
    static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

    fn prob_emit_x_or_y() -> LogProb {
        LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
    }

    struct TestEmissionParams {
        x: &'static [u8],
        y: &'static [u8],
    }

    impl EmissionParameters for TestEmissionParams {
        fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
            if self.x[i] == self.y[j] {
                XYEmission::Match(LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST))
            } else {
                XYEmission::Mismatch(LogProb::from(PROB_ILLUMINA_SUBST / Prob(3.0)))
            }
        }

        fn prob_emit_x(&self, _: usize) -> LogProb {
            prob_emit_x_or_y()
        }

        fn prob_emit_y(&self, _: usize) -> LogProb {
            prob_emit_x_or_y()
        }

        fn len_x(&self) -> usize {
            self.x.len()
        }

        fn len_y(&self) -> usize {
            self.y.len()
        }
    }

    struct TestSingleGapParams;

    impl GapParameters for TestSingleGapParams {
        fn prob_gap_x(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_INS)
        }

        fn prob_gap_y(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_DEL)
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }
    }

    impl StartEndGapParameters for TestSingleGapParams {
        fn free_start_gap_x(&self) -> bool {
            false
        }

        fn free_end_gap_x(&self) -> bool {
            false
        }
    }

    pub struct SemiglobalGapParams;

    impl GapParameters for SemiglobalGapParams {
        fn prob_gap_x(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_INS)
        }

        fn prob_gap_y(&self) -> LogProb {
            LogProb::from(PROB_ILLUMINA_DEL)
        }

        fn prob_gap_x_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }

        fn prob_gap_y_extend(&self) -> LogProb {
            LogProb::ln_zero()
        }
    }

    impl StartEndGapParameters for SemiglobalGapParams {
        fn free_start_gap_x(&self) -> bool {
            true
        }

        fn free_end_gap_x(&self) -> bool {
            true
        }
    }

    const EMIT_MATCH: LogProb = LogProb(-0.0021022080918701985);
    const EMIT_GAP_X: LogProb = LogProb(-0.0021022080918701985);
    const EMIT_GAP_Y: LogProb = LogProb(-0.0021022080918701985);
    const T_MATCH: LogProb = LogProb(-7.900031205113962e-06);
    const T_GAP_X: LogProb = LogProb(-12.785891140783116);
    const T_GAP_Y: LogProb = LogProb(-12.186270018233994);

    #[test]
    fn impossible_global_alignment() {
        let x = b"AAA";
        let y = b"A";
        let emission_params = TestEmissionParams { x, y };

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&TestSingleGapParams, &emission_params, None);
        assert_eq!(p, LogProb::ln_zero());
    }

    #[test]
    fn test_interleave_gaps_y() {
        let x = b"ACGTACGTACGT";
        let y = b"AGAGAG";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestSingleGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let n_matches = 6.;
        let n_insertions = 6.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_Y * n_insertions
                + *T_GAP_Y * n_insertions
                + (1. - *PROB_ILLUMINA_DEL).ln() * n_insertions,
        );

        let p_max = LogProb(*T_GAP_Y * n_insertions);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_interleave_gaps_x() {
        let x = b"AGAGAG";
        let y = b"ACGTACGTACGT";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestSingleGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let n_matches = 6.;
        let n_insertions = 6.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_X * n_insertions
                + *T_GAP_X * n_insertions
                + (1. - *PROB_ILLUMINA_INS).ln() * n_insertions,
        );

        let p_max = LogProb(*T_GAP_X * n_insertions);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_same() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestSingleGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);
        let n = x.len() as f64;
        let p_most_likely_path = LogProb(*EMIT_MATCH * n + *T_MATCH * (n - 1.));
        let p_max = LogProb(*EMIT_MATCH * n);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.001);
        assert_relative_eq!(*p, *p_max, epsilon = 0.001);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_x() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCTGATCGATCT";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestSingleGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let n_matches = 17.;
        let n_insertions = 2.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH * (n_matches - n_insertions)
                + *EMIT_GAP_X * n_insertions
                + *T_GAP_X * n_insertions
                + (1. - *PROB_ILLUMINA_INS).ln(),
        );

        let p_max = LogProb(*T_GAP_X * 2.);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_gap_y() {
        let x = b"AGCTCGATCTGATCGATCT";
        let y = b"AGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestSingleGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let n_matches = 17.;
        let n_deletions = 2.;

        let p_most_likely_path = LogProb(
            *EMIT_MATCH * n_matches
                + *T_MATCH * (n_matches - n_deletions)
                + *EMIT_GAP_Y * n_deletions
                + *T_GAP_Y * n_deletions
                + (1. - *PROB_ILLUMINA_DEL).ln(),
        );

        let p_max = LogProb(*T_GAP_Y * 2.);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 0.01);
        assert_relative_eq!(*p, *p_max, epsilon = 0.1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_mismatch() {
        let x = b"AGCTCGAGCGATCGATC";
        let y = b"TGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestSingleGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let n = x.len() as f64;
        let p_most_likely_path = LogProb(
            *EMIT_MATCH * (n - 2.) + *T_MATCH * (n - 1.) + (*PROB_ILLUMINA_SUBST / 3.).ln() * 2.,
        );
        let p_max = LogProb((*PROB_ILLUMINA_SUBST / 3.).ln() * 2.);
        assert!(*p <= 0.0);
        assert_relative_eq!(*p_most_likely_path, *p, epsilon = 1e-4);
        assert_relative_eq!(*p, *p_max, epsilon = 1e-1);
        assert!(*p <= *p_max);
    }

    #[test]
    fn test_banded() {
        let x = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC\
ATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTAT\
CTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGT";

        let y = b"GGGTATGCACGCGATAGCATTGCGAGATGCTGGAGCTGGAGCACCCTATGTCGC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = SemiglobalGapParams;

        let mut pair_hmm = PairHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let p_banded = pair_hmm.prob_related(&gap_params, &emission_params, Some(2));

        assert_relative_eq!(*p, *p_banded, epsilon = 1e-7);
    }
}
