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

use crate::stats::error_profiles::{
    Base, EmissionParameters, ErrorProfile, GapParameters, HomopolymerParameters,
    StartEndGapParameters,
};
use crate::stats::{LogProb, Prob};
use itertools::Itertools;
use num_traits::Zero;
use smallvec::SmallVec;
use num_traits::real::Real;

fn cache_probs<F>(f: F) -> SmallVec<[LogProb; 4]>
where
    F: Fn(Base) -> LogProb,
{
    Base::values().iter().map(|&base| f(base)).collect()
}

/// A pair Hidden Markov Model for comparing sequences x and y as described by
/// Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). Biological Sequence Analysis.
/// Current Topics in Genome Analysis 2008. http://doi.org/10.1017/CBO9780511790492.
#[derive(Debug, Clone)]
pub struct PairHHMM {
    fmm: [[Vec<LogProb>; 4]; 2],
    fgx: [Vec<LogProb>; 2],
    fgy: [Vec<LogProb>; 2],
    fhx: [[Vec<LogProb>; 4]; 2],
    fhy: [[Vec<LogProb>; 4]; 2],
    min_edit_dist: [Vec<usize>; 2],
    prob_cols: Vec<LogProb>,
}

impl Default for PairHHMM {
    fn default() -> Self {
        PairHHMM {
            fmm: [
                [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
                [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            ],
            fgx: [Vec::new(), Vec::new()],
            fgy: [Vec::new(), Vec::new()],
            fhx: [
                [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
                [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            ],
            fhy: [
                [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
                [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            ],
            min_edit_dist: [Vec::new(), Vec::new()],
            prob_cols: Vec::new(),
        }
    }
}

const NUM_BASES: usize = 4;

impl PairHHMM {
    pub fn new() -> Self {
        Default::default()
    }

    /// Calculate the probability of sequence x being related to y via any alignment.
    ///
    /// # Arguments
    ///
    /// * `error_profile` - bundled parameters for emissions, gaps and homopolymer runs.
    /// * `start_end_params` - e.g. global or semi-global alignment, see `Locality`.
    /// * `max_edit_dist` - maximum edit distance to consider; if not `None`, perform banded alignment.
    pub fn prob_related<E, G>(
        &mut self,
        error_profile: &E,
        start_end_params: &G,
        max_edit_dist: Option<usize>,
    ) -> LogProb
    where
        G: StartEndGapParameters,
        E: ErrorProfile,
    {
        let emission_params = error_profile.emission_parameters();
        let gap_params = error_profile.gap_parameters();
        let homopolymer_params = error_profile.homopolymer_parameters();
        for k in 0..2 {
            self.prob_cols.clear();
            let l = emission_params.len_y() + 1;
            self.fgx[k].reset_and_resize(l, LogProb::ln_zero());
            self.fgy[k].reset_and_resize(l, LogProb::ln_zero());
            for b in 0..NUM_BASES {
                self.fhx[k][b].reset_and_resize(l, LogProb::ln_zero());
                self.fhy[k][b].reset_and_resize(l, LogProb::ln_zero());
                self.fmm[k][b].reset_and_resize(l, LogProb::ln_zero());
            }
            self.min_edit_dist[k].reset_and_resize(l, usize::MAX);

            if start_end_params.free_end_gap_x() {
                let c = (emission_params.len_x() * 3).saturating_sub(self.prob_cols.capacity());
                self.prob_cols.reserve_exact(c);
            }
        }

        // cache probs
        let prob_no_gap_x_extend = gap_params.prob_gap_x_extend().ln_one_minus_exp();
        let prob_no_gap_y_extend = gap_params.prob_gap_y_extend().ln_one_minus_exp();
        let prob_open_gap_x = gap_params.prob_gap_x();
        let prob_open_gap_y = gap_params.prob_gap_y();
        let prob_extend_gap_x = gap_params.prob_gap_x_extend();
        let prob_extend_gap_y = gap_params.prob_gap_y_extend();
        let extend_gap_in_y = prob_extend_gap_y != LogProb::ln_zero();
        let extend_gap_in_x = prob_extend_gap_x != LogProb::ln_zero();

        let no_gap_no_hop_probs = cache_probs(|base| {
            homopolymer_params
                .prob_hop_x(base)
                .ln_add_exp(homopolymer_params.prob_hop_y(base))
                .ln_add_exp(gap_params.prob_gap_x())
                .ln_add_exp(gap_params.prob_gap_y())
                .ln_one_minus_exp()
        });

        let no_gap_no_hop_x_extend_probs = cache_probs(|base| {
            homopolymer_params
                .prob_hop_x_extend(base)
                .ln_add_exp(gap_params.prob_gap_x())
                .ln_add_exp(gap_params.prob_gap_y())
                .ln_one_minus_exp()
        });

        let no_gap_no_hop_y_extend_probs = cache_probs(|base| {
            homopolymer_params
                .prob_hop_y_extend(base)
                .ln_add_exp(gap_params.prob_gap_x())
                .ln_add_exp(gap_params.prob_gap_y())
                .ln_one_minus_exp()
        });

        let start_hop_x_probs = cache_probs(|base| homopolymer_params.prob_hop_x(base));
        let start_hop_y_probs = cache_probs(|base| homopolymer_params.prob_hop_y(base));
        let extend_hop_x_probs = cache_probs(|base| homopolymer_params.prob_hop_x_extend(base));
        let extend_hop_y_probs = cache_probs(|base| homopolymer_params.prob_hop_y_extend(base));

        let mut prev = 0;
        let mut curr = 1;
        for b in 0..NUM_BASES {
            self.fmm[prev][b][0] = LogProb::ln_one();
        }
        println!("{}\t{}", -1, self.fmm[curr][0].iter().map(|p|Prob::from(*p).0).map(|v| format!("{:>.8}", v)).join("\t"));
        // iterate over x
        for i in 0..emission_params.len_x() {
            // allow alignment to start from offset in x (if prob_start_gap_x is set accordingly)
            for b in 0..NUM_BASES {
                self.fmm[prev][b][0] =
                    self.fmm[prev][b][0].ln_add_exp(start_end_params.prob_start_gap_x(i));
            }

            if start_end_params.free_start_gap_x() {
                self.min_edit_dist[prev][0] = 0;
            }

            let emit_x_and_gap_probs =
                cache_probs(|base| emission_params.prob_emit_x_and_gap(base, i));

            let j_min = if extend_gap_in_x {
                0
            } else {
                // The final upper triangle in the dynamic programming matrix
                // will be only zeros if gap extension is not allowed on x.
                // Hence we can omit it.
                emission_params
                    .len_y()
                    .saturating_sub(emission_params.len_x() - i + 1)
            };

            let j_max = if extend_gap_in_x {
                emission_params.len_y()
            } else {
                // The initial lower triangle in the dynamic programming matrix
                // will be only zeros if gap extension is not allowed on x.
                // Hence we can omit it.
                cmp::min(i + 1, emission_params.len_y())
            };

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

                let (
                    match_mismatch_probs,
                    gap_in_x_prob,
                    gap_in_y_prob,
                    hop_in_x_probs,
                    hop_in_y_probs,
                    min_edit_dist,
                ) = {
                    // TODO can be shortened using a macro which makes use of `concat_idents!`
                    let fmm_curr = &self.fmm[curr];
                    let fmm_prev = &self.fmm[prev];
                    let fgx_prev = &self.fgx[prev];
                    let fgy_curr = &self.fgy[curr];
                    let fgy_prev = &self.fgy[prev];
                    let fhx_curr = &self.fhx[curr];
                    let fhx_prev = &self.fhx[prev];
                    let fhy_curr = &self.fhy[curr];
                    let fhy_prev = &self.fhy[prev];

                    // match or mismatch
                    let match_mismatch_probs: SmallVec<[LogProb; NUM_BASES]> = Base::values()
                        .iter()
                        .map(|&b1| {
                            let b1_ = b1 as usize;
                            emission_params.prob_emit_x_and_y(b1, i, j).prob()
                                + LogProb::ln_sum_exp(
                                    &(0..NUM_BASES)
                                        .flat_map(|b2| {
                                            vec![
                                                // coming from any of the match states
                                                no_gap_no_hop_probs[b2] - LogProb::from(Prob(4.0))
                                                    + fmm_prev[b1_][j_minus_one],
                                                // coming from any of the hop x states (!= b1)
                                                if b2 != b1_ {
                                                    no_gap_no_hop_x_extend_probs[b2]
                                                        - LogProb::from(Prob(3.0))
                                                        + fhx_prev[b2][j_minus_one]
                                                } else {
                                                    LogProb::zero()
                                                },
                                                // coming from any of the hop y states (!= b1)
                                                if b2 != b1_ {
                                                    no_gap_no_hop_y_extend_probs[b2]
                                                        - LogProb::from(Prob(3.0))
                                                        + fhy_prev[b2][j_minus_one]
                                                } else {
                                                    LogProb::zero()
                                                },
                                            ]
                                        })
                                        .chain(vec![
                                            // coming from one of the gap states
                                            prob_no_gap_x_extend - LogProb::from(Prob(4.0))
                                                + fgx_prev[j_minus_one],
                                            prob_no_gap_y_extend - LogProb::from(Prob(4.0))
                                                + fgy_prev[j_minus_one],
                                        ])
                                        .collect_vec(),
                                )
                        })
                        .collect();

                    // gap in y
                    let mut prob_gap_in_y = emit_x_and_gap_probs[0]
                        + LogProb::ln_sum_exp(
                            &(0..NUM_BASES)
                                .flat_map(|b| {
                                    vec![
                                        // open gap from one of the match states
                                        prob_open_gap_y + fmm_prev[b][j_],
                                        // open gap from one of the hop x states
                                        prob_open_gap_y + fhx_prev[b][j_],
                                        // open gap from one of the hop y states
                                        prob_open_gap_y + fhy_prev[b][j_],
                                    ]
                                })
                                .collect_vec(),
                        );
                    if extend_gap_in_y {
                        prob_gap_in_y = prob_gap_in_y.ln_add_exp(
                            // extend gap
                            prob_extend_gap_y + fgx_prev[j_],
                        );
                    }

                    // gap in x
                    let mut prob_gap_in_x = emission_params.prob_emit_gap_and_y(Base::A, j) // any base will do
                        + LogProb::ln_sum_exp(
                            &(0..NUM_BASES)
                                .flat_map(|b| {
                                    vec![
                                        // open gap from one of the match states
                                        prob_open_gap_x + fmm_curr[b][j_minus_one],
                                        // open gap from one of the hop x states
                                        prob_open_gap_x + fhx_curr[b][j_minus_one],
                                        // open gap from one of the hop y states
                                        prob_open_gap_x + fhy_curr[b][j_minus_one],
                                    ]
                                })
                                .collect_vec(),
                        );
                    if extend_gap_in_x {
                        prob_gap_in_x = prob_gap_in_x.ln_add_exp(
                            // extend gap
                            prob_extend_gap_x + fgy_curr[j_minus_one],
                        );
                    }

                    // hop in y
                    let hop_in_y_probs: SmallVec<[LogProb; NUM_BASES]> = Base::values()
                        .iter()
                        .map(|&base| {
                            let b = base as usize;
                            emission_params.prob_emit_hop_and_y(base, j)
                                + (start_hop_y_probs[b] + fmm_prev[b][j_]).ln_add_exp(
                                    if extend_hop_y_probs[b] != LogProb::zero() {
                                        extend_hop_y_probs[b] + fhy_prev[b][j_]
                                    } else {
                                        LogProb::zero()
                                    },
                                )
                        })
                        .collect();

                    // hop in x
                    let hop_in_x_probs: SmallVec<[LogProb; NUM_BASES]> = Base::values()
                        .iter()
                        .map(|&base| {
                            let b = base as usize;
                            emission_params.prob_emit_x_and_hop(base, i)
                                + (start_hop_x_probs[b] + fmm_curr[b][j_minus_one]).ln_add_exp(
                                    if extend_hop_x_probs[b] != LogProb::zero() {
                                        extend_hop_x_probs[b] + fhx_curr[b][j_minus_one]
                                    } else {
                                        LogProb::zero()
                                    },
                                )
                        })
                        .collect();

                    // calculate minimal number of mismatches
                    let is_match = Base::values()
                        .iter()
                        .any(|&base| emission_params.prob_emit_x_and_y(base, i, j).is_match());
                    let min_edit_dist = if max_edit_dist.is_some() {
                        cmp::min(
                            if is_match {
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

                    (
                        match_mismatch_probs,
                        prob_gap_in_x,
                        prob_gap_in_y,
                        hop_in_x_probs,
                        hop_in_y_probs,
                        min_edit_dist,
                    )
                };
                let m = match_mismatch_probs.iter().cloned().collect_vec();
                let hx = hop_in_x_probs.iter().cloned().collect_vec();
                let hy = hop_in_y_probs.iter().cloned().collect_vec();
                let match_prob = LogProb::ln_sum_exp(&m).exp();
                let hx_prob = LogProb::ln_sum_exp(&hx).exp();
                let hy_prob = LogProb::ln_sum_exp(&hy).exp();
                let gx_prob = gap_in_x_prob.exp();
                let gy_prob = gap_in_y_prob.exp();
                // dbg!(((i, j), match_prob, hx_prob, hy_prob, gx_prob, gy_prob));

                for &base in &Base::values() {
                    let b = base as usize;
                    self.fmm[curr][b][j_] = match_mismatch_probs[b];
                    self.fhx[curr][b][j_] = hop_in_x_probs[b];
                    self.fhy[curr][b][j_] = hop_in_y_probs[b];
                }
                self.fgx[curr][j_] = gap_in_y_prob;
                self.fgy[curr][j_] = gap_in_x_prob;
                if max_edit_dist.is_some() {
                    self.min_edit_dist[curr][j_] = min_edit_dist;
                }
            }

            if start_end_params.free_end_gap_x() {
                // Cache column probabilities or simply record the last probability.
                // We can put all of them in one array since we simply have to sum in the end.
                // This is also good for numerical stability.
                for &base in &Base::values() {
                    self.prob_cols
                        .push(self.fmm[curr][base as usize].last().unwrap().clone());
                    self.prob_cols
                        .push(self.fhx[curr][base as usize].last().unwrap().clone());
                    self.prob_cols
                        .push(self.fhy[curr][base as usize].last().unwrap().clone());
                }
                self.prob_cols.push(self.fgx[curr].last().unwrap().clone());
                // TODO check removing this (we don't want open gaps in x):
                self.prob_cols.push(self.fgy[curr].last().unwrap().clone());
            }
            println!("{}\t{}", i, self.fmm[curr][0].iter().map(|p| Prob::from(*p).0).map(|v| format!("{:>.8}", v)).join("\t"));
            // next column
            mem::swap(&mut curr, &mut prev);
            // reset next column to zeros
            for &base in &Base::values() {
                self.fmm[curr][base as usize].reset(LogProb::ln_zero());
            }
        }

        let p = if start_end_params.free_end_gap_x() {
            LogProb::ln_sum_exp(&self.prob_cols)
        } else {
            LogProb::ln_sum_exp(
                &(0..NUM_BASES)
                    .flat_map(|b| {
                        vec![
                            *self.fmm[prev][b].last().unwrap(),
                            *self.fhx[prev][b].last().unwrap(),
                            *self.fhy[prev][b].last().unwrap(),
                        ]
                    })
                    .chain(vec![
                        *self.fgx[prev].last().unwrap(),
                        *self.fgy[prev].last().unwrap(),
                    ])
                    .collect_vec(),
            )
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

pub trait Reset<T: Copy> {
    fn reset(&mut self, value: T);
}

impl<T: Copy> Reset<T> for [T] {
    fn reset(&mut self, value: T) {
        for v in self {
            *v = value;
        }
    }
}

pub trait ResetResize<T: Copy> {
    fn reset_and_resize(&mut self, new_len: usize, value: T);
}

impl<T: Copy> ResetResize<T> for Vec<T> {
    fn reset_and_resize(&mut self, new_len: usize, value: T) {
        self.reset(value);
        self.resize(new_len, value);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::error_profiles::illumina::IlluminaR1ErrorProfile;
    use crate::stats::error_profiles::illumina::{
        PROB_ILLUMINA_DEL, PROB_ILLUMINA_INS, PROB_ILLUMINA_SUBST,
    };
    use crate::stats::error_profiles::Locality;
    use crate::stats::Prob;

    #[test]
    fn test_same() {
        let x = b"AAAA";
        let y = b"AAAA";

        let error_profile = IlluminaR1ErrorProfile { x, y };
        let start_end_gap_params = Locality::Global;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&error_profile, &start_end_gap_params, None);
        let l = x.len() as i32;
        dbg!(p.exp(), (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(l));
        assert!(*p <= 0.0);
        assert_relative_eq!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(l),
            epsilon = 0.001
        );
    }

    #[test]
    fn test_substitution() {
        let x = b"AAAA";
        let y = b"AAAC";

        let error_profile = IlluminaR1ErrorProfile { x, y };
        let start_end_gap_params = Locality::Global;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&error_profile, &start_end_gap_params, None);
        let l = x.len() as i32;
        dbg!(p.exp(), (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(l - 1) * PROB_ILLUMINA_SUBST.0);
        assert!(*p <= 0.0);
        assert_relative_eq!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(l - 1) * PROB_ILLUMINA_SUBST.0,
            epsilon = 0.001
        );
    }

    #[test]
    fn test_homopolymer() {
        let x = b"AGCCCT";
        let y = b"AGCT";

        let error_profile = IlluminaR1ErrorProfile { x, y };
        let start_end_gap_params = Locality::Global;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&error_profile, &start_end_gap_params, None);

        assert!(*p <= 0.0);
        dbg!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(4) * 0.07 * 0.07
        );
        assert_relative_eq!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(4) * 0.07 * 0.07,
            epsilon = 1e-11
        );
    }

    #[test]
    fn test_insertion() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCTGATCGATCT";

        let error_profile = IlluminaR1ErrorProfile { x, y };
        let start_end_gap_params = Locality::Global;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&error_profile, &start_end_gap_params, None);

        dbg!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(17) * PROB_ILLUMINA_INS.powi(2)
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(17) * PROB_ILLUMINA_INS.powi(2),
            epsilon = 1e-11
        );
    }

    #[test]
    fn test_deletion() {
        let x = b"AGCTCGATCTGATCGATCT";
        let y = b"AGCTCGATCGATCGATC";

        let error_profile = IlluminaR1ErrorProfile { x, y };
        let start_end_gap_params = Locality::Global;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&error_profile, &start_end_gap_params, None);

        dbg!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(17) * PROB_ILLUMINA_DEL.powi(2)
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(17) * PROB_ILLUMINA_DEL.powi(2),
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_mismatch() {
        let x = b"AGCTCGAGCGATCGATC";
        let y = b"TGCTCGATCGATCGATC";

        let error_profile = IlluminaR1ErrorProfile { x, y };
        let start_end_gap_params = Locality::Global;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&error_profile, &start_end_gap_params, None);

        dbg!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(15) * (PROB_ILLUMINA_SUBST / Prob(3.0)).powi(2)
        );
        assert!(*p <= 0.0);
        assert_relative_eq!(
            p.exp(),
            (Prob(1.0) - PROB_ILLUMINA_SUBST).powi(15) * (PROB_ILLUMINA_SUBST / Prob(3.0)).powi(2),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_large() {
        let x = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC\
ATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTAT\
CTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGT";

        let y = b"GGGTATGCACGCGATAGCATTGCGAGATGCTGGAGCTGGAGCACCCTATGTCGC";

        let error_profile = IlluminaR1ErrorProfile { x, y };
        let start_end_gap_params = Locality::Semiglobal;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&error_profile, &start_end_gap_params, None);

        let p_banded = pair_hmm.prob_related(&error_profile, &start_end_gap_params, Some(2));
        dbg!(p, p_banded);
        assert_relative_eq!(*p, *p_banded, epsilon = 1e-7);
    }
}
