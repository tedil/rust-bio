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

use crate::stats::LogProb;
use itertools::Itertools;
use num_traits::Zero;
use smallvec::SmallVec;

fn cache_probs<F>(f: F) -> SmallVec<[LogProb; 4]>
where
    F: Fn(Base) -> LogProb,
{
    Base::values().iter().map(|&base| f(base)).collect()
}

/// Trait for parametrization of `PairHMM` gap behavior.
pub trait GapParameters {
    /// Probability to open gap in x.
    fn prob_gap_x(&self) -> LogProb;

    /// Probability to open gap in y.
    fn prob_gap_y(&self) -> LogProb;

    /// Probability to extend gap in x.
    fn prob_gap_x_extend(&self) -> LogProb;

    /// Probability to extend gap in y.
    fn prob_gap_y_extend(&self) -> LogProb;
}

/// Trait for parametrization of `PairHMM` gap behavior.
pub trait HomopolymerParameters {
    /// Probability to start a homopolymer run extension in x.
    fn prob_hop_x(&self, state: Base) -> LogProb;

    /// Probability to start a homopolymer run extension in y.
    fn prob_hop_y(&self, state: Base) -> LogProb;

    /// Probability to extend homopolymer extension in x.
    fn prob_hop_x_extend(&self, state: Base) -> LogProb;

    /// Probability to extend homopolymer extension in y.
    fn prob_hop_y_extend(&self, state: Base) -> LogProb;
}

/// Trait for parametrization of `PairHMM` start and end gap behavior.
/// This trait can be used to implement global and semiglobal alignments.
///
/// * global: methods return `false` and `LogProb::ln_zero()`.
/// * semiglobal: methods return `true` and `LogProb::ln_one()`.
pub trait StartEndGapParameters {
    /// Probability to start at x[i]. This can be left unchanged if you use `free_start_gap_x` and
    /// `free_end_gap_x`.
    #[inline]
    #[allow(unused_variables)]
    fn prob_start_gap_x(&self, i: usize) -> LogProb {
        if self.free_start_gap_x() {
            LogProb::ln_one()
        } else {
            // For global alignment, this has to return 0.0.
            LogProb::ln_zero()
        }
    }

    /// Allow free start gap in x.
    fn free_start_gap_x(&self) -> bool;

    /// Allow free end gap in x.
    fn free_end_gap_x(&self) -> bool;
}

pub enum XYEmission {
    Match(LogProb),
    Mismatch(LogProb),
}

impl XYEmission {
    pub fn prob(&self) -> LogProb {
        match self {
            &XYEmission::Match(p) => p,
            &XYEmission::Mismatch(p) => p,
        }
    }

    pub fn is_match(&self) -> bool {
        match self {
            &XYEmission::Match(_) => true,
            &XYEmission::Mismatch(_) => false,
        }
    }
}

/// Trait for parametrization of `PairHMM` emission behavior.
pub trait EmissionParameters {
    /// Emission probability for (x[i], y[j]).
    /// Returns a tuple with probability and a boolean indicating whether emissions match
    /// (e.g., are the same DNA alphabet letter).
    fn prob_emit_x_and_y(&self, state: Base, i: usize, j: usize) -> XYEmission;

    /// Emission probability for (x[i], -).
    fn prob_emit_x_and_gap(&self, state: Base, i: usize) -> LogProb;

    /// Emission probability for (-, y[j]).
    fn prob_emit_gap_and_y(&self, state: Base, j: usize) -> LogProb;

    /// Emission probability for (x[i], +).
    fn prob_emit_x_and_hop(&self, state: Base, i: usize) -> LogProb;

    /// Emission probability for (+, y[j]).
    fn prob_emit_hop_and_y(&self, state: Base, j: usize) -> LogProb;

    fn len_x(&self) -> usize;

    fn len_y(&self) -> usize;
}

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

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    const fn values() -> [Base; NUM_BASES] {
        [Base::A, Base::C, Base::G, Base::T]
    }
}

impl PairHHMM {
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
    pub fn prob_related<G, H, E>(
        &mut self,
        gap_params: &G,
        homopolymer_params: &H,
        emission_params: &E,
        max_edit_dist: Option<usize>,
    ) -> LogProb
    where
        G: GapParameters + StartEndGapParameters,
        H: HomopolymerParameters,
        E: EmissionParameters,
    {
        for k in 0..2 {
            self.fgx[k].clear();
            self.fgy[k].clear();
            for b in 0..NUM_BASES {
                self.fhx[k][b].clear();
                self.fhy[k][b].clear();
                self.fmm[k][b].clear();
            }
            self.min_edit_dist[k].clear();
            self.prob_cols.clear();

            let l = emission_params.len_y() + 1;
            self.fgx[k].resize(l, LogProb::ln_zero());
            self.fgy[k].resize(l, LogProb::ln_zero());
            for b in 0..NUM_BASES {
                self.fhx[k][b].resize(l, LogProb::ln_zero());
                self.fhy[k][b].resize(l, LogProb::ln_zero());
                self.fmm[k][b].resize(l, LogProb::ln_zero());
            }
            self.min_edit_dist[k].resize(l, usize::MAX);

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
        let prob_open_gap_x = gap_params.prob_gap_x();
        let prob_open_gap_y = gap_params.prob_gap_y();
        let prob_extend_gap_x = gap_params.prob_gap_x_extend();
        let prob_extend_gap_y = gap_params.prob_gap_y_extend();
        let extend_gap_in_y = prob_extend_gap_y != LogProb::ln_zero();
        let extend_gap_in_x = prob_extend_gap_x != LogProb::ln_zero();

        let no_hop_probs = cache_probs(|base| {
            homopolymer_params
                .prob_hop_x(base)
                .ln_add_exp(homopolymer_params.prob_hop_y(base))
                .ln_one_minus_exp()
        });
        let no_hop_x_extend_probs = cache_probs(|base| {
            homopolymer_params
                .prob_hop_x_extend(base)
                .ln_one_minus_exp()
        });
        let no_hop_y_extend_probs = cache_probs(|base| {
            homopolymer_params
                .prob_hop_y_extend(base)
                .ln_one_minus_exp()
        });
        let start_hop_x_probs = cache_probs(|base| homopolymer_params.prob_hop_x(base));
        let start_hop_y_probs = cache_probs(|base| homopolymer_params.prob_hop_y(base));
        let extend_hop_x_probs = cache_probs(|base| homopolymer_params.prob_hop_x_extend(base));
        let extend_hop_y_probs = cache_probs(|base| homopolymer_params.prob_hop_y_extend(base));

        let mut prev = 0;
        let mut curr = 1;
        for b in 0..NUM_BASES {
            self.fmm[prev][0][b] = LogProb::ln_one();
        }

        // iterate over x
        for i in 0..emission_params.len_x() {
            // allow alignment to start from offset in x (if prob_start_gap_x is set accordingly)
            for b in 0..NUM_BASES {
                self.fmm[prev][0][b] =
                    self.fmm[prev][0][b].ln_add_exp(gap_params.prob_start_gap_x(i));
            }

            if gap_params.free_start_gap_x() {
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
                        .map(|&base| {
                            emission_params.prob_emit_x_and_y(base, i, j).prob()
                                + LogProb::ln_sum_exp(&[
                                    // coming from one of the NUM_BASES match states
                                    prob_no_gap + no_hop_probs[0] + fmm_prev[j_minus_one][0],
                                    prob_no_gap + no_hop_probs[1] + fmm_prev[j_minus_one][1],
                                    prob_no_gap + no_hop_probs[2] + fmm_prev[j_minus_one][2],
                                    prob_no_gap + no_hop_probs[3] + fmm_prev[j_minus_one][3],
                                    // coming from one of the gap states
                                    prob_no_gap_x_extend + fgx_prev[j_minus_one],
                                    prob_no_gap_y_extend + fgy_prev[j_minus_one],
                                    // coming from one of the hop states
                                    prob_no_gap
                                        + no_hop_x_extend_probs[0]
                                        + fhx_prev[j_minus_one][0],
                                    prob_no_gap
                                        + no_hop_x_extend_probs[1]
                                        + fhx_prev[j_minus_one][1],
                                    prob_no_gap
                                        + no_hop_x_extend_probs[2]
                                        + fhx_prev[j_minus_one][2],
                                    prob_no_gap
                                        + no_hop_x_extend_probs[3]
                                        + fhx_prev[j_minus_one][3],
                                    prob_no_gap
                                        + no_hop_y_extend_probs[0]
                                        + fhy_prev[j_minus_one][0],
                                    prob_no_gap
                                        + no_hop_y_extend_probs[1]
                                        + fhy_prev[j_minus_one][1],
                                    prob_no_gap
                                        + no_hop_y_extend_probs[2]
                                        + fhy_prev[j_minus_one][2],
                                    prob_no_gap
                                        + no_hop_y_extend_probs[3]
                                        + fhy_prev[j_minus_one][3],
                                ])
                        })
                        .collect();

                    // gap in y
                    // TODO for each base: prob_gap_in_y[base] = ...
                    let mut prob_gap_in_y = emit_x_and_gap_probs[0]
                        + LogProb::ln_sum_exp(&[
                            // open gap from one of the match states
                            prob_open_gap_y + fmm_prev[j_][0],
                            prob_open_gap_y + fmm_prev[j_][1],
                            prob_open_gap_y + fmm_prev[j_][2],
                            prob_open_gap_y + fmm_prev[j_][3],
                            // open gap from one of the hop x states
                            prob_open_gap_y + fhx_prev[j_][0],
                            prob_open_gap_y + fhx_prev[j_][1],
                            prob_open_gap_y + fhx_prev[j_][2],
                            prob_open_gap_y + fhx_prev[j_][3],
                            // open gap from one of the hop y states
                            prob_open_gap_y + fhy_prev[j_][0],
                            prob_open_gap_y + fhy_prev[j_][1],
                            prob_open_gap_y + fhy_prev[j_][2],
                            prob_open_gap_y + fhy_prev[j_][3],
                        ]);
                    if extend_gap_in_y {
                        prob_gap_in_y = prob_gap_in_y.ln_add_exp(
                            // extend gap
                            prob_extend_gap_y + fgx_prev[j_],
                        );
                    }

                    // gap in x
                    // TODO for each base: prob_gap_in_x[base] = ...
                    let mut prob_gap_in_x = emission_params.prob_emit_gap_and_y(Base::A, j)
                        + LogProb::ln_sum_exp(&[
                            // open gap from one of the match states
                            prob_open_gap_x + fmm_curr[j_minus_one][0],
                            prob_open_gap_x + fmm_curr[j_minus_one][1],
                            prob_open_gap_x + fmm_curr[j_minus_one][2],
                            prob_open_gap_x + fmm_curr[j_minus_one][3],
                            // open gap from one of the hop x states
                            prob_open_gap_x + fhx_curr[j_minus_one][0],
                            prob_open_gap_x + fhx_curr[j_minus_one][1],
                            prob_open_gap_x + fhx_curr[j_minus_one][2],
                            prob_open_gap_x + fhx_curr[j_minus_one][3],
                            // open gap from one of the hop y states
                            prob_open_gap_x + fhy_curr[j_minus_one][0],
                            prob_open_gap_x + fhy_curr[j_minus_one][1],
                            prob_open_gap_x + fhy_curr[j_minus_one][2],
                            prob_open_gap_x + fhy_curr[j_minus_one][3],
                        ]);
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
                                +
                                // start hop
                                start_hop_y_probs[b] + fmm_prev[j_][b]
                                + if extend_hop_y_probs[b] != LogProb::zero() {
                                extend_hop_y_probs[b] + fhx_prev[j_][b]
                            } else { LogProb::zero() }
                        })
                        .collect();

                    // hop in x
                    let hop_in_x_probs: SmallVec<[LogProb; NUM_BASES]> = Base::values()
                        .iter()
                        .map(|&base| {
                            let b = base as usize;
                            emission_params.prob_emit_x_and_hop(base, i)
                                +
                                // start hop
                                start_hop_x_probs[b] + fmm_curr[j_minus_one][b]
                                + if extend_hop_x_probs[b] != LogProb::zero() {
                                extend_hop_x_probs[b] + fhy_curr[j_minus_one][b]
                            } else { LogProb::zero() }
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

                for &base in &Base::values() {
                    let b = base as usize;
                    self.fmm[curr][j_][b] = match_mismatch_probs[b];
                    self.fhx[curr][j_][b] = hop_in_x_probs[b];
                    self.fhy[curr][j_][b] = hop_in_y_probs[b];
                }
                self.fgx[curr][j_] = gap_in_y_prob;
                self.fgy[curr][j_] = gap_in_x_prob;
                if max_edit_dist.is_some() {
                    self.min_edit_dist[curr][j_] = min_edit_dist;
                }
            }

            if gap_params.free_end_gap_x() {
                // Cache column probabilities or simply record the last probability.
                // We can put all of them in one array since we simply have to sum in the end.
                // This is also good for numerical stability.
                for &base in &Base::values() {
                    self.prob_cols
                        .push(self.fmm[curr][base as usize].last().unwrap().clone());
                }
                self.prob_cols.push(self.fgx[curr].last().unwrap().clone());
                // TODO check removing this (we don't want open gaps in x):
                self.prob_cols.push(self.fgy[curr].last().unwrap().clone());
            }

            // next column
            mem::swap(&mut curr, &mut prev);
            // reset next column to zeros
            for v in &mut self.fmm[curr] {
                *v = vec![
                    LogProb::ln_zero(),
                    LogProb::ln_zero(),
                    LogProb::ln_zero(),
                    LogProb::ln_zero(),
                ];
            }
        }

        let p = if gap_params.free_end_gap_x() {
            LogProb::ln_sum_exp(&self.prob_cols)
        } else {
            LogProb::ln_sum_exp(&[
                *self.fmm[prev][0].last().unwrap(),
                *self.fmm[prev][1].last().unwrap(),
                *self.fmm[prev][2].last().unwrap(),
                *self.fmm[prev][3].last().unwrap(),
                *self.fgx[prev].last().unwrap(),
                *self.fgy[prev].last().unwrap(),
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
    use super::*;
    use crate::stats::{LogProb, Prob};

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
        fn prob_emit_x_and_y(&self, i: usize, j: usize) -> XYEmission {
            if self.x[i] == self.y[j] {
                XYEmission::Match(LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST))
            } else {
                XYEmission::Mismatch(LogProb::from(PROB_ILLUMINA_SUBST / Prob(3.0)))
            }
        }

        fn prob_emit_x_and_gap(&self, _: usize) -> LogProb {
            prob_emit_x_or_y()
        }

        fn prob_emit_gap_and_y(&self, _: usize) -> LogProb {
            prob_emit_x_or_y()
        }

        fn len_x(&self) -> usize {
            self.x.len()
        }

        fn len_y(&self) -> usize {
            self.y.len()
        }
    }

    struct TestGapParams;

    impl GapParameters for TestGapParams {
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

    impl StartEndGapParameters for TestGapParams {
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

    #[test]
    fn test_same() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(*p, 0.0, epsilon = 0.1);
    }

    #[test]
    fn test_insertion() {
        let x = b"AGCTCGATCGATCGATC";
        let y = b"AGCTCGATCTGATCGATCT";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(p.exp(), PROB_ILLUMINA_INS.powi(2), epsilon = 1e-11);
    }

    #[test]
    fn test_deletion() {
        let x = b"AGCTCGATCTGATCGATCT";
        let y = b"AGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(p.exp(), PROB_ILLUMINA_DEL.powi(2), epsilon = 1e-10);
    }

    #[test]
    fn test_mismatch() {
        let x = b"AGCTCGAGCGATCGATC";
        let y = b"TGCTCGATCGATCGATC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = TestGapParams;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        assert!(*p <= 0.0);
        assert_relative_eq!(
            p.exp(),
            (PROB_ILLUMINA_SUBST / Prob(3.0)).powi(2),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_large() {
        let x = b"GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGC\
ATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTAT\
CTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGT";

        let y = b"GGGTATGCACGCGATAGCATTGCGAGATGCTGGAGCTGGAGCACCCTATGTCGC";

        let emission_params = TestEmissionParams { x, y };
        let gap_params = SemiglobalGapParams;

        let mut pair_hmm = PairHHMM::new();
        let p = pair_hmm.prob_related(&gap_params, &emission_params, None);

        let p_banded = pair_hmm.prob_related(&gap_params, &emission_params, Some(2));

        assert_relative_eq!(*p, *p_banded, epsilon = 1e-7);
    }
}
