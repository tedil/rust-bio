use crate::stats::LogProb;

pub mod illumina;

const NUM_BASES: usize = 4;

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    pub const fn values() -> [Base; NUM_BASES] {
        [Base::A, Base::C, Base::G, Base::T]
    }
}

pub trait ErrorProfile {
    type EmissionParams: EmissionParameters;
    type GapParams: GapParameters;
    type HomopolymerParams: HomopolymerParameters;
    fn emission_parameters(&self) -> Self::EmissionParams;
    fn gap_parameters(&self) -> Self::GapParams;
    fn homopolymer_parameters(&self) -> Self::HomopolymerParams;
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
pub enum Locality {
    Global,
    Semiglobal,
}

impl StartEndGapParameters for Locality {
    fn free_start_gap_x(&self) -> bool {
        match self {
            Locality::Global => false,
            Locality::Semiglobal => true,
        }
    }

    fn free_end_gap_x(&self) -> bool {
        match self {
            Locality::Global => false,
            Locality::Semiglobal => true,
        }
    }
}
