use crate::stats::error_profiles::{
    Base, EmissionParameters, ErrorProfile, GapParameters, HomopolymerParameters, XYEmission,
};
use crate::stats::{LogProb, Prob};
use num_traits::Zero;

pub struct IlluminaR1ErrorProfile {
    pub x: &'static [u8],
    pub y: &'static [u8],
}

impl ErrorProfile for IlluminaR1ErrorProfile {
    type EmissionParams = IlluminaR1EmissionParams;
    type GapParams = IlluminaR1GapParams;
    type HomopolymerParams = IlluminaR1HomopolymerParams;

    fn emission_parameters(&self) -> Self::EmissionParams {
        IlluminaR1EmissionParams {
            x: self.x,
            y: self.y,
        }
    }

    fn gap_parameters(&self) -> Self::GapParams {
        IlluminaR1GapParams
    }

    fn homopolymer_parameters(&self) -> Self::HomopolymerParams {
        IlluminaR1HomopolymerParams
    }
}
// Single base insertion and deletion rates for R1 according to Schirmer et al.
// BMC Bioinformatics 2016, 10.1186/s12859-016-0976-y
pub static PROB_ILLUMINA_INS: Prob = Prob(2.8e-6);
pub static PROB_ILLUMINA_DEL: Prob = Prob(5.1e-6);
pub static PROB_ILLUMINA_SUBST: Prob = Prob(0.0021);

pub struct IlluminaR1EmissionParams {
    x: &'static [u8],
    y: &'static [u8],
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl EmissionParameters for IlluminaR1EmissionParams {
    fn prob_emit_x_and_y(&self, state: Base, i: usize, j: usize) -> XYEmission {
        if self.x[i] == self.y[j] {
            if BASES[state as usize] == self.x[i] {
                XYEmission::Match(LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST))
            } else {
                XYEmission::Match(LogProb::zero())
            }
        } else {
            if BASES[state as usize] == self.x[i] || BASES[state as usize] == self.y[j] {
                XYEmission::Mismatch(LogProb::from(PROB_ILLUMINA_SUBST))
            } else {
                XYEmission::Mismatch(LogProb::zero())
            }
        }
    }

    fn prob_emit_x_and_gap(&self, _state: Base, _: usize) -> LogProb {
        LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
    }

    fn prob_emit_gap_and_y(&self, _state: Base, _: usize) -> LogProb {
        LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
    }

    fn prob_emit_x_and_hop(&self, state: Base, i: usize) -> LogProb {
        if BASES[state as usize] == self.x[i] {
            LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
        } else {
            LogProb::zero()
        }
    }

    fn prob_emit_hop_and_y(&self, state: Base, j: usize) -> LogProb {
        if BASES[state as usize] == self.y[j] {
            LogProb::from(Prob(1.0) - PROB_ILLUMINA_SUBST)
        } else {
            LogProb::zero()
        }
    }

    fn len_x(&self) -> usize {
        self.x.len()
    }

    fn len_y(&self) -> usize {
        self.y.len()
    }
}

pub struct IlluminaR1HomopolymerParams;

impl HomopolymerParameters for IlluminaR1HomopolymerParams {
    fn prob_hop_x(&self, state: Base) -> LogProb {
        // match state {
        //     Base::A => LogProb::from(Prob(0.11)),
        //     Base::C => LogProb::from(Prob(0.07)),
        //     Base::G => LogProb::from(Prob(0.08)),
        //     Base::T => LogProb::from(Prob(0.12)),
        // }
        LogProb::from(Prob(1e-5))
    }

    fn prob_hop_y(&self, state: Base) -> LogProb {
        // match state {
        //     Base::A => LogProb::from(Prob(0.12)),
        //     Base::C => LogProb::from(Prob(0.13)),
        //     Base::G => LogProb::from(Prob(0.11)),
        //     Base::T => LogProb::from(Prob(0.11)),
        // }
        LogProb::from(Prob(1e-5))
    }

    fn prob_hop_x_extend(&self, state: Base) -> LogProb {
        match state {
            Base::A => LogProb::from(Prob(0.11)),
            Base::C => LogProb::from(Prob(0.07)),
            Base::G => LogProb::from(Prob(0.08)),
            Base::T => LogProb::from(Prob(0.12)),
        }
    }

    fn prob_hop_y_extend(&self, state: Base) -> LogProb {
        match state {
            Base::A => LogProb::from(Prob(0.12)),
            Base::C => LogProb::from(Prob(0.13)),
            Base::G => LogProb::from(Prob(0.11)),
            Base::T => LogProb::from(Prob(0.11)),
        }
    }
}

pub struct IlluminaR1GapParams;

impl GapParameters for IlluminaR1GapParams {
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
