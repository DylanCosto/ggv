//! Structured answers about alignments, computed straight from the data with no rendering.
//!
//! This is the layer agents talk to most: "what is at chr7:55,191,822?" comes back as numbers
//! (depth, per-allele strand counts, mapping quality, clipping), not as a screenshot.

pub mod chunks;
mod junctions;
mod pileup;
mod repeat;
mod summary;
mod sv;

pub use pileup::{
    AltAllele, BaseCounts, FetchReference, PileupOptions, ReadPredicate, ReadReference, ReadTest, PositionStats, ReadCounts, StrandCount, pileup, pileup_with_reads,
};
pub use ggv_model::splice_motif as motif;
pub use junctions::{Junction, JunctionOptions, Junctions, junctions};
pub use repeat::{RepeatAllele, RepeatLengths, RepeatRead, repeat_flank, repeat_lengths};
pub use summary::{RegionSummary, SummaryOptions, VariableSite, summarize};
pub use sv::{SvEvent, SvKind, SvSite, SvSupport, sv_events, sv_support};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Io(#[from] ggv_io::Error),

    #[error("pileup: {0}")]
    Pileup(#[from] ggv_io::pasta::pileup::Error),

    #[error("writing output: {0}")]
    Output(#[source] std::io::Error),
}

pub type Result<T> = std::result::Result<T, Error>;

/// A read's haplotype (`HP`): an integer tag, or a string holding a number (some tools write
/// `HP:Z:1`). Anything else is no haplotype; htslib's integer reading gave 0 for a string,
/// which showed as a group "HP 0".
pub fn haplotype(read: &ggv_io::pasta::pileup::Read) -> Option<i64> {
    match read.string_tag(*b"HP") {
        Some(s) => std::str::from_utf8(s).ok()?.trim().parse().ok(),
        None => read.integer_tag(*b"HP"),
    }
}
