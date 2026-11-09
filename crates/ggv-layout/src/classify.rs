//! What a read is, for colouring: its pair's orientation and insert size, and the methylation
//! likelihood at each called base. Shared by the GPU view and the headless renderer, so a
//! picture from `ggv render` and the app colour the same reads the same way.

use crate::AlignedRead;

/// How a read's pair sits, relative to the expected forward-reverse (LR) layout and the usual
/// insert size.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PairClass {
    /// The mate is on another contig.
    OtherContig,
    /// Reverse-forward (outward facing): duplications, circular templates.
    RightLeft,
    /// Both forward: inversions.
    LeftLeft,
    /// Both reverse: inversions.
    RightRight,
    /// Insert larger than nearly all others in view: deletions.
    Large,
    /// Insert smaller than nearly all others: insertions.
    Small,
}

/// The insert sizes beyond which a pair is unusual (0.5th and 99.5th percentiles), from the
/// reads in view; `None` with under 100 pairs (too few to say).
pub fn insert_bounds(reads: &[AlignedRead]) -> Option<(u32, u32)> {
    let mut sizes: Vec<u32> = reads
        .iter()
        .filter(|r| r.mate.is_some_and(|m| m.same_contig) && r.template_length != 0)
        .map(|r| r.template_length.unsigned_abs())
        .collect();
    if sizes.len() < 100 {
        return None;
    }
    sizes.sort_unstable();
    let at = |q: f64| sizes[((sizes.len() - 1) as f64 * q) as usize];
    Some((at(0.005), at(0.995)))
}

/// The read's pair class, or `None` for an unpaired read or an ordinary pair.
pub fn pair_class(r: &AlignedRead, bounds: Option<(u32, u32)>) -> Option<PairClass> {
    let mate = r.mate?;
    if !mate.same_contig {
        return Some(PairClass::OtherContig);
    }
    let read_is_left = (r.start as i64) <= mate.pos;
    let (left_rev, right_rev) = if read_is_left { (r.reverse, mate.reverse) } else { (mate.reverse, r.reverse) };
    // Mates that overlap (short inserts, adapter read-through) can start in either order, so
    // start positions say nothing about RL there; only apart mates can be outward-facing.
    let overlapping = (r.start as i64 - mate.pos).unsigned_abs() < r.end - r.start;
    match (left_rev, right_rev) {
        (false, true) => {}
        (true, false) if overlapping => {}
        (true, false) => return Some(PairClass::RightLeft),
        (false, false) => return Some(PairClass::LeftLeft),
        (true, true) => return Some(PairClass::RightRight),
    }
    let (low, high) = bounds?;
    let tlen = r.template_length.unsigned_abs();
    if tlen > high {
        Some(PairClass::Large)
    } else if tlen < low && tlen > 0 {
        Some(PairClass::Small)
    } else {
        None
    }
}

/// Each called base's modification likelihood (all codes summed, capped at 255), in position
/// order: the calls for one base sit together in `ModCall` order.
pub fn mod_likelihoods(r: &AlignedRead) -> impl Iterator<Item = (u64, u8)> + '_ {
    let mut i = 0;
    std::iter::from_fn(move || {
        let first = r.mods.get(i)?;
        let mut sum = 0u32;
        while let Some(m) = r.mods.get(i).filter(|m| m.offset == first.offset) {
            sum += u32::from(m.prob);
            i += 1;
        }
        Some((first.pos(r.start), sum.min(255) as u8))
    })
}

/// IGV's two-colour modification scale as (red, green, blue, alpha): red for likely modified,
/// blue for likely not, more opaque the surer the call.
pub fn mod_rgba(likelihood: u8) -> (u8, u8, u8, u8) {
    if likelihood >= 128 {
        (220, 20, 20, ((u32::from(likelihood) - 128) * 2).min(255) as u8)
    } else {
        (20, 20, 220, ((128 - u32::from(likelihood)) * 2).min(255) as u8)
    }
}

/// A colour for a tag value (`CB`, `PS`...): the same value always gets the same colour, in the
/// app and in headless pictures. Twelve well-separated hues; values are hashed (FNV-1a).
pub fn tag_rgb(value: &str) -> (u8, u8, u8) {
    const PALETTE: [(u8, u8, u8); 12] = [
        (228, 26, 28), (55, 126, 184), (77, 175, 74), (152, 78, 163), (255, 127, 0), (166, 86, 40),
        (247, 129, 191), (0, 158, 160), (190, 174, 30), (102, 102, 230), (140, 200, 90), (200, 60, 120),
    ];
    let mut h: u32 = 0x811c_9dc5;
    for b in value.bytes() {
        h = (h ^ u32::from(b)).wrapping_mul(0x0100_0193);
    }
    PALETTE[(h % PALETTE.len() as u32) as usize]
}
