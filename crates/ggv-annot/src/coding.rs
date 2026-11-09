//! What a base or a change means for a transcript: where it falls (exon, intron, UTR, splice
//! site), its HGVS `c.` position, the codon and amino acid, and the consequence of a change
//! (`c.1234A>G p.Lys412Arg`, missense).
//!
//! Positions are 0-based genomic; transcript positions follow the direction of transcription.
//! Reference bases come from a caller-supplied lookup, so this works from a loaded view, a
//! FASTA, or a test string alike. Everything is O(exons) per call, with no allocation beyond
//! the strings returned.

use crate::{Strand, Transcript, without_copy_suffix};

/// The standard genetic code, indexed by `16 * b1 + 4 * b2 + b3` with T=0, C=1, A=2, G=3.
const CODE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/// The vertebrate mitochondrial code (NCBI table 2): TGA Trp, ATA Met, AGA/AGG stop.
const CODE_MITO: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// One-letter amino acid for a codon (`*` stop, `X` if any base is unknown), standard code.
pub fn translate(codon: [u8; 3]) -> u8 {
    translate_with(1, codon)
}

/// Like [`translate`], with NCBI genetic code `table` (2 is vertebrate mitochondrial; others
/// read as standard).
pub fn translate_with(table: u8, codon: [u8; 3]) -> u8 {
    let code = if table == 2 { CODE_MITO } else { CODE };
    match (base_index(codon[0]), base_index(codon[1]), base_index(codon[2])) {
        (Some(a), Some(b), Some(c)) => code[16 * a + 4 * b + c],
        _ => b'X',
    }
}

/// The amino acid codon `k` of `t` encodes, read as `codon` (on the transcript's strand): the
/// transcript's genetic code, and `U` for a TGA the annotation marks as selenocysteine.
pub fn amino_acid(t: &Transcript, k: u64, codon: [u8; 3]) -> u8 {
    let aa = translate_with(t.translation.table, codon);
    if aa == b'*' && codon.eq_ignore_ascii_case(b"TGA") && is_selenocysteine(t, k) { b'U' } else { aa }
}

fn is_selenocysteine(t: &Transcript, k: u64) -> bool {
    let sec = &t.translation.selenocysteine;
    !sec.is_empty() && codon_positions(t, k).is_some_and(|cp| cp.iter().any(|p| sec.iter().any(|&(s, e)| (s..e).contains(p))))
}

/// Three-letter name of a one-letter amino acid (`Ter` for stop).
pub fn three_letter(aa: u8) -> &'static str {
    match aa {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'Q' => "Gln",
        b'E' => "Glu",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        b'U' => "Sec",
        b'*' => "Ter",
        _ => "Xaa",
    }
}

pub fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn reverse(t: &Transcript) -> bool {
    t.strand == Strand::Reverse
}

/// Exons in transcription order.
fn exons_tx(t: &Transcript) -> impl Iterator<Item = (u64, u64)> + '_ {
    let n = t.exons.len();
    (0..n).map(move |i| if reverse(t) { t.exons[n - 1 - i] } else { t.exons[i] })
}

/// Position within the spliced transcript (0-based, from the 5' end), if `pos` is exonic.
pub fn tx_offset(t: &Transcript, pos: u64) -> Option<u64> {
    let mut acc = 0;
    for (s, e) in exons_tx(t) {
        if (s..e).contains(&pos) {
            return Some(acc + if reverse(t) { e - 1 - pos } else { pos - s });
        }
        acc += e - s;
    }
    None
}

/// The genomic position at transcript offset `off`.
pub fn tx_pos(t: &Transcript, mut off: u64) -> Option<u64> {
    for (s, e) in exons_tx(t) {
        let len = e - s;
        if off < len {
            return Some(if reverse(t) { e - 1 - off } else { s + off });
        }
        off -= len;
    }
    None
}

/// The coding span in transcript offsets: first coding base, one past the last.
fn cds_tx(t: &Transcript) -> Option<(u64, u64)> {
    let (cs, ce) = t.cds?;
    let (first, last) = if reverse(t) { (ce - 1, cs) } else { (cs, ce - 1) };
    Some((tx_offset(t, first)?, tx_offset(t, last)? + 1))
}

/// The coding span read in codons: from the first complete codon (past the phase of an
/// incomplete start) to one past the last coding base, in transcript offsets.
fn frame_tx(t: &Transcript) -> Option<(u64, u64)> {
    let (a, b) = cds_tx(t)?;
    let a = a + u64::from(t.translation.phase);
    (a < b).then_some((a, b))
}

/// Number of whole codons in the coding sequence.
fn codon_count(t: &Transcript) -> u64 {
    frame_tx(t).map_or(0, |(a, b)| (b - a) / 3)
}

/// 0-based index into the coding sequence from its first complete codon, if `pos` is a coding
/// base (none for the bases before that codon in an incomplete coding sequence).
pub fn coding_index(t: &Transcript, pos: u64) -> Option<u64> {
    let (a, b) = frame_tx(t)?;
    let o = tx_offset(t, pos)?;
    (a..b).contains(&o).then(|| o - a)
}

/// Genomic positions of codon `k` (0-based), in transcription order.
pub fn codon_positions(t: &Transcript, k: u64) -> Option<[u64; 3]> {
    let (a, b) = frame_tx(t)?;
    let at = |i: u64| (a + 3 * k + i < b).then(|| tx_pos(t, a + 3 * k + i)).flatten();
    Some([at(0)?, at(1)?, at(2)?])
}

/// Where a position falls in a transcript.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Location {
    /// `exon 3/11` or `intron 2/10`, numbered in transcription order.
    pub feature: String,
    /// HGVS `c.` position without the prefix: `123`, `-14`, `*20`, `123+5`, `124-2`; for
    /// non-coding transcripts the transcript position (`n.`).
    pub c: String,
    pub coding: bool,
    pub kind: Region,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Region {
    Utr5,
    Coding,
    Utr3,
    NonCodingExon,
    SpliceDonor,
    SpliceAcceptor,
    /// Intronic, 3–8 bases from an exon.
    SpliceRegion,
    Intron,
}

impl Region {
    pub fn name(self) -> &'static str {
        match self {
            Self::Utr5 => "5' UTR",
            Self::Coding => "coding",
            Self::Utr3 => "3' UTR",
            Self::NonCodingExon => "non-coding exon",
            Self::SpliceDonor => "splice donor",
            Self::SpliceAcceptor => "splice acceptor",
            Self::SpliceRegion => "splice region",
            Self::Intron => "intron",
        }
    }
}

/// The `c.` number for an exonic transcript offset.
fn c_number(t: &Transcript, off: u64) -> String {
    match cds_tx(t) {
        Some((a, _)) if off < a => format!("-{}", a - off),
        Some((_, b)) if off >= b => format!("*{}", off - b + 1),
        Some((a, _)) => (off - a + 1).to_string(),
        None => (off + 1).to_string(),
    }
}

/// Where `pos` is in `t`, or `None` outside the transcript.
pub fn locate(t: &Transcript, pos: u64) -> Option<Location> {
    if !(t.start..t.end).contains(&pos) {
        return None;
    }
    let n = t.exons.len();
    let exons: Vec<(u64, u64)> = exons_tx(t).collect();
    if let Some(i) = exons.iter().position(|&(s, e)| (s..e).contains(&pos)) {
        let off = tx_offset(t, pos)?;
        let kind = match cds_tx(t) {
            None => Region::NonCodingExon,
            Some((a, _)) if off < a => Region::Utr5,
            Some((_, b)) if off >= b => Region::Utr3,
            Some(_) => Region::Coding,
        };
        return Some(Location { feature: format!("exon {}/{n}", i + 1), c: c_number(t, off), coding: kind == Region::Coding, kind });
    }
    // Intronic: between exon i and i + 1 in transcription order.
    for i in 0..n.saturating_sub(1) {
        let (a, b) = (exons[i], exons[i + 1]);
        let (lo, hi) = if reverse(t) { (b.1, a.0) } else { (a.1, b.0) };
        if !(lo..hi).contains(&pos) {
            continue;
        }
        // Distance from the end of exon i (donor side) and to the start of exon i + 1.
        let (after_donor, before_acceptor) = if reverse(t) { (a.0 - pos, pos - b.1 + 1) } else { (pos - a.1 + 1, b.0 - pos) };
        let last_of_i = tx_offset(t, if reverse(t) { a.0 } else { a.1 - 1 })?;
        let first_of_next = last_of_i + 1;
        let c = if after_donor <= before_acceptor {
            format!("{}+{after_donor}", c_number(t, last_of_i))
        } else {
            format!("{}-{before_acceptor}", c_number(t, first_of_next))
        };
        let kind = match (after_donor, before_acceptor) {
            (d, _) if d <= 2 => Region::SpliceDonor,
            (_, d) if d <= 2 => Region::SpliceAcceptor,
            (d, e) if d <= 8 || e <= 8 => Region::SpliceRegion,
            _ => Region::Intron,
        };
        return Some(Location { feature: format!("intron {}/{}", i + 1, n - 1), c, coding: false, kind });
    }
    None
}

/// What a change does to a transcript.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Consequence {
    pub transcript: String,
    pub gene: String,
    pub location: Location,
    /// e.g. `missense`, `synonymous`, `stop gained`, `frameshift`, `splice donor`, `intron`.
    pub effect: String,
    /// `c.1234A>G` (bases on the transcript's strand).
    pub hgvs_c: String,
    /// `p.Lys412Arg`, `p.Lys412=`, `p.Arg13Ter` for coding substitutions; `p.Phe508del`,
    /// `p.Gly5AlafsTer3` and the like for coding indels.
    pub hgvs_p: Option<String>,
    /// Reference and changed codon, upper case for the changed base.
    pub codons: Option<(String, String)>,
    /// Why the result may be wrong or incomplete: a REF that disagrees with the reference, an
    /// effect that isn't predicted, an annotation missing what's needed.
    pub warnings: Vec<String>,
}

/// The effect when the protein change can't be named with confidence.
pub const UNCERTAIN: &str = "coding (protein effect uncertain)";
/// The effect of a change in the stop codon that leaves it a stop where it was.
const STOP_RETAINED: &str = "stop retained";

/// Whether an allele is a base sequence (A, C, G, T, N), not a symbolic allele (`<DEL>`), a
/// breakend (`G]chr2:100]`), a missing allele (`.`) or a spanning deletion (`*`).
pub fn is_sequence(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

/// How many REF bases are checked against the reference.
const REF_CHECK: usize = 1000;

fn internal_stop_warning(k: u64) -> String {
    format!(
        "codon {} is a stop inside the coding sequence (a selenocysteine or other translation exception the annotation doesn't give?): the protein change isn't predicted",
        k + 1
    )
}

/// The consequence of changing `reference` (as written in a VCF, starting at `pos`) to `alt`.
/// `None` outside the transcript, or when either allele isn't a sequence (see
/// [`is_sequence`]): a symbolic allele says nothing about which bases change.
/// `base(p)` gives the reference base at genomic `p` (needed for the codon of a coding
/// substitution; `None` leaves the protein change unknown).
pub fn consequence(t: &Transcript, pos: u64, reference: &str, alt: &str, base: impl Fn(u64) -> Option<u8>) -> Option<Consequence> {
    // ALT equal to REF is no change: nothing to name (`c.1799T>T` isn't HGVS).
    if !is_sequence(reference) || !is_sequence(alt) || reference.eq_ignore_ascii_case(alt) {
        return None;
    }
    let rc = |s: &str| -> String { s.bytes().rev().map(|b| complement(b) as char).collect() };
    let strand = |s: &str| if reverse(t) { rc(s) } else { s.to_ascii_uppercase() };
    let (r, a) = (reference.to_ascii_uppercase(), alt.to_ascii_uppercase());

    let mut warnings = Vec::new();
    // The consequence is read from REF: say so when the reference disagrees.
    let found: Option<String> = (0..r.len().min(REF_CHECK) as u64).map(|i| base(pos + i).map(|b| b.to_ascii_uppercase() as char)).collect();
    if let Some(found) = found
        && !r.starts_with(&found)
    {
        warnings.push(format!("REF {r} does not match the reference here ({found}): this consequence assumes the VCF's REF"));
    }

    // Trim the bases REF and ALT share (the VCF anchor, and the padding of non-minimal or split
    // multi-allelic records: `ATCTT>AT` is `TCT>-`, `ACC>ATC` is `C>T`), first then last, so
    // the change starts where it acts and the same change gets one name.
    let (rb, ab) = (r.as_bytes(), a.as_bytes());
    let pre = rb.iter().zip(ab).take_while(|(x, y)| x == y).count();
    let suf = rb[pre..].iter().rev().zip(ab[pre..].iter().rev()).take_while(|(x, y)| x == y).count();
    let (pos, r, a) = (pos + pre as u64, r[pre..r.len() - suf].to_owned(), a[pre..a.len() - suf].to_owned());
    // HGVS places an indel in a repeat at its most 3' position on the transcript: shift it
    // through the repeat (up the genome on +, down it on −), within its exon or intron.
    let Shifted { pos, r, a, first, across } = shift_3prime(t, pos, r, a, &base);
    let del = r.len() as u64;
    // A change reaching past an end of the transcript: what it does isn't worked out.
    if del > 0 && (pos < t.start || pos + del > t.end) {
        return past_ends(t, pos, &r, &strand(&a), warnings);
    }
    // A pure insertion sits between pos - 1 and pos; locate it at the base before.
    let anchor = if r.is_empty() { pos.saturating_sub(1) } else { pos };
    let location = locate(t, anchor)?;
    if let Some(exon) = across {
        let (here, there) = if exon { ("intron", "exon") } else { ("exon", "intron") };
        warnings.push(format!(
            "the change is in a repeat running on across an exon-intron junction: it is named in the {here} it starts in (not shifted across the junction), but the same sequence change is equally consistent with one in the next {there}"
        ));
    }
    if t.exons_from_span {
        warnings.push("the annotation gives no exons for this transcript: its whole span is taken as one exon".into());
    }
    if t.translation.phase > 0 && t.cds.is_some() && location.kind == Region::Coding {
        warnings.push(format!(
            "the coding sequence is incomplete at its 5' end (phase {}): codons are read from the first complete one, and c./p. numbers count from the annotated start",
            t.translation.phase
        ));
    }
    let mut out = Consequence {
        transcript: without_copy_suffix(&t.id).to_owned(),
        gene: t.gene_name.to_string(),
        effect: location.kind.name().to_owned(),
        hgvs_c: String::new(),
        hgvs_p: None,
        codons: None,
        location,
        warnings,
    };
    let prefix = if t.cds.is_some() { "c." } else { "n." };

    if r.len() == 1 && a.len() == 1 {
        out.hgvs_c = format!("{prefix}{}{}>{}", out.location.c, strand(&r), strand(&a));
        if out.location.kind != Region::Coding {
            return Some(out);
        }
        let Some(ci) = coding_index(t, pos) else {
            // Before the first complete codon of an incomplete coding sequence.
            out.hgvs_p = Some("p.?".into());
            out.effect = "coding (incomplete codon)".into();
            return Some(out);
        };
        let k = ci / 3;
        let phase = (ci % 3) as usize;
        let Some(cp) = codon_positions(t, k) else { return Some(out) };
        let mut codon = [0u8; 3];
        for (i, &p) in cp.iter().enumerate() {
            let b = if p == pos { r.as_bytes()[0] } else { base(p).map(|b| b.to_ascii_uppercase()).unwrap_or(b'N') };
            codon[i] = if reverse(t) { complement(b) } else { b };
        }
        let mut changed = codon;
        changed[phase] = if reverse(t) { complement(a.as_bytes()[0]) } else { a.as_bytes()[0] };
        let (ra, aa) = (amino_acid(t, k, codon), amino_acid(t, k, changed));
        let show = |c: [u8; 3]| -> String {
            c.iter().enumerate().map(|(i, &b)| if i == phase { b as char } else { (b as char).to_ascii_lowercase() }).collect()
        };
        out.codons = Some((show(codon), show(changed)));
        let n = k + 1;
        if ra == b'X' || aa == b'X' {
            out.hgvs_p = Some("p.?".into());
            out.effect = "coding (codon unknown)".into();
        } else if ra == b'*' && k + 1 < codon_count(t) {
            // A stop before the last codon: the annotation is missing how it is read.
            out.hgvs_p = Some("p.?".into());
            out.effect = UNCERTAIN.into();
            out.warnings.push(internal_stop_warning(k));
        } else {
            out.effect = match (ra, aa) {
                (x, y) if x == y => "synonymous",
                (_, b'*') => "stop gained",
                (b'*', _) => "stop lost",
                _ if k == 0 && ra == b'M' && t.translation.phase == 0 => "start lost",
                _ => "missense",
            }
            .into();
            out.hgvs_p = Some(if ra == aa {
                format!("p.{}{n}=", three_letter(ra))
            } else if ra == b'*' {
                // A stop read on: named with its extension, as for longer changes (one code path).
                protein_change(t, pos, &r, &strand(&a), &base).map_or_else(|| format!("p.Ter{n}{}", three_letter(aa)), |pc| pc.p)
            } else {
                format!("p.{}{n}{}", three_letter(ra), three_letter(aa))
            });
        }
        return Some(out);
    }

    // Insertions and deletions (and longer substitutions).
    let ins = a.len() as i64;
    let end_loc = if del > 1 { locate(t, pos + del - 1) } else { None };
    let span = match (&end_loc, reverse(t)) {
        (Some(e), false) => format!("{}_{}", out.location.c, e.c),
        (Some(e), true) => format!("{}_{}", e.c, out.location.c),
        (None, _) => out.location.c.clone(),
    };
    // An insertion repeating the bases just before it (on the transcript) is a dup of them: on
    // + those end at pos - 1; on − they start at pos (3' is down the genome).
    let dup = if del == 0 {
        let n = a.len() as u64;
        let span = if reverse(t) { pos..pos + n } else { pos.saturating_sub(n)..pos };
        // Those bases must lie in the insertion's own exon or intron (as in HGVS, a dup in an
        // intron is one too): near an exon's edge the genomic span would take in intron, which
        // isn't on the transcript.
        let in_one = n > 0 && shift_bounds(t, pos, 0).is_some_and(|(lo, hi)| lo <= span.start && span.end <= hi);
        let repeated = span.start + n == span.end
            && in_one
            && span.clone().zip(a.bytes()).all(|(p, b)| base(p).is_some_and(|x| x.eq_ignore_ascii_case(&b)));
        if repeated { locate(t, span.start).zip(locate(t, span.end - 1)).map(|(lo, hi)| (lo, hi, span)) } else { None }
    } else {
        None
    };
    out.hgvs_c = match (del, ins) {
        (0, _) => {
            let n = a.len() as u64;
            match &dup {
                Some((lo, hi, _)) => {
                    let (x, y) = if reverse(t) { (&hi.c, &lo.c) } else { (&lo.c, &hi.c) };
                    if n == 1 { format!("{prefix}{x}dup") } else { format!("{prefix}{x}_{y}dup") }
                }
                None => {
                    // Between the anchor and the next base on the transcript.
                    let next = locate(t, pos).map_or_else(|| "?".to_owned(), |l| l.c);
                    let (x, y) = if reverse(t) { (next, out.location.c.clone()) } else { (out.location.c.clone(), next) };
                    format!("{prefix}{x}_{y}ins{}", strand(&a))
                }
            }
        }
        (_, 0) => format!("{prefix}{span}del"),
        _ => format!("{prefix}{span}delins{}", strand(&a)),
    };
    // Where an insertion is read. Shifted 3' up to the end of its exon (or intron), beside the
    // junction, it is still the same change one base 5', inside that exon (or intron): read
    // it there, where its effect matches its name (an exon's last base duplicated lengthens
    // the exon). Anywhere else it is read where it is named.
    let (mut ins_at, mut ins_seq) = (pos, a.clone());
    if del == 0
        && pos != first
        && let Some((rs, re)) = shift_bounds(t, pos, 0)
    {
        let mut d = a.clone().into_bytes();
        let at = if reverse(t) && pos == rs && rs != t.start {
            d.rotate_left(1);
            pos + 1
        } else if !reverse(t) && pos == re && re != t.end {
            d.rotate_right(1);
            pos - 1
        } else {
            pos
        };
        // Not into a start codon it otherwise leaves alone (an exon ending at c.3).
        if !alters_start(t, at, 0, &base) || alters_start(t, pos, 0, &base) {
            (ins_at, ins_seq) = (at, String::from_utf8(d).unwrap_or_default());
        }
    }
    let coding = if del == 0 {
        // An insertion changes the coding sequence only with coding bases on both sides.
        let sides = |pos: u64| {
            let (before, after) = if reverse(t) { (Some(pos), pos.checked_sub(1)) } else { (pos.checked_sub(1), Some(pos)) };
            let kind = |p: Option<u64>| p.and_then(|p| locate(t, p)).map(|l| l.kind);
            (kind(before), kind(after))
        };
        match sides(ins_at) {
            (Some(Region::Coding), Some(Region::Coding)) => true,
            (Some(Region::Coding), Some(k)) | (Some(k), Some(Region::Coding)) => {
                out.effect = k.name().to_owned();
                if !matches!(k, Region::Utr5 | Region::Utr3) {
                    out.warnings.push(
                        "the insertion is at an exon-intron junction: its effect on splicing and the protein isn't predicted".into(),
                    );
                }
                false
            }
            // Otherwise what it is inserted into: as the base 5' of it on the transcript (on −
            // the base `location` gives is the one 3' of it).
            (Some(k), _) => {
                out.effect = k.name().to_owned();
                false
            }
            _ => false,
        }
    } else {
        out.location.kind == Region::Coding || end_loc.as_ref().is_some_and(|l| l.kind == Region::Coding)
    };
    if coding {
        let net = ins - del as i64;
        out.effect = if net % 3 != 0 {
            "frameshift".into()
        } else if net > 0 {
            "inframe insertion".into()
        } else if net < 0 {
            "inframe deletion".into()
        } else {
            "inframe substitution".into()
        };
    }
    // A deletion taking any splice-site base counts as one, wherever its ends fall.
    let splice = if del > 0 { splice_hit(t, pos, pos + del) } else { None };
    if let Some(k) = splice {
        out.effect = k.name().to_owned();
        if t.exons.iter().any(|&(s, e)| (pos < s && s < pos + del) || (pos < e && e < pos + del)) {
            out.warnings.push(
                "the deletion spans an exon-intron boundary, taking exonic and splice-site bases: its effect on the protein isn't predicted"
                    .into(),
            );
        }
    } else if alters_start(t, ins_at, del, &base) {
        if before_start(t, first, del) {
            // Named on the start codon by the 3' rule, but the same change made in the 5' UTR
            // (an A more or less in a run ending at c.1): the ATG and all after it are kept.
            out.effect = Region::Utr5.name().into();
            out.hgvs_p = None;
        } else {
            // Any change to the start codon (as for an SNV there): what is translated isn't known.
            out.effect = "start lost".into();
            out.hgvs_p = Some("p.?".into());
        }
    } else if let Some(pc) = protein_change(t, ins_at, &r, &strand(&ins_seq), &base) {
        match pc.effect {
            Some(e) if e == UNCERTAIN || e == STOP_RETAINED => out.effect = e.into(),
            Some(e) if out.effect.starts_with("inframe") => out.effect = e.into(),
            _ => {}
        }
        out.warnings.extend(pc.warning);
        out.hgvs_p = Some(pc.p);
    }
    Some(out)
}

/// The splice site (the two intron bases at either end of an intron) that genomic `[lo, hi)`
/// takes in, donor first.
fn splice_hit(t: &Transcript, lo: u64, hi: u64) -> Option<Region> {
    let mut acceptor = false;
    for w in t.exons.windows(2) {
        let (is, ie) = (w[0].1, w[1].0);
        let left = (is, (is + 2).min(ie));
        let right = (ie.saturating_sub(2).max(is), ie);
        let hit = |(s, e): (u64, u64)| lo < e && s < hi;
        // On + the intron's left end is the donor; on − its right end.
        let (donor, acc) = if reverse(t) { (right, left) } else { (left, right) };
        if hit(donor) {
            return Some(Region::SpliceDonor);
        }
        acceptor |= hit(acc);
    }
    acceptor.then_some(Region::SpliceAcceptor)
}

/// Whether deleting `del` bases at `pos` (or inserting before `pos`, `del` 0) changes the Met
/// start codon of a complete coding sequence: takes one of its bases, or goes between two.
fn alters_start(t: &Transcript, pos: u64, del: u64, base: &impl Fn(u64) -> Option<u8>) -> bool {
    let Some(cp) = codon_positions(t, 0).filter(|_| t.translation.phase == 0) else { return false };
    let mut codon = [0u8; 3];
    for (c, &p) in codon.iter_mut().zip(&cp) {
        let b = base(p).map_or(b'N', |b| b.to_ascii_uppercase());
        *c = if reverse(t) { complement(b) } else { b };
    }
    if amino_acid(t, 0, codon) != b'M' {
        return false;
    }
    if del > 0 {
        cp.iter().any(|p| (pos..pos + del).contains(p))
    } else {
        pos.checked_sub(1).is_some_and(|b| cp.contains(&b)) && cp.contains(&pos)
    }
}

/// Whether deleting `del` bases at `pos` (or inserting before `pos` on the transcript, `del`
/// 0) is wholly in the 5' UTR: before the start codon, which it leaves as it was.
fn before_start(t: &Transcript, pos: u64, del: u64) -> bool {
    let utr5 = |p: u64| locate(t, p).is_some_and(|l| l.kind == Region::Utr5);
    match del {
        0 if reverse(t) => utr5(pos),
        0 => pos.checked_sub(1).is_some_and(utr5),
        _ => (pos..pos + del).all(utr5),
    }
}

/// A deletion or longer change reaching past an end of `t`: named by the part on the
/// transcript, `?` for the end outside it, with no effect predicted. `None` if it misses `t`.
fn past_ends(t: &Transcript, pos: u64, r: &str, ins: &str, mut warnings: Vec<String>) -> Option<Consequence> {
    let end = pos + r.len() as u64;
    let (lo, hi) = (pos.max(t.start), end.min(t.end));
    if lo >= hi {
        return None;
    }
    let (l, h) = (locate(t, lo)?, locate(t, hi - 1)?);
    let x = if pos < t.start { "?".to_owned() } else { l.c.clone() };
    let y = if end > t.end { "?".to_owned() } else { h.c.clone() };
    let (x, y) = if reverse(t) { (y, x) } else { (x, y) };
    let prefix = if t.cds.is_some() { "c." } else { "n." };
    let hgvs_c = if ins.is_empty() { format!("{prefix}{x}_{y}del") } else { format!("{prefix}{x}_{y}delins{ins}") };
    warnings.push("the change reaches past an end of the transcript: its effect isn't predicted".into());
    Some(Consequence {
        transcript: without_copy_suffix(&t.id).to_owned(),
        gene: t.gene_name.to_string(),
        effect: "transcript boundary (effect unknown)".into(),
        hgvs_c,
        hgvs_p: None,
        codons: None,
        // The 5'-most base of the transcript it takes.
        location: if reverse(t) { h } else { l },
        warnings,
    })
}

/// How far a frameshift or stop loss is followed looking for the new stop codon.
const FRAMESHIFT_CODONS: u64 = 3000;

/// A named protein change, with the effect it implies when the change's length doesn't say.
struct Protein {
    p: String,
    /// `stop gained`, `stop lost`, `missense`, `synonymous` or [`UNCERTAIN`].
    effect: Option<&'static str>,
    warning: Option<String>,
}

/// The protein change of an indel or longer substitution that starts in the coding sequence
/// and stays on the exon (it may run on into the 3' UTR): `p.Phe508del`, `p.Phe4dup`,
/// `p.Glu3_Ile4delinsVal`, `p.Gly5AlafsTer3`, `p.Ter3GlnextTer2`, `p.?` when the start codon
/// changes. `r` is the deleted genomic sequence at `pos` (shifted 3'), `ins` the inserted bases
/// on the transcript's strand.
fn protein_change(t: &Transcript, pos: u64, r: &str, ins: &str, base: &impl Fn(u64) -> Option<u8>) -> Option<Protein> {
    let (cds_start, _) = frame_tx(t)?;
    let del = r.len() as u64;
    let ins = ins.as_bytes();
    let n_ins = ins.len() as u64;
    // The first coding index the change touches, in transcription order.
    let i0 = if del > 0 {
        let (first, last) = if reverse(t) { (pos + del - 1, pos) } else { (pos, pos + del - 1) };
        // All of it on one stretch of exon, starting in the coding sequence.
        let i = coding_index(t, first)?;
        (tx_offset(t, first)?.abs_diff(tx_offset(t, last)?) == del - 1).then_some(i)?
    } else {
        // An insertion between genomic pos - 1 and pos.
        let (before, after) = if reverse(t) { (pos, pos.checked_sub(1)?) } else { (pos.checked_sub(1)?, pos) };
        coding_index(t, after)?;
        coding_index(t, before)? + 1
    };
    // Bases from the start of the coding sequence, running on into the 3' UTR.
    let reference = |i: u64| -> Option<u8> {
        let b = base(tx_pos(t, cds_start + i)?)?.to_ascii_uppercase();
        Some(if reverse(t) { complement(b) } else { b })
    };
    // Where base `i` of the changed sequence comes from in the reference (none if inserted).
    let src = |i: u64| -> Option<u64> {
        if i < i0 {
            Some(i)
        } else if i < i0 + n_ins {
            None
        } else {
            Some(i - n_ins + del)
        }
    };
    let changed = |i: u64| -> Option<u8> {
        match src(i) {
            Some(j) => reference(j),
            None => Some(ins[(i - i0) as usize].to_ascii_uppercase()),
        }
    };
    let ref_aa = |k: u64| -> Option<u8> { Some(amino_acid(t, k, [reference(3 * k)?, reference(3 * k + 1)?, reference(3 * k + 2)?])) };
    let alt_aa = |k: u64| -> Option<u8> {
        let codon = [changed(3 * k)?, changed(3 * k + 1)?, changed(3 * k + 2)?];
        // A reference codon carried over whole keeps its reading (selenocysteine, say).
        Some(match (src(3 * k), src(3 * k + 1), src(3 * k + 2)) {
            (Some(x), Some(y), Some(z)) if x % 3 == 0 && y == x + 1 && z == x + 2 => amino_acid(t, x / 3, codon),
            _ => translate_with(t.translation.table, codon),
        })
    };
    let name = |a: u8, k: u64| format!("{}{}", three_letter(a), k + 1);
    let k0 = i0 / 3;
    let plain = |p: String| Some(Protein { p, effect: None, warning: None });
    if k0 == 0 {
        return plain("p.?".into());
    }
    let codons = codon_count(t);
    let uncertain = |k: u64| Some(Protein { p: "p.?".into(), effect: Some(UNCERTAIN), warning: Some(internal_stop_warning(k)) });
    // The stop at codon `k` read on as `y`: an extension to the next stop in frame.
    let extension = |k: u64, y: u8| {
        let ter = (k + 1..k + FRAMESHIFT_CODONS).find_map(|j| match alt_aa(j) {
            Some(b'*') => Some(Some(j - k)),
            None => Some(None),
            _ => None,
        });
        let ter = ter.flatten().map_or_else(|| "?".to_owned(), |n| n.to_string());
        Some(Protein { p: format!("p.{}{}extTer{ter}", name(b'*', k), three_letter(y)), effect: Some("stop lost"), warning: None })
    };
    // The protein ends at its stop codon `k` as before: nothing changed.
    let retained = |k: u64| Some(Protein { p: format!("p.{}=", name(b'*', k)), effect: Some(STOP_RETAINED), warning: None });
    let net = n_ins as i64 - del as i64;

    if net % 3 != 0 {
        // Frameshift: the first amino acid that differs, then how far to the new stop.
        let mut k = k0;
        let (x, y) = loop {
            let (x, y) = (ref_aa(k)?, alt_aa(k)?);
            if x != y {
                break (x, y);
            }
            if x == b'*' {
                // Read through to the stop unchanged: the shift comes after it.
                return if k + 1 < codons { uncertain(k) } else { retained(k) };
            }
            if k + 1 >= codons {
                return None;
            }
            k += 1;
        };
        if x == b'X' || y == b'X' {
            return None;
        }
        if x == b'*' {
            return if k + 1 < codons { uncertain(k) } else { extension(k, y) };
        }
        if y == b'*' {
            return Some(Protein { p: format!("p.{}Ter", name(x, k)), effect: Some("stop gained"), warning: None });
        }
        let ter = (k + 1..k + FRAMESHIFT_CODONS).find_map(|j| match alt_aa(j) {
            Some(b'*') => Some(Some(j - k + 1)),
            None => Some(None),
            _ => None,
        });
        let ter = ter.flatten().map_or_else(|| "?".to_owned(), |n| n.to_string());
        return plain(format!("p.{}{}fsTer{ter}", name(x, k), three_letter(y)));
    }

    // In frame: translate a window around the change and trim what's shared at each end; the
    // prefix goes first, which leaves the change at its most 3' (C-terminal) position.
    let last = (i0 + del.max(1) - 1) / 3;
    let ks = k0.saturating_sub(10);
    let ke = (last + 10).min(codons.saturating_sub(1));
    let ref_win: Vec<u8> = (ks..=ke).map(ref_aa).collect::<Option<_>>()?;
    let alt_end = (ke as i64 + net / 3).max(ks as i64) as u64;
    let alt_win: Vec<u8> = (ks..=alt_end).map(alt_aa).collect::<Option<_>>()?;
    if ref_win.contains(&b'X') || alt_win.contains(&b'X') {
        return plain("p.?".into());
    }
    let pre = ref_win.iter().zip(&alt_win).take_while(|(a, b)| a == b).count();
    // Bases added or lost after the stop codon, which is read as before: no protein change.
    if net != 0
        && let Some(j) = ref_win[..pre].iter().position(|&a| a == b'*')
        && ks + j as u64 + 1 >= codons
    {
        return retained(ks + j as u64);
    }
    let room = ref_win.len().min(alt_win.len()) - pre;
    let suf = ref_win.iter().rev().zip(alt_win.iter().rev()).take(room).take_while(|(a, b)| a == b).count();
    let gone = &ref_win[pre..ref_win.len() - suf];
    let mut added = &alt_win[pre..alt_win.len() - suf];
    if let Some(j) = gone.iter().position(|&a| a == b'*') {
        let k = ks + (pre + j) as u64;
        if k + 1 < codons {
            return uncertain(k);
        }
        if added.contains(&b'*') {
            return None;
        }
        // The stop codon itself is the first change: read on from what replaces it.
        if j == 0
            && let Some(y) = alt_aa(k).filter(|&y| y != b'*' && y != b'X')
        {
            return extension(k, y);
        }
        return Some(Protein { p: "p.?".into(), effect: Some("stop lost"), warning: None });
    }
    let at = |j: usize| name(ref_win[j], ks + j as u64);
    let span = |a: usize, b: usize| if a == b { at(a) } else { format!("{}_{}", at(a), at(b)) };
    let letters = |s: &[u8]| s.iter().map(|&a| three_letter(a)).collect::<String>();
    // Translation ends at the first stop: nothing after it is named (`p.Cys114Ter`, not
    // `p.Cys114delinsTerPro`), and a residue that becomes it is a nonsense substitution.
    if let Some(s) = added.iter().position(|&a| a == b'*') {
        added = &added[..=s];
        if s == 0 && !gone.is_empty() {
            return Some(Protein { p: format!("p.{}Ter", at(pre)), effect: Some("stop gained"), warning: None });
        }
    }
    let mut effect = added.contains(&b'*').then_some("stop gained");
    let p = match (gone.len(), added.len()) {
        // Silent: named by the codons it changes, as for an SNV (`p.=` is the whole protein).
        (0, 0) => {
            effect = Some("synonymous");
            format!("p.{}=", span((k0 - ks) as usize, (last.min(ke) - ks) as usize))
        }
        // One amino acid for another is a substitution, whatever the bases (an MNV such as
        // KRAS GGT>GAC is p.Gly12Asp, missense): `delins` is for more than one.
        (1, 1) => {
            effect = effect.or(Some("missense"));
            format!("p.{}{}", at(pre), three_letter(added[0]))
        }
        (_, 0) => format!("p.{}del", span(pre, pre + gone.len() - 1)),
        (0, n) if pre >= n && ref_win[pre - n..pre] == *added => format!("p.{}dup", span(pre - n, pre - 1)),
        (0, _) if pre >= 1 && pre < ref_win.len() => format!("p.{}_{}ins{}", at(pre - 1), at(pre), letters(added)),
        (0, _) => return None,
        (n, _) => format!("p.{}delins{}", span(pre, pre + n - 1), letters(added)),
    };
    Some(Protein { p, effect, warning: None })
}

/// The genomic range of the exon or intron a change at `pos` (`len` bases, or 0 for an
/// insertion before `pos`) is in, or anywhere when it lies wholly outside the transcript.
/// `None` when it straddles a boundary.
fn shift_bounds(t: &Transcript, pos: u64, len: u64) -> Option<(u64, u64)> {
    // For an insertion, the base before it on the transcript.
    let (s, e) = match (len, reverse(t)) {
        (0, true) => (pos, pos + 1),
        (0, false) => (pos.checked_sub(1)?, pos),
        _ => (pos, pos + len),
    };
    if e <= t.start || s >= t.end {
        return Some((0, u64::MAX));
    }
    let mut prev: Option<u64> = None;
    for &(xs, xe) in &t.exons {
        if let Some(p) = prev
            && s >= p
            && e <= xs
        {
            return Some((p, xs));
        }
        if s >= xs && e <= xe {
            return Some((xs, xe));
        }
        prev = Some(xe);
    }
    None
}

/// A pure insertion or deletion placed by [`shift_3prime`].
struct Shifted {
    /// Where it is named: its most 3' place in the exon or intron it starts in.
    pos: u64,
    r: String,
    a: String,
    /// Its most 5' place on the transcript (the same change): where it "starts".
    first: u64,
    /// The repeat it is in runs on over a junction into the next exon or intron (3'): the
    /// same change could be read there too. `Some(true)` when that next region is an exon.
    across: Option<bool>,
}

/// Places a pure deletion (`alt` empty) or insertion (`reference` empty) at `pos`, whatever
/// padding the VCF gave it: all the places the repeat it sits in allows are found on the
/// genome; the change starts at the most 5' of them on the transcript (left on +, right on −)
/// and is then moved to its most 3' place (HGVS's 3' rule) without leaving that exon or intron
/// (as Mutalyzer and biocommons hgvs do by default): a repeat running over an exon-intron
/// junction doesn't carry an exonic change into the splice site, or an intronic one into the
/// exon. Other changes, or unknown bases, stay.
fn shift_3prime(t: &Transcript, mut pos: u64, r: String, a: String, base: &impl Fn(u64) -> Option<u8>) -> Shifted {
    const LIMIT: u64 = 1000;
    let up = |p: u64| base(p).map(|b| b.to_ascii_uppercase());
    if r.is_empty() == a.is_empty() {
        return Shifted { pos, r, a, first: pos, across: None };
    }
    let deleting = a.is_empty();
    let mut d: Vec<u8> = if deleting { r.into_bytes() } else { a.into_bytes() };
    let len = if deleting { d.len() as u64 } else { 0 };
    // The repeat's left end: one starting point for every padding of the same change.
    for _ in 0..LIMIT {
        if pos.checked_sub(1).and_then(up) != Some(*d.last().expect("non-empty")) {
            break;
        }
        d.rotate_right(1);
        pos -= 1;
    }
    let lo = pos;
    let (mut hi, n) = (lo, d.len() as u64);
    while hi - lo < LIMIT && up(hi + len) == Some(d[((hi - lo) % n) as usize]) {
        hi += 1;
    }
    // Every place in lo..=hi is the same change. Its most 5' place on the transcript, kept on
    // it if any place is (a deletion wholly on it; an insertion beside a transcript base, the
    // one before it on the transcript).
    let first = if reverse(t) {
        let last = if deleting { t.end.saturating_sub(len) } else { t.end.saturating_sub(1) };
        hi.min(last).max(lo)
    } else {
        let start = if deleting { t.start } else { t.start + 1 };
        lo.max(start).min(hi)
    };
    // 3' (up the genome on +, down it on −) from there, within its exon or intron. An insertion
    // whose most 5' place is right at a junction is between the two: it goes with the one 3'
    // of it, where the 3' rule takes it (an intron's first base duplicated is c.N+1dup).
    let region = match shift_bounds(t, first, len) {
        Some((rs, _)) if !deleting && reverse(t) && first == rs && rs != t.start => shift_bounds(t, first - 1, 0),
        Some((_, re)) if !deleting && !reverse(t) && first == re && re != t.end => shift_bounds(t, first + 1, 0),
        region => region,
    };
    let (to, across) = match region {
        // Already across a junction (a deletion taking both sides): it stays.
        None => (first, None),
        Some((rs, _)) if reverse(t) => {
            // An insertion keeps a transcript base on its 3' side.
            let floor = if deleting { rs } else { rs.max(t.start + 1) };
            let next = (rs != t.start && rs > 0).then(|| locate(t, rs - 1)).flatten();
            (lo.max(floor).min(first), next.filter(|_| lo < rs).map(|l| !l.feature.starts_with("intron")))
        }
        Some((_, re)) => {
            let ceil = if deleting { re.saturating_sub(len) } else { re };
            let next = (re != t.end).then(|| locate(t, re)).flatten();
            let beyond = if deleting { hi + len > re } else { hi > re };
            (hi.min(ceil).max(first), next.filter(|_| beyond).map(|l| !l.feature.starts_with("intron")))
        }
    };
    d.rotate_left(((to - lo) % n) as usize);
    let d = String::from_utf8(d).unwrap_or_default();
    let (r, a) = if deleting { (d, String::new()) } else { (String::new(), d) };
    Shifted { pos: to, r, a, first, across }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A forward transcript: exons [100,110) [200,215) [300,320); CDS 105..305.
    /// Reverse transcript: the same exons on the minus strand.
    fn tx(strand: Strand) -> Transcript {
        Transcript {
            contig: 0,
            start: 100,
            end: 320,
            strand,
            id: "T1".into(),
            gene_id: "G1".into(),
            gene_name: "GENE".into(),
            biotype: "protein_coding".into(),
            exons: vec![(100, 110), (200, 215), (300, 320)],
            cds: Some((105, 305)),
            tags: 0,
            representative: true,
            translation: Default::default(),
            exons_from_span: false,
        }
    }

    /// Reference: ATG at 105..108, then filler; the genome is "N" elsewhere.
    fn genome(p: u64) -> Option<u8> {
        // Coding (forward): 105..110 = ATGGC, 200..215 = AGAAATTTGGGCCCT, 300..305 = TAAGG.
        let s: &[(u64, &[u8])] = &[(100, b"CCCCCATGGC"), (200, b"AGAAATTTGGGCCCT"), (300, b"TAAGGCCCCCCCCCCCCCCC")];
        s.iter().find_map(|(st, seq)| (p >= *st && p < st + seq.len() as u64).then(|| seq[(p - st) as usize]))
    }

    #[test]
    fn genetic_code() {
        assert_eq!(translate(*b"ATG"), b'M');
        assert_eq!(translate(*b"TAA"), b'*');
        assert_eq!(translate(*b"TGG"), b'W');
        assert_eq!(translate(*b"GGC"), b'G');
        assert_eq!(translate(*b"AAA"), b'K');
        assert_eq!(translate(*b"NAA"), b'X');
        assert_eq!(three_letter(b'K'), "Lys");
    }

    #[test]
    fn forward_positions_and_codons() {
        let t = tx(Strand::Forward);
        assert_eq!(coding_index(&t, 105), Some(0));
        assert_eq!(coding_index(&t, 109), Some(4));
        assert_eq!(coding_index(&t, 200), Some(5), "exon 2 continues the coding sequence");
        assert_eq!(codon_positions(&t, 1), Some([108, 109, 200]), "a codon split across an intron");
        assert_eq!(locate(&t, 103).unwrap().c, "-2");
        assert_eq!(locate(&t, 105).unwrap().c, "1");
        assert_eq!(locate(&t, 306).unwrap().c, "*2");
        assert_eq!(locate(&t, 306).unwrap().kind, Region::Utr3);
        let i = locate(&t, 111).unwrap();
        assert_eq!((i.c.as_str(), i.kind, i.feature.as_str()), ("5+2", Region::SpliceDonor, "intron 1/2"));
        let i = locate(&t, 198).unwrap();
        assert_eq!((i.c.as_str(), i.kind), ("6-2", Region::SpliceAcceptor));
        assert_eq!(locate(&t, 150).unwrap().kind, Region::Intron);
        assert_eq!(locate(&t, 115).unwrap().kind, Region::SpliceRegion);
        assert!(locate(&t, 99).is_none());
    }

    #[test]
    fn forward_substitutions() {
        let t = tx(Strand::Forward);
        // ATG -> ACG at the second base: start lost (Met1Thr).
        let c = consequence(&t, 106, "T", "C", genome).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.2T>C", Some("p.Met1Thr"), "start lost"));
        // Codon 2 = G C A (108, 109, 200) = GCA Ala; 200 A>G gives GCG Ala: synonymous.
        let c = consequence(&t, 200, "A", "G", genome).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.6A>G", Some("p.Ala2="), "synonymous"));
        assert_eq!(c.codons, Some(("gcA".into(), "gcG".into())));
        // Codon 3 = 201..204 = GAA Glu; 201 G>T gives TAA: stop gained.
        let c = consequence(&t, 201, "G", "T", genome).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Glu3Ter"), "stop gained"));
        // Codon 4 = 204..207 = ATT Ile; 204 A>C -> CTT Leu: missense.
        let c = consequence(&t, 204, "A", "C", genome).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Ile4Leu"), "missense"));
        // Intronic and UTR changes get a c. position but no protein change.
        let c = consequence(&t, 111, "C", "T", |_| None).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str(), c.hgvs_p), ("c.5+2C>T", "splice donor", None));
    }

    #[test]
    fn reverse_strand() {
        let t = tx(Strand::Reverse);
        // Transcription starts at 319; CDS starts at 304 (the 5' end of the minus-strand CDS).
        assert_eq!(coding_index(&t, 304), Some(0));
        assert_eq!(coding_index(&t, 300), Some(4));
        assert_eq!(coding_index(&t, 214), Some(5));
        assert_eq!(codon_positions(&t, 1), Some([301, 300, 214]));
        assert_eq!(locate(&t, 310).unwrap().c, "-6");
        // Codon 1 on the minus strand: complement of 304, 303, 302 = rc(TAAGG)[..3]:
        // 304 G->C, 303 G->C, 302 A->T: "CCT" Pro. 303 G>A on + is C>T on -: CTT Leu.
        let c = consequence(&t, 303, "G", "A", genome).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.2C>T", Some("p.Pro1Leu"), "missense"));
        // Deleting a G from the GG at 303..305 on −: 3' on the transcript is the lower one.
        let c = consequence(&t, 303, "GG", "G", genome).unwrap();
        assert_eq!(c.hgvs_c, "c.2del");
        // Intron between exon 1 (300..320) and exon 2 (200..215), 2 bases after exon 1's end.
        let i = locate(&t, 298).unwrap();
        assert_eq!((i.c.as_str(), i.kind), ("5+2", Region::SpliceDonor));
    }

    #[test]
    fn indels() {
        let t = tx(Strand::Forward);
        // Delete one A of the AAA run at 202..205: HGVS names the most 3' one (c.10).
        let c = consequence(&t, 201, "GA", "G", genome).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10del", "frameshift"));
        // Without the reference it can't shift, and says where the VCF put it.
        assert_eq!(consequence(&t, 201, "GA", "G", |_| None).unwrap().hgvs_c, "c.8del");
        // Inserting an A into the run shifts to after its last A, where it repeats that A: a dup.
        let c = consequence(&t, 201, "G", "GA", genome).unwrap();
        assert_eq!(c.hgvs_c, "c.10dup");
        // Delete 3 coding bases: inframe.
        let c = consequence(&t, 201, "GAAA", "G", genome).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.8_10del", "inframe deletion"));
        // Insert 2 bases after 201: frameshift.
        let c = consequence(&t, 201, "G", "GTT", genome).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.7_8insTT", "frameshift"));
        // A frameshift that never meets a stop on this transcript.
        let c = consequence(&t, 201, "GA", "G", genome).unwrap();
        assert_eq!(c.hgvs_p.as_deref(), Some("p.Ile4PhefsTer?"));
        assert_eq!(consequence(&t, 201, "GAAA", "G", genome).unwrap().hgvs_p.as_deref(), Some("p.Glu3_Ile4delinsVal"));
        // Symbolic alleles and breakends are not sequence: no consequence.
        for alt in ["<DEL>", "<DUP:TANDEM>", "G]chr2:100]", ".", "*"] {
            assert!(consequence(&t, 201, "G", alt, genome).is_none(), "{alt}");
        }
    }

    #[test]
    fn indel_protein_changes() {
        // One exon: ATG AAA TTC TTT GGC TGA (Met Lys Phe Phe Gly Ter), then GTAAAAAAAAAA.
        const SEQ: &[u8] = b"ATGAAATTCTTTGGCTGAGTAAAAAAAAAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let fwd = Transcript { exons: vec![(0, 30)], cds: Some((0, 18)), start: 0, end: 30, ..tx(Strand::Forward) };
        let p = |t: &Transcript, pos, r, a| consequence(t, pos, r, a, seq).unwrap();
        // Losing one Phe of the two (like CFTR F508del): the 3'-most is named.
        let c = p(&fwd, 8, "CTTT", "C");
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.10_12del", Some("p.Phe4del"), "inframe deletion"));
        // A third Phe is a duplication of the second, in c. as in p.
        let c = p(&fwd, 11, "T", "TTTT");
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref()), ("c.10_12dup", Some("p.Phe4dup")));
        // Two Lys inserted between Phe4 and Gly5.
        assert_eq!(p(&fwd, 11, "T", "TAAAAAA").hgvs_p.as_deref(), Some("p.Phe4_Gly5insLysLys"));
        // In-frame insertion of a stop: a stop gained.
        let c = p(&fwd, 11, "T", "TTAA");
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Phe4_Gly5insTer"), "stop gained"));
        // Frameshift: GGC -> GCT Ala, then GAG Glu, TAA stop: the third codon from the change.
        assert_eq!(p(&fwd, 11, "TG", "T").hgvs_p.as_deref(), Some("p.Gly5AlafsTer3"));
        // A frameshift whose first new codon is a stop: TA before GGC makes TAG.
        assert_eq!(p(&fwd, 11, "T", "TTA").hgvs_p.as_deref(), Some("p.Gly5Ter"));
        // An inserted T shifts to the end of the TT run; Phe3 is unchanged (TTC -> TTT).
        assert_eq!(p(&fwd, 5, "A", "AT").hgvs_p.as_deref(), Some("p.Phe4LeufsTer?"));
        // Anything touching the start codon.
        assert_eq!(p(&fwd, 0, "ATG", "A").hgvs_p.as_deref(), Some("p.?"));
        // The same changes read on the minus strand of the reverse complement.
        let rc: Vec<u8> = SEQ.iter().rev().map(|&b| complement(b)).collect();
        let rseq = |q: u64| rc.get(q as usize).copied();
        let rev = Transcript { exons: vec![(0, 30)], cds: Some((12, 30)), start: 0, end: 30, ..tx(Strand::Reverse) };
        // CTTT on + at 8..12 is AAAG at 18..22 on the reverse complement: delete AAA after G at 17.
        let c = consequence(&rev, 17, "GAAA", "G", rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref()), ("c.10_12del", Some("p.Phe4del")));
        // And the duplication: AAA inserted before the AAA run (3' on the minus strand).
        let anchor = (rc[17] as char).to_string();
        let c = consequence(&rev, 17, &anchor, &format!("{anchor}AAA"), rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref()), ("c.10_12dup", Some("p.Phe4dup")));
        // An insertion that repeats nothing stays an insertion.
        assert_eq!(p(&fwd, 11, "T", "TAAAAAA").hgvs_c, "c.12_13insAAAAAA");
        // Just inside an exon, in a repeat reaching back into the intron (...TTTAG|CTTT, AGC
        // inserted after the C): the change starts in the intron, so it is named there, at
        // its end, and not a dup (the bases before it are TAG).
        const SPLICED: &[u8] = b"ATGAAATTCTGTAAGTTTAGCTTTGGCTGAAAAAAAAAAA";
        let spliced = |q: u64| SPLICED.get(q as usize).copied();
        let two = Transcript { exons: vec![(0, 10), (20, 40)], cds: Some((0, 35)), start: 0, end: 40, ..tx(Strand::Forward) };
        let c = consequence(&two, 20, "C", "CAGC", spliced).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.11-1_11insCAG", "splice acceptor"));
        assert!(c.warnings.iter().any(|w| w.ends_with("next exon")), "{:?}", c.warnings);
    }

    /// One exon over all of `seq`, coding over `cds`, on `strand`.
    fn single(len: u64, cds: Option<(u64, u64)>, strand: Strand) -> Transcript {
        Transcript { exons: vec![(0, len)], cds, start: 0, end: len, ..tx(strand) }
    }

    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| complement(b)).collect()
    }

    #[test]
    fn selenocysteine_and_internal_stops() {
        // ATG TGA AAA TAA: codon 2 is a TGA inside the coding sequence.
        const SEQ: &[u8] = b"ATGTGAAAATAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let mut t = single(12, Some((0, 12)), Strand::Forward);
        // Not marked: the TGA>TAA change is not called synonymous.
        let c = consequence(&t, 4, "G", "A", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.?"), UNCERTAIN));
        assert!(c.warnings.iter().any(|w| w.contains("codon 2 is a stop")), "{:?}", c.warnings);
        // Marked as selenocysteine (a GTF Selenocysteine row): Sec to stop.
        t.translation.selenocysteine = vec![(3, 6)];
        let c = consequence(&t, 4, "G", "A", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Sec2Ter"), "stop gained"));
        assert!(c.warnings.is_empty());
        let c = consequence(&t, 5, "A", "G", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Sec2Trp"), "missense"));
        // Indels read through it: Lys3 deleted, the Sec before it unchanged.
        assert_eq!(consequence(&t, 5, "AAAA", "A", seq).unwrap().hgvs_p.as_deref(), Some("p.Lys3del"));
        assert_eq!(amino_acid(&t, 1, *b"TGA"), b'U');
        assert_eq!(three_letter(b'U'), "Sec");
    }

    #[test]
    fn mitochondrial_code() {
        // ATG TGA TAA: TGA is Trp in the vertebrate mitochondrial code.
        const SEQ: &[u8] = b"ATGTGATAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let mut t = single(9, Some((0, 9)), Strand::Forward);
        assert_eq!(consequence(&t, 4, "G", "A", seq).unwrap().hgvs_p.as_deref(), Some("p.?"), "standard code: internal stop");
        t.translation.table = 2;
        let c = consequence(&t, 4, "G", "A", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Trp2Ter"), "stop gained"));
        assert_eq!(translate_with(2, *b"ATA"), b'M');
        assert_eq!(translate_with(2, *b"AGA"), b'*');
        assert_eq!(translate_with(1, *b"AGA"), b'R');
    }

    #[test]
    fn incomplete_start_uses_the_phase() {
        // Phase 1: the first complete codon is TGG at 1..4, then AAA, TAA.
        const SEQ: &[u8] = b"ATGGAAATAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let mut t = single(10, Some((0, 10)), Strand::Forward);
        t.translation.phase = 1;
        let c = consequence(&t, 2, "G", "A", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Trp1Ter"), "stop gained"));
        assert_eq!(c.codons, Some(("tGg".into(), "tAg".into())));
        assert!(c.warnings.iter().any(|w| w.contains("phase 1")));
        // The base before the first complete codon has no codon.
        let c = consequence(&t, 0, "A", "G", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.?"), "coding (incomplete codon)"));
        assert_eq!(codon_positions(&t, 0), Some([1, 2, 3]));
        // The same on the minus strand.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let mut t = single(10, Some((0, 10)), Strand::Reverse);
        t.translation.phase = 1;
        assert_eq!(codon_positions(&t, 0), Some([8, 7, 6]));
        let c = consequence(&t, 7, "C", "T", rseq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Trp1Ter"), "stop gained"));
    }

    #[test]
    fn changes_past_the_transcript_ends() {
        let t = tx(Strand::Forward);
        // Deleting 99..102 takes the transcript's first two bases.
        let c = consequence(&t, 98, "AAAA", "A", |_| None).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.?_-4del", "transcript boundary (effect unknown)"));
        assert!(!c.warnings.is_empty());
        // And its last two, on the minus strand (its 5' end).
        let t = tx(Strand::Reverse);
        let c = consequence(&t, 317, "AAAA", "A", |_| None).unwrap();
        assert_eq!(c.hgvs_c, "c.?_-14del");
        // Wholly outside: nothing.
        assert!(consequence(&t, 400, "AAAA", "A", |_| None).is_none());
    }

    #[test]
    fn stop_loss_by_longer_changes() {
        // ATG AAA TGA, then CAA TAA in the 3' UTR.
        const SEQ: &[u8] = b"ATGAAATGACAATAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(15, Some((0, 9)), Strand::Forward);
        let c = consequence(&t, 6, "TGA", "CAA", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Ter3GlnextTer2"), "stop lost"));
        // Deleting the stop codon reads on into the UTR.
        let c = consequence(&t, 5, "ATGA", "A", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Ter3GlnextTer1"), "stop lost"));
        // A frameshift in the stop codon: an extension with no stop in sight.
        let c = consequence(&t, 6, "TG", "T", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Ter3TyrextTer?"), "frameshift"));
        // The same MNV on the minus strand.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = single(15, Some((6, 15)), Strand::Reverse);
        let c = consequence(&t, 6, "TCA", "TTG", rseq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Ter3GlnextTer2"), "stop lost"));
    }

    #[test]
    fn indels_stay_in_their_exon_when_shifted() {
        // Exon 0..9 ends in AAA, and the A run goes on past the transcript.
        const SEQ: &[u8] = b"ATGCCCAAAAAAAAACCCCC";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(9, None, Strand::Forward);
        let c = consequence(&t, 6, "AA", "A", seq).unwrap();
        assert_eq!(c.hgvs_c, "n.9del");
        let c = consequence(&t, 6, "A", "AA", seq).unwrap();
        assert_eq!(c.hgvs_c, "n.9dup");
        // Minus strand: the run goes on below the exon.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = Transcript { exons: vec![(11, 20)], start: 11, end: 20, ..single(20, None, Strand::Reverse) };
        assert_eq!(consequence(&t, 12, "TT", "T", rseq).unwrap().hgvs_c, "n.9del");
    }

    #[test]
    fn deletions_into_a_splice_site() {
        // Coding exons 0..9 and 19..28; chrS:7 GTAGTCC>G takes coding 7..9 and intron 9..13.
        const SEQ: &[u8] = b"ATGCCCGTAGTCCCGATACGCTAAATAAGGGGG";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 9), (19, 28)], cds: Some((0, 28)), start: 0, end: 28, ..tx(Strand::Forward) };
        let c = consequence(&t, 6, "GTAGTCC", "G", seq).unwrap();
        assert_eq!((c.effect.as_str(), c.hgvs_p.as_deref()), ("splice donor", None));
        assert!(c.warnings.iter().any(|w| w.contains("exon-intron boundary")));
        // Wholly intronic, over the acceptor from inside: still a splice acceptor.
        assert_eq!(consequence(&t, 14, "CGATA", "C", seq).unwrap().effect, "splice acceptor");
        // The mirror image on the minus strand.
        const REV: &[u8] = b"CCCCCTTATTTAGCGTATCGGGACTACGGGCAT";
        let rseq = |p: u64| REV.get(p as usize).copied();
        let t = Transcript { exons: vec![(5, 14), (24, 33)], cds: Some((5, 33)), start: 5, end: 33, ..tx(Strand::Reverse) };
        let c = consequence(&t, 19, "GGGACTA", "G", rseq).unwrap();
        assert_eq!((c.effect.as_str(), c.hgvs_p.as_deref()), ("splice donor", None));
    }

    #[test]
    fn insertions_beside_the_coding_sequence() {
        const SEQ: &[u8] = b"ATGAAATGACAATAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(15, Some((0, 9)), Strand::Forward);
        // Right after the stop codon: 3' UTR, not a frameshift.
        let c = consequence(&t, 8, "A", "AT", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str(), c.hgvs_p.as_deref()), ("c.9_*1insT", "3' UTR", None));
        // Inside the stop codon: a frameshift.
        assert_eq!(consequence(&t, 7, "G", "GT", seq).unwrap().effect, "frameshift");
        // Minus strand: after the stop is down the genome.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = single(15, Some((6, 15)), Strand::Reverse);
        let c = consequence(&t, 5, "G", "GA", rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.9_*1insT", "3' UTR"));
        // At an exon's edge inside the coding sequence: the junction, flagged.
        let t = tx(Strand::Forward);
        let c = consequence(&t, 109, "C", "CTT", |_| None).unwrap();
        assert_eq!(c.effect, "splice donor");
        assert!(!c.warnings.is_empty());
    }

    #[test]
    fn padded_alleles_are_named_as_their_minimal_change() {
        // ATG AAA TTC TTT GGC TGA, as in `indel_protein_changes`.
        const SEQ: &[u8] = b"ATGAAATTCTTTGGCTGAGTAAAAAAAAAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(30, Some((0, 18)), Strand::Forward);
        let p = |pos, r, a| consequence(&t, pos, r, a, seq).unwrap();
        // Shared bases after the change too (as in a split multi-allelic record): one name.
        let c = p(8, "CTTTG", "CG");
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.10_12del", Some("p.Phe4del"), "inframe deletion"));
        assert_eq!(p(11, "TG", "TTTTG").hgvs_c, "c.10_12dup");
        assert_eq!(p(3, "AAA", "AA").hgvs_c, "c.6del");
        // A padded SNV is an SNV: AAA>AGA is c.5A>G, Lys2Arg.
        let c = p(3, "AAA", "AGA");
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.5A>G", Some("p.Lys2Arg"), "missense"));
        assert_eq!(c.codons, Some(("aAa".into(), "aGa".into())));
    }

    #[test]
    fn a_multi_base_change_of_one_amino_acid_is_missense() {
        const SEQ: &[u8] = b"ATGAAATTCTTTGGCTGAGTAAAAAAAAAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(30, Some((0, 18)), Strand::Forward);
        // AAA>GTA: Lys2Val, an MNV (like KRAS GGT>GAC, p.Gly12Asp).
        let c = consequence(&t, 3, "AA", "GT", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.4_5delinsGT", Some("p.Lys2Val"), "missense"));
        // To a stop: stop gained.
        let c = consequence(&t, 3, "AAA", "TAG", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Lys2Ter"), "stop gained"));
        // Two amino acids changed stay a delins; a silent MNV is synonymous.
        assert_eq!(consequence(&t, 3, "AAATT", "GTAAA", seq).unwrap().hgvs_p.as_deref(), Some("p.Lys2_Phe3delinsValAsn"));
        let c = consequence(&t, 8, "CTTT", "TTTC", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Phe3_Phe4="), "synonymous"));
        // The same MNV on the minus strand: + 3..5 AA>GT is 25..27 TT>AC on the reverse complement.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = single(30, Some((12, 30)), Strand::Reverse);
        let c = consequence(&t, 25, "TT", "AC", rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.4_5delinsGT", Some("p.Lys2Val"), "missense"));
    }

    #[test]
    fn intronic_insertions_repeating_their_neighbours_are_dups() {
        // Exons 0..10 and 30..40, coding throughout; the intron has TTTTT at 15..20, then CAGCAG.
        const SEQ: &[u8] = b"ATGAAACCCGGTAAGTTTTTCAGCAGCTAGATGAAACCCG";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 10), (30, 40)], cds: Some((0, 40)), start: 0, end: 40, ..tx(Strand::Forward) };
        // A T into the run shifts to its end and repeats the last T: c.10+10dup, not ins.
        let c = consequence(&t, 14, "G", "GT", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10+10dup", "intron"));
        // A CAG unit, shifted 3' through CAGCAG.
        assert_eq!(consequence(&t, 19, "T", "TCAG", seq).unwrap().hgvs_c, "c.11-6_11-4dup");
        // An insertion repeating nothing stays one.
        assert_eq!(consequence(&t, 14, "G", "GA", seq).unwrap().hgvs_c, "c.10+5_10+6insA");
        // The minus strand: the same T insertion read on the reverse complement.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = Transcript { strand: Strand::Reverse, ..t };
        assert_eq!(consequence(&t, 24, "A", "AA", rseq).unwrap().hgvs_c, "c.10+10dup");
    }

    #[test]
    fn any_change_to_the_start_codon_is_start_lost() {
        // 5' UTR CC, then ATG AAA TTC TTT GGC TGA.
        const SEQ: &[u8] = b"CCATGAAATTCTTTGGCTGAGTAAAAAAAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(30, Some((2, 20)), Strand::Forward);
        let c = |pos, r, a| {
            let c = consequence(&t, pos, r, a, seq).unwrap();
            (c.hgvs_c, c.hgvs_p, c.effect)
        };
        let lost = |hgvs: &str| (hgvs.to_owned(), Some("p.?".to_owned()), "start lost".to_owned());
        assert_eq!(c(2, "AT", "GC"), lost("c.1_2delinsGC"), "an MNV");
        assert_eq!(c(2, "ATGA", "A"), lost("c.2_4del"), "an in-frame deletion");
        assert_eq!(c(2, "AT", "A"), lost("c.2del"), "a frameshift");
        assert_eq!(c(3, "T", "TCC"), lost("c.2_3insCC"), "an insertion inside it");
        assert_eq!(c(1, "CA", "TC"), lost("c.-1_1delinsTC"), "from the 5' UTR");
        // Beside it, not in it: not a start loss.
        assert_eq!(c(4, "G", "GCC").2, "frameshift");
        assert_eq!(c(1, "C", "CTT").2, "5' UTR");
        // The minus strand: + 2..4 AT>GC is GC>AT... read on the reverse complement.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = single(30, Some((10, 28)), Strand::Reverse);
        let c = consequence(&t, 26, "AT", "GC", rseq).unwrap();
        assert_eq!((c.hgvs_c, c.hgvs_p, c.effect), lost("c.1_2delinsGC"));
        let c = consequence(&t, 27, "TG", "GA", rseq).unwrap();
        assert_eq!((c.hgvs_c, c.hgvs_p, c.effect), lost("c.-1_1delinsTC"));
    }

    #[test]
    fn an_indel_in_a_repeat_across_a_junction_stays_in_its_exon_whatever_the_padding() {
        // Like EGFR c.240: exon 1 ends ...AAAG, intron 1 starts GTTGG; exon 2 at 20..40.
        const SEQ: &[u8] = b"ATGCCCAAAGGTTGGCCCAGCCCAAATTTGGGCCCTAAAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 10), (20, 40)], cds: Some((0, 37)), start: 0, end: 40, ..tx(Strand::Forward) };
        let across = |c: &Consequence, there: &str| c.warnings.iter().any(|w| w.contains("across an exon-intron junction") && w.ends_with(there));
        // One G of GG deleted: left-aligned, with the other anchor, or padded. It starts in the
        // exon and isn't shifted out of it: the exon's last base, a frameshift.
        for (pos, r, a) in [(8, "AG", "A"), (9, "GG", "G"), (8, "AGG", "AG")] {
            let c = consequence(&t, pos, r, a, seq).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10del", "frameshift"), "{r}>{a}");
            assert!(c.hgvs_p.is_some() && across(&c, "next intron"), "{c:?}");
        }
        // One G more: a dup of the exon's G, read in the exon.
        for (pos, r, a) in [(8, "A", "AG"), (9, "G", "GG"), (8, "AG", "AGG")] {
            let c = consequence(&t, pos, r, a, seq).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10dup", "frameshift"), "{r}>{a}");
            assert!(c.hgvs_p.is_some() && across(&c, "next intron"), "{c:?}");
        }
        // GG into the exon's last G (EGFR A>AGG): an exonic insertion at the exon's end.
        let c = consequence(&t, 8, "A", "AGG", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10_10+1insGG", "frameshift"));
        assert!(c.hgvs_p.is_some(), "{c:?}");
        // The mirror image on the minus strand (like KRAS c.450): the same.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 20), (30, 40)], cds: Some((3, 40)), start: 0, end: 40, ..tx(Strand::Reverse) };
        let (x, y) = (r[28] as char, r[29] as char);
        for (pos, r, a) in [(29, "CC".to_owned(), "C".to_owned()), (28, format!("{x}{y}C"), format!("{x}{y}")), (28, format!("{x}C"), x.to_string())] {
            let c = consequence(&t, pos, &r, &a, rseq).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10del", "frameshift"), "{pos} {r}>{a}");
            assert!(across(&c, "next intron"));
        }
        for (pos, r, a) in [(28, x.to_string(), format!("{x}C")), (29, "C".to_owned(), "CC".to_owned())] {
            let c = consequence(&t, pos, &r, &a, rseq).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10dup", "frameshift"), "{pos} {r}>{a}");
        }
        // An exon's last base repeated, the intron's GT after it: a dup in the exon, read there
        // (a frameshift), with no repeat over the junction to warn of.
        const LAST: &[u8] = b"ATGCCCAAACGTAAGCCCAGCCCAAATTTGGGCCCTAAAA";
        let seq = |p: u64| LAST.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 10), (20, 40)], cds: Some((0, 37)), start: 0, end: 40, ..tx(Strand::Forward) };
        let c = consequence(&t, 9, "C", "CC", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10dup", "frameshift"));
        assert!(c.hgvs_p.is_some() && c.warnings.is_empty(), "{c:?}");
    }

    #[test]
    fn an_in_frame_deletion_at_an_exon_end_stays_exonic() {
        // Like NM_003288.4 c.474_476del: exon 1 ends ...CGAG, intron 1 starts GAGT.
        const SEQ: &[u8] = b"ATGAAACGAGGAGTAAGCCCAGCCCAAATTTGGGCCCTAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 10), (21, 40)], cds: Some((0, 40)), start: 0, end: 40, ..tx(Strand::Forward) };
        for (pos, r, a) in [(6, "CGAG", "C"), (9, "GGAG", "G"), (6, "CGAGG", "CG")] {
            let c = consequence(&t, pos, r, a, seq).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.8_10del", Some("p.Gly4del"), "inframe deletion"), "{r}>{a}");
        }
    }

    #[test]
    fn a_dup_in_an_intron_beside_the_exon_is_named_in_the_intron() {
        // Like ATM c.6976-4_6976-1dup: intron 1 (10..20) ends ...TACAG, exon 2 starts AACAAT.
        const SEQ: &[u8] = b"ATGAAACCCGGTAAGTACAGAACAATGGCCCTAGAAAAAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 10), (20, 40)], cds: Some((0, 34)), start: 0, end: 40, ..tx(Strand::Forward) };
        for (r, a) in [("T", "TACAG"), ("TAC", "TACAGAC")] {
            let c = consequence(&t, 15, r, a, seq).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.11-4_11-1dup", None, "splice acceptor"), "{r}>{a}");
            assert!(c.warnings.iter().any(|w| w.contains("named in the intron") && w.ends_with("next exon")), "{:?}", c.warnings);
        }
        // The minus strand, on the reverse complement.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 20), (30, 40)], cds: Some((6, 40)), start: 0, end: 40, ..tx(Strand::Reverse) };
        let anchor = (r[23] as char).to_string();
        let c = consequence(&t, 23, &anchor, &format!("{anchor}CTGT"), rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.11-4_11-1dup", "splice acceptor"));
        // A G inserted right after the exon, before the intron's GG: it can't go 5' into the
        // exon, so it is the intron's (c.10+2dup), not an insertion at the junction.
        const DONOR: &[u8] = b"ATGAAACCCCGGTAAGTTTTAACAATGGCCCTAGAAAAAA";
        let seq = |p: u64| DONOR.get(p as usize).copied();
        let t = Transcript { exons: vec![(0, 10), (20, 40)], cds: Some((0, 34)), start: 0, end: 40, ..tx(Strand::Forward) };
        let c = consequence(&t, 9, "C", "CG", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.10+2dup", "splice donor"));
        assert!(c.warnings.is_empty(), "{:?}", c.warnings);
    }

    #[test]
    fn a_dup_of_the_stop_codon_after_it_changes_no_protein() {
        // ATG AAA TGA, then the 3' UTR (like DUSP15 c.706_708dup, DIDO1 c.1689dup, TRIB3
        // c.1077_*2dup): the copy goes after the stop, which still ends the protein.
        let check = |seq: &[u8], pos: u64, r: &str, a: &str, name: &str| {
            let t = single(seq.len() as u64, Some((0, 9)), Strand::Forward);
            let c = consequence(&t, pos, r, a, |p| seq.get(p as usize).copied()).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), (name, None, "3' UTR"), "+ {r}>{a}");
            // The same on the minus strand: REF/ALT reverse-complemented, anchored on the left.
            let rs = rc(seq);
            let n = seq.len() as u64;
            let t = single(n, Some((n - 9, n)), Strand::Reverse);
            let rpos = n - pos - r.len() as u64;
            let (rr, ra) = (String::from_utf8(rc(r.as_bytes())).unwrap(), String::from_utf8(rc(a.as_bytes())).unwrap());
            let anchor = rs[rpos as usize - 1] as char;
            let (rr, ra) = (format!("{anchor}{}", &rr[..rr.len() - 1]), format!("{anchor}{}", &ra[..ra.len() - 1]));
            let c = consequence(&t, rpos - 1, &rr, &ra, |p| rs.get(p as usize).copied()).unwrap();
            assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), (name, None, "3' UTR"), "- {rr}>{ra}");
        };
        check(b"ATGAAATGACCCCC", 8, "A", "ATGA", "c.7_9dup");
        check(b"ATGAAATAACCCCC", 8, "A", "AA", "c.9dup");
        check(b"ATGAAATAGGACCCCC", 10, "A", "AGGA", "c.9_*2dup");
        check(b"ATGAAATAAGCCCCC", 9, "G", "GAG", "c.9_*1dup");
    }

    #[test]
    fn a_dup_beside_the_start_codon_keeps_its_protein_change() {
        // 5' UTR CC, then ATG GGT GCG AAA TAA (like NM_018478.3 c.3_5dup p.Gly2dup).
        const SEQ: &[u8] = b"CCATGGGTGCGAAATAACCC";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(20, Some((2, 17)), Strand::Forward);
        let c = consequence(&t, 4, "G", "GGGG", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.3_5dup", Some("p.Gly2dup"), "inframe insertion"));
        // ATG AAG CTC ... on the minus strand (like NM_001011878.3): c.2_4dup copies TGA after
        // c.4, ATG ATG AAG: Met1 duplicated; c.3_6dup is a frameshift from Leu3.
        const PLUS: &[u8] = b"CCATGAAGCTCGGCTGGTAGCCC";
        let r = rc(PLUS);
        let rseq = |p: u64| r.get(p as usize).copied();
        let n = PLUS.len() as u64;
        let t = single(n, Some((3, n - 2)), Strand::Reverse);
        // TGA inserted after + 5 (c.4) is TCA after n - 7 on the reverse complement.
        let at = |i: u64| (r[i as usize] as char).to_string();
        let c = consequence(&t, n - 7, &at(n - 7), &format!("{}TCA", at(n - 7)), rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.2_4dup", Some("p.Met1dup"), "inframe insertion"));
        let c = consequence(&t, n - 9, &at(n - 9), &format!("{}CTTC", at(n - 9)), rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.effect.as_str()), ("c.3_6dup", "frameshift"));
        assert!(c.hgvs_p.as_deref().is_some_and(|p| p.starts_with("p.Leu3Glufs")), "{c:?}");
    }

    #[test]
    fn an_indel_in_a_run_from_the_5utr_into_the_atg_keeps_the_start() {
        // Like BRCA2: 5' UTR ...GTAAAA, then ATG: the A run ends at c.1.
        const SEQ: &[u8] = b"GTAAAAATGCCTGGCTAACCC";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(21, Some((6, 18)), Strand::Forward);
        let c = |pos, r, a| {
            let c = consequence(&t, pos, r, a, seq).unwrap();
            (c.hgvs_c, c.hgvs_p, c.effect)
        };
        let utr = |n: &str| (n.to_owned(), None, "5' UTR".to_owned());
        assert_eq!(c(1, "T", "TA"), utr("c.1dup"));
        assert_eq!(c(5, "A", "AA"), utr("c.1dup"));
        assert_eq!(c(1, "TA", "T"), utr("c.1del"));
        assert_eq!(c(5, "AA", "A"), utr("c.1del"));
        assert_eq!(c(1, "TAA", "T"), utr("c.-1_1del"));
        // Into the ATG itself, or taking more than the run: start lost.
        assert_eq!(c(6, "A", "AC").2, "start lost");
        assert_eq!(c(5, "AAT", "A").2, "start lost");
        // KRAS: the minus strand, AAA|ATG on the transcript.
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = single(21, Some((3, 15)), Strand::Reverse);
        for (pos, rr, aa, name) in [(14, "T", "TT", "c.1dup"), (14, "TT", "T", "c.1del"), (17, "T", "TT", "c.1dup")] {
            let c = consequence(&t, pos, rr, aa, rseq).unwrap();
            assert_eq!((c.hgvs_c, c.hgvs_p, c.effect), utr(name), "{pos} {rr}>{aa}");
        }
    }

    #[test]
    fn a_stop_lost_by_an_snv_names_its_extension() {
        // ATG AAA TGA, then CAA TAA in the 3' UTR: TGA>CGA reads on Arg Gln, then stops.
        const SEQ: &[u8] = b"ATGAAATGACAATAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(15, Some((0, 9)), Strand::Forward);
        let c = consequence(&t, 6, "T", "C", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Ter3ArgextTer2"), "stop lost"));
        // As the same change written as an MNV is.
        assert_eq!(consequence(&t, 6, "TG", "CG", seq).unwrap().hgvs_p.as_deref(), Some("p.Ter3ArgextTer2"));
        let r = rc(SEQ);
        let rseq = |p: u64| r.get(p as usize).copied();
        let t = single(15, Some((6, 15)), Strand::Reverse);
        let c = consequence(&t, 8, "A", "G", rseq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.7T>C", Some("p.Ter3ArgextTer2"), "stop lost"));
    }

    #[test]
    fn in_frame_changes_stop_naming_at_the_first_stop() {
        // ATG AAA TGC TTT GGC TGA (Met Lys Cys Phe Gly Ter).
        const SEQ: &[u8] = b"ATGAAATGCTTTGGCTGAGTAAAAAAAAAA";
        let seq = |p: u64| SEQ.get(p as usize).copied();
        let t = single(30, Some((0, 18)), Strand::Forward);
        // ACC into TG|C (like UBE2C c.341_342insACC): TGA CCC, so Cys3 becomes the stop.
        let c = consequence(&t, 7, "G", "GACC", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.8_9insACC", Some("p.Cys3Ter"), "stop gained"));
        // Lys2 replaced by Val, a stop and Pro: nothing after the stop.
        let c = consequence(&t, 3, "AAA", "GTGTAGCCC", seq).unwrap();
        assert_eq!((c.hgvs_p.as_deref(), c.effect.as_str()), (Some("p.Lys2delinsValTer"), "stop gained"));
        // A silent MNV names its residue, as a silent SNV does: CTA>TTG is Leu either way.
        const LEU: &[u8] = b"ATGCTAAAATAA";
        let seq = |p: u64| LEU.get(p as usize).copied();
        let t = single(12, Some((0, 12)), Strand::Forward);
        let c = consequence(&t, 3, "CTA", "TTG", seq).unwrap();
        assert_eq!((c.hgvs_c.as_str(), c.hgvs_p.as_deref(), c.effect.as_str()), ("c.4_6delinsTTG", Some("p.Leu2="), "synonymous"));
    }

    #[test]
    fn ref_mismatch_is_flagged() {
        let t = tx(Strand::Forward);
        // 204 is A; a VCF saying C still gets a consequence, but with a warning.
        let c = consequence(&t, 204, "C", "T", genome).unwrap();
        assert!(c.warnings.iter().any(|w| w.contains("REF C does not match the reference here (A)")), "{:?}", c.warnings);
        assert!(consequence(&t, 204, "A", "T", genome).unwrap().warnings.is_empty());
    }
}
