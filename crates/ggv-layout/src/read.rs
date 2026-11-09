use ggv_io::pasta::pileup::{Read, mods::BaseMods};
use serde::Serialize;

const M: u32 = 0;
const I: u32 = 1;
const D: u32 = 2;
const N: u32 = 3;
const S: u32 = 4;
const EQ: u32 = 7;
const X: u32 = 8;

/// A gapless aligned stretch (CIGAR `M`, `=`, `X`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Block {
    /// 0-based reference start.
    pub ref_start: u64,
    pub query_start: u32,
    pub len: u32,
}

/// A deletion, reference skip, or insertion. For an insertion, `pos` is the reference position
/// the inserted bases come before.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Indel {
    pub pos: u64,
    pub len: u32,
}

/// A read base that differs from the reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Mismatch {
    pub pos: u64,
    /// `A`, `C`, `G`, `T` or `N`.
    pub base: u8,
    pub quality: u8,
}

/// A base-modification call (MM/ML tags) at an aligned reference position.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct ModCall {
    /// Bases from the read's start ([`ModCall::pos`]): eight bytes a call rather than sixteen,
    /// for reads with calls on every A and C.
    pub offset: u32,
    /// The modification code: `m` (5mC), `h` (5hmC), `a` (6mA), … (`?` for a ChEBI number).
    pub code: u8,
    /// The canonical base it modifies (`C`, `A`, …).
    pub canonical: u8,
    /// Likelihood of the modification, 0–255 (ML; 255 when ML is absent).
    pub prob: u8,
}

impl ModCall {
    /// The call's reference position, for a read starting at `start`.
    pub fn pos(&self, start: u64) -> u64 {
        start + u64::from(self.offset)
    }
}

/// One entry of the `SA` tag: another alignment of the same read.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct SplitAlignment {
    pub contig: String,
    /// 0-based.
    pub pos: u64,
    pub reverse: bool,
    pub mapq: u8,
    pub cigar: String,
}

/// Parses an `SA:Z` value: `contig,pos,strand,CIGAR,mapq,NM;` repeated (pos 1-based).
pub fn parse_sa(sa: &str) -> Vec<SplitAlignment> {
    sa.split(';')
        .filter_map(|entry| {
            let f: Vec<&str> = entry.split(',').collect();
            (f.len() >= 5).then_some(())?;
            Some(SplitAlignment {
                contig: f[0].to_owned(),
                pos: f[1].parse::<u64>().ok()?.checked_sub(1)?,
                reverse: f[2] == "-",
                mapq: f[4].parse().ok()?,
                cigar: f[3].to_owned(),
            })
        })
        .collect()
}

/// Where a read's mate is, for paired reads with a mapped mate.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Mate {
    /// Whether the mate maps to the same contig.
    pub same_contig: bool,
    /// The mate's contig, as an index into the BAM header's reference sequences.
    pub contig_id: i32,
    /// The mate's 0-based position.
    pub pos: i64,
    pub reverse: bool,
}

/// One alignment, decoded into what a renderer draws.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct AlignedRead {
    pub name: String,
    /// 0-based start of the aligned part.
    pub start: u64,
    /// 0-based exclusive end of the aligned part.
    pub end: u64,
    pub reverse: bool,
    pub mapq: u8,
    pub flags: u16,
    pub blocks: Vec<Block>,
    pub deletions: Vec<Indel>,
    pub ref_skips: Vec<Indel>,
    pub insertions: Vec<Indel>,
    /// Empty when no reference was given.
    pub mismatches: Vec<Mismatch>,
    pub soft_clip_left: u32,
    pub soft_clip_right: u32,
    /// Soft-clipped bases that differ from the reference where the clip would align (just
    /// before `start` / after `end`); empty without a reference.
    pub clip_mismatches: Vec<Mismatch>,
    /// Haplotype (`HP` tag), e.g. 1 or 2 for phased long reads.
    pub haplotype: Option<u8>,
    /// Read group (`RG` tag).
    pub read_group: Option<Box<str>>,
    /// Template length (`TLEN`): signed insert size for pairs.
    pub template_length: i32,
    pub mate: Option<Mate>,
    /// Other alignments of this read (`SA` tag): split reads and chimeras.
    pub supplementary: Vec<SplitAlignment>,
    /// Base-modification calls on aligned bases, in position order; empty without MM tags.
    /// Where several codes are called at one base (5mC and 5hmC), each is listed.
    pub mods: Vec<ModCall>,
    /// Aligned positions whose base is below the minimum base quality (or has none), in order:
    /// they count as not observed, as in a pileup. Empty unless a minimum was asked for.
    pub low_quality: Vec<u64>,
    /// The value of the tag the view was loaded to show (`ViewOptions::tag`), as text; `None`
    /// when not asked for or the read lacks it.
    pub tag: Option<Box<str>>,
    /// The read's bases (query order, letters), kept only for reads loaded without a reference:
    /// with no mismatches to go by, grouping and listings by base read them here.
    pub bases: Option<Box<[u8]>>,
}

impl AlignedRead {
    /// Whether the base at `pos` was left out for low base quality.
    pub fn is_low_quality(&self, pos: u64) -> bool {
        self.low_quality.binary_search(&pos).is_ok()
    }

    /// The read base aligned at `pos` from [`Self::bases`], when kept: `None` without them, off
    /// the aligned blocks, or for a base the read doesn't store (`=` or no sequence).
    pub fn base_at(&self, pos: u64) -> Option<u8> {
        let bases = self.bases.as_deref()?;
        let b = self.blocks.iter().find(|b| b.ref_start <= pos && pos < b.ref_start + u64::from(b.len))?;
        bases.get((b.query_start as u64 + pos - b.ref_start) as usize).copied().filter(|&c| c != b'=')
    }

    /// Whether the read is paired with a mapped mate.
    pub fn is_paired(&self) -> bool {
        self.flags & 0x1 != 0 && self.flags & 0x8 == 0
    }
}

/// A read's bases as letters, in query order (empty for a read with no stored sequence).
pub fn read_bases(read: &Read) -> Box<[u8]> {
    (0..read.sequence_len()).map(|i| base_letter(read.base(i))).collect()
}

/// Decodes a read's CIGAR into blocks and indels, and finds mismatches against `reference`
/// (its bases and the 0-based position of its first base).
pub fn decode(read: &Read, reference: Option<(&[u8], u64)>) -> AlignedRead {
    decode_with(read, reference, 0)
}

/// [`decode`], treating aligned bases below `min_base_quality` (or without a quality) as not
/// observed: they are listed in [`AlignedRead::low_quality`] and never count as mismatches,
/// exactly as `ggv query pileup` leaves them out of depth and allele counts.
pub fn decode_with(read: &Read, reference: Option<(&[u8], u64)>, min_base_quality: u8) -> AlignedRead {
    decode_parts(read, reference, min_base_quality, true)
}

/// [`decode_with`], leaving out the base-modification calls unless `mods`: reads with calls on
/// every base (ONT 6mA, 5mC and 5hmC on all A and C) hold several times more calls than
/// aligned bases, so a view that doesn't draw them shouldn't keep them.
pub fn decode_parts(read: &Read, reference: Option<(&[u8], u64)>, min_base_quality: u8, mods: bool) -> AlignedRead {
    let start = read.position().max(0) as u64;
    let mut ref_pos = start;
    let mut query_pos: u32 = 0;

    let mut out = AlignedRead {
        name: String::from_utf8_lossy(read.name().strip_suffix(&[0]).unwrap_or(read.name()))
            .into_owned(),
        start,
        end: start,
        reverse: read.is_reverse_complemented(),
        mapq: read.mapping_quality(),
        flags: read.flags(),
        blocks: Vec::new(),
        deletions: Vec::new(),
        ref_skips: Vec::new(),
        insertions: Vec::new(),
        mismatches: Vec::new(),
        soft_clip_left: 0,
        soft_clip_right: 0,
        clip_mismatches: Vec::new(),
        haplotype: ggv_query::haplotype(read).and_then(|v| u8::try_from(v).ok()),
        read_group: read.string_tag(*b"RG").map(|v| String::from_utf8_lossy(v).into()),
        template_length: read.template_length(),
        mate: None,
        supplementary: read.string_tag(*b"SA").map(|v| parse_sa(&String::from_utf8_lossy(v))).unwrap_or_default(),
        mods: Vec::new(),
        low_quality: Vec::new(),
        tag: None,
        bases: None,
    };
    let flags = read.flags();
    if flags & 0x1 != 0 && flags & 0x8 == 0 && read.mate_reference_sequence_id() >= 0 {
        out.mate = Some(Mate {
            same_contig: read.mate_reference_sequence_id() == read.reference_sequence_id(),
            contig_id: read.mate_reference_sequence_id(),
            pos: read.mate_position(),
            reverse: flags & 0x20 != 0,
        });
    }

    let ops = read.cigar();
    let has_seq = read.sequence_len() > 0;
    let quals = read.quality_scores();

    for (i, &c) in ops.iter().enumerate() {
        let (op, len) = (c & 0xf, c >> 4);
        match op {
            M | EQ | X => {
                out.blocks.push(Block { ref_start: ref_pos, query_start: query_pos, len });
                if min_base_quality > 0 {
                    for k in 0..len {
                        let q = quals.get((query_pos + k) as usize).copied().unwrap_or(0xff);
                        if q == 0xff || q < min_base_quality {
                            out.low_quality.push(ref_pos + u64::from(k));
                        }
                    }
                }
                if let (Some((seq, offset)), true) = (reference, has_seq) {
                    for k in 0..len {
                        let rp = ref_pos + u64::from(k);
                        let qp = (query_pos + k) as usize;
                        if min_base_quality > 0 && out.low_quality.binary_search(&rp).is_ok() {
                            continue;
                        }
                        let Some(&ref_base) = rp.checked_sub(offset).and_then(|o| seq.get(o as usize))
                        else {
                            continue;
                        };
                        let base = if qp < read.sequence_len() { base_letter(read.base(qp)) } else { b'N' };
                        if base == b'=' || ref_base == b'N' || base == ref_base {
                            continue;
                        }
                        out.mismatches.push(Mismatch {
                            pos: rp,
                            base,
                            quality: quals.get(qp).copied().unwrap_or(0xff),
                        });
                    }
                }
                ref_pos += u64::from(len);
                query_pos += len;
            }
            I => {
                out.insertions.push(Indel { pos: ref_pos, len });
                query_pos += len;
            }
            D => {
                out.deletions.push(Indel { pos: ref_pos, len });
                ref_pos += u64::from(len);
            }
            N => {
                out.ref_skips.push(Indel { pos: ref_pos, len });
                ref_pos += u64::from(len);
            }
            S => {
                let left = i == 0 || ops[..i].iter().all(|c| c & 0xf == 5);
                if left {
                    out.soft_clip_left = len;
                } else {
                    out.soft_clip_right = len;
                }
                // Where the clipped bases would sit if aligned: just before the start (left
                // clip) or from the current reference position (right clip).
                if let (Some((seq, offset)), true) = (reference, has_seq) {
                    let first_ref = if left { i64::try_from(start).unwrap_or(0) - i64::from(len) } else { ref_pos as i64 };
                    for k in 0..len {
                        let rp = first_ref + i64::from(k);
                        let qp = (query_pos + k) as usize;
                        let Some(&ref_base) = u64::try_from(rp)
                            .ok()
                            .and_then(|rp| rp.checked_sub(offset))
                            .and_then(|o| seq.get(o as usize))
                        else {
                            continue;
                        };
                        let base = if qp < read.sequence_len() { base_letter(read.base(qp)) } else { b'N' };
                        if base != b'=' && ref_base != b'N' && base != ref_base {
                            out.clip_mismatches.push(Mismatch {
                                pos: rp as u64,
                                base,
                                quality: quals.get(qp).copied().unwrap_or(0xff),
                            });
                        }
                    }
                }
                query_pos += len;
            }
            _ => {} // H, P
        }
    }

    out.end = ref_pos.max(start + 1);
    if mods {
        out.mods = mod_calls(read, &out.blocks, start);
    }
    out
}

/// Walks MM/ML along the read (as HTSlib's `bam_mods_at_next_pos` does) and keeps the calls on
/// bases inside aligned blocks, at their reference positions.
fn mod_calls(read: &Read, blocks: &[Block], start: u64) -> Vec<ModCall> {
    if read.tag(*b"MM").is_none() && read.tag(*b"Mm").is_none() {
        return Vec::new();
    }
    let Ok(mut state) = BaseMods::parse(read) else { return Vec::new() };
    if state.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut calls = Vec::new();
    let mut block = 0;
    for qpos in 0..read.sequence_len() as u32 {
        if state.next_position(read, &mut calls).is_err() {
            break;
        }
        if calls.is_empty() {
            continue;
        }
        // Advance to the block holding this query position, if any (blocks are in query order).
        while block < blocks.len() && blocks[block].query_start + blocks[block].len <= qpos {
            block += 1;
        }
        let Some(b) = blocks.get(block).filter(|b| b.query_start <= qpos) else { continue };
        let pos = b.ref_start + u64::from(qpos - b.query_start);
        for c in &calls {
            let code = u8::try_from(c.modified_base).ok().filter(u8::is_ascii_alphabetic).unwrap_or(b'?');
            let prob = if c.qual < 0 { 255 } else { c.qual.clamp(0, 255) as u8 };
            out.push(ModCall { offset: (pos - start) as u32, code, canonical: c.canonical_base, prob });
        }
    }
    // Growth by doubling leaves up to half the vector spare, and these are the largest part of
    // a read with calls on every base.
    out.shrink_to_fit();
    out
}

fn base_letter(code: u8) -> u8 {
    match code {
        0 => b'=',
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Encodes a BAM record (without `block_size`) the way a BAM file stores it.
    fn record(pos: i32, flag: u16, cigar: &[(u32, u32)], seq: &[u8], qual: u8) -> Vec<u8> {
        let code = |b: u8| match b {
            b'=' => 0u8,
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        };
        let name = b"r1\0";
        let mut v = Vec::new();
        v.extend(0i32.to_le_bytes()); // refID
        v.extend(pos.to_le_bytes());
        v.push(name.len() as u8);
        v.push(60); // mapq
        v.extend(0u16.to_le_bytes()); // bin
        v.extend((cigar.len() as u16).to_le_bytes());
        v.extend(flag.to_le_bytes());
        v.extend((seq.len() as u32).to_le_bytes());
        v.extend((-1i32).to_le_bytes());
        v.extend((-1i32).to_le_bytes());
        v.extend(0i32.to_le_bytes());
        v.extend(name);
        for (len, op) in cigar {
            v.extend((len << 4 | op).to_le_bytes());
        }
        for pair in seq.chunks(2) {
            v.push(code(pair[0]) << 4 | pair.get(1).map_or(0, |b| code(*b)));
        }
        v.extend(std::iter::repeat_n(qual, seq.len()));
        v
    }

    #[test]
    fn parses_sa_tags() {
        let sa = parse_sa("chr2,1001,-,50S100M,60,0;chr1,5,+,100M50S,7,2;");
        assert_eq!(sa.len(), 2);
        assert_eq!((sa[0].contig.as_str(), sa[0].pos, sa[0].reverse, sa[0].mapq), ("chr2", 1000, true, 60));
        assert_eq!((sa[1].contig.as_str(), sa[1].pos, sa[1].reverse, sa[1].cigar.as_str()), ("chr1", 4, false, "100M50S"));
        assert!(parse_sa("").is_empty());
        assert!(parse_sa("chr1,0,+,10M,60,0;").is_empty(), "position 0 is invalid (1-based)");
    }

    #[test]
    fn haplotype_tags_of_any_type() {
        let base = record(100, 0, &[(4, 0)], b"ACGT", 30);
        let with = |tag: &[u8]| {
            let mut r = base.clone();
            r.extend_from_slice(tag);
            decode(&Read::new(r).unwrap(), None).haplotype
        };
        assert_eq!(with(b"HPC\x02"), Some(2), "HP:i (as C)");
        assert_eq!(with(b"HPZ1\0"), Some(1), "HP:Z:1");
        assert_eq!(with(b"HPZone\0"), None, "a non-numeric string is no haplotype, not 0");
        assert_eq!(with(b""), None);
    }

    #[test]
        fn decodes_cigar_and_mismatches() {
        // 2S 4M 1I 3M 2D 3M at 0-based 100.
        let cigar = [(2, S), (4, M), (1, I), (3, M), (2, D), (3, M)];
        let seq = b"NNACGAGAAAAAC";
        let read = Read::new(record(100, 0x10, &cigar, seq, 30)).unwrap();
        let reference = b"ACGTAAAAAAAAAAAA"; // starts at 100

        let r = decode(&read, Some((reference, 100)));

        assert_eq!(r.name, "r1");
        assert!(r.reverse);
        assert_eq!((r.start, r.end), (100, 112));
        assert_eq!(r.soft_clip_left, 2);
        assert_eq!(r.soft_clip_right, 0);
        assert_eq!(
            r.blocks,
            vec![
                Block { ref_start: 100, query_start: 2, len: 4 },
                Block { ref_start: 104, query_start: 7, len: 3 },
                Block { ref_start: 109, query_start: 10, len: 3 },
            ]
        );
        assert_eq!(r.insertions, vec![Indel { pos: 104, len: 1 }]);
        assert_eq!(r.deletions, vec![Indel { pos: 107, len: 2 }]);
        assert_eq!(
            r.mismatches,
            vec![
                Mismatch { pos: 103, base: b'A', quality: 30 },
                Mismatch { pos: 111, base: b'C', quality: 30 },
            ]
        );

        // Without a reference there is nothing to compare against.
        assert!(decode(&read, None).mismatches.is_empty());
    }
}

/// A tag's value as text: strings and characters as they are, numbers in decimal, arrays as
/// comma-separated numbers (as `samtools view` prints them).
pub fn tag_text(read: &Read, tag: [u8; 2]) -> Option<Box<str>> {
    let (ty, v) = read.tag(tag)?;
    let int = |ty: u8, b: &[u8]| -> Option<i64> {
        Some(match ty {
            b'c' => i64::from(*b.first()? as i8),
            b'C' => i64::from(*b.first()?),
            b's' => i64::from(i16::from_le_bytes(b.get(..2)?.try_into().ok()?)),
            b'S' => i64::from(u16::from_le_bytes(b.get(..2)?.try_into().ok()?)),
            b'i' => i64::from(i32::from_le_bytes(b.get(..4)?.try_into().ok()?)),
            b'I' => i64::from(u32::from_le_bytes(b.get(..4)?.try_into().ok()?)),
            _ => return None,
        })
    };
    let text = match ty {
        b'Z' | b'H' => String::from_utf8_lossy(v.strip_suffix(&[0]).unwrap_or(v)).into_owned(),
        b'A' => (*v.first()? as char).to_string(),
        b'f' => f32::from_le_bytes(v.get(..4)?.try_into().ok()?).to_string(),
        b'B' => {
            let sub = *v.first()?;
            let n = u32::from_le_bytes(v.get(1..5)?.try_into().ok()?) as usize;
            let size = match sub {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                _ => 4,
            };
            let items: Vec<String> = (0..n.min(64))
                .filter_map(|k| {
                    let b = v.get(5 + k * size..5 + (k + 1) * size)?;
                    if sub == b'f' { Some(f32::from_le_bytes(b.try_into().ok()?).to_string()) } else { int(sub, b).map(|x| x.to_string()) }
                })
                .collect();
            items.join(",")
        }
        ty => int(ty, v)?.to_string(),
    };
    Some(text.into())
}

/// A two-letter tag name (`CB`, `PS`, `HP`...).
pub fn parse_tag_name(s: &str) -> Option<[u8; 2]> {
    let b = s.trim().as_bytes();
    (b.len() == 2 && b[0].is_ascii_alphabetic() && b[1].is_ascii_alphanumeric()).then(|| [b[0], b[1]])
}
