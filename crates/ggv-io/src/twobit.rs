//! UCSC `.2bit` references: random access without unpacking (UCSC's FASTA downloads are plain
//! gzip, which can't be read at a position). Format: <https://genome.ucsc.edu/FAQ/FAQformat.html#format7>.
//!
//! A header and an index of sequence names to offsets; per sequence, its length, blocks of `N`,
//! soft-mask blocks (ignored: bases come back uppercase, as everywhere in GGV) and the bases at
//! two bits each (T, C, A, G = 0-3, first base in the high bits). Version 1 files (genomes over
//! 4 GB) use 64-bit offsets.

use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use ggv_model::Contig;

/// One sequence's layout, read on first use.
struct Layout {
    len: u64,
    /// Sorted, non-overlapping (start, end) runs of N.
    n_blocks: Vec<(u64, u64)>,
    /// File offset of the packed bases.
    dna: u64,
}

pub struct TwoBit {
    path: PathBuf,
    file: File,
    little: bool,
    contigs: Vec<Contig>,
    offsets: HashMap<String, u64>,
    layouts: HashMap<String, Layout>,
}

fn err(path: &Path, msg: impl Into<String>) -> crate::Error {
    crate::Error::Io { path: path.to_owned(), cause: std::io::Error::other(msg.into()) }
}

impl TwoBit {
    pub fn open(path: &Path) -> crate::Result<Self> {
        let io = |e: std::io::Error| crate::Error::Io { path: path.to_owned(), cause: e };
        let mut file = File::open(path).map_err(io)?;
        let mut head = [0u8; 16];
        file.read_exact(&mut head).map_err(io)?;
        let little = match head[..4] {
            [0x43, 0x27, 0x41, 0x1A] => true,
            [0x1A, 0x41, 0x27, 0x43] => false,
            _ => return Err(err(path, "not a .2bit file")),
        };
        let u32_at = |b: &[u8]| if little { u32::from_le_bytes(b.try_into().unwrap()) } else { u32::from_be_bytes(b.try_into().unwrap()) };
        let version = u32_at(&head[4..8]);
        if version > 1 {
            return Err(err(path, format!(".2bit version {version} is not supported")));
        }
        let count = u32_at(&head[8..12]) as usize;
        let mut me = Self { path: path.to_owned(), file, little, contigs: Vec::with_capacity(count), offsets: HashMap::with_capacity(count), layouts: HashMap::new() };
        let mut names = Vec::with_capacity(count);
        for _ in 0..count {
            let n = me.read_u8()? as usize;
            let mut name = vec![0u8; n];
            me.file.read_exact(&mut name).map_err(io)?;
            let offset = if version == 1 { me.read_u64()? } else { u64::from(me.read_u32()?) };
            let name = String::from_utf8_lossy(&name).into_owned();
            me.offsets.insert(name.clone(), offset);
            names.push((name, offset));
        }
        // Lengths are at each sequence's offset (one small read each).
        for (name, offset) in names {
            me.file.seek(SeekFrom::Start(offset)).map_err(io)?;
            let len = u64::from(me.read_u32()?);
            me.contigs.push(Contig { name, length: len });
        }
        Ok(me)
    }

    pub fn contigs(&self) -> &[Contig] {
        &self.contigs
    }

    fn io(&self) -> impl Fn(std::io::Error) -> crate::Error + '_ {
        move |e| crate::Error::Io { path: self.path.clone(), cause: e }
    }

    fn read_u8(&mut self) -> crate::Result<u8> {
        let mut b = [0u8; 1];
        self.file.read_exact(&mut b).map_err(|e| crate::Error::Io { path: self.path.clone(), cause: e })?;
        Ok(b[0])
    }

    fn read_u32(&mut self) -> crate::Result<u32> {
        let mut b = [0u8; 4];
        self.file.read_exact(&mut b).map_err(|e| crate::Error::Io { path: self.path.clone(), cause: e })?;
        Ok(if self.little { u32::from_le_bytes(b) } else { u32::from_be_bytes(b) })
    }

    fn read_u64(&mut self) -> crate::Result<u64> {
        let mut b = [0u8; 8];
        self.file.read_exact(&mut b).map_err(|e| crate::Error::Io { path: self.path.clone(), cause: e })?;
        Ok(if self.little { u64::from_le_bytes(b) } else { u64::from_be_bytes(b) })
    }

    fn read_u32s(&mut self, n: usize) -> crate::Result<Vec<u32>> {
        let mut raw = vec![0u8; n * 4];
        self.file.read_exact(&mut raw).map_err(|e| crate::Error::Io { path: self.path.clone(), cause: e })?;
        Ok(raw.chunks_exact(4).map(|c| if self.little { u32::from_le_bytes(c.try_into().unwrap()) } else { u32::from_be_bytes(c.try_into().unwrap()) }).collect())
    }

    fn layout(&mut self, name: &str) -> crate::Result<&Layout> {
        if !self.layouts.contains_key(name) {
            let offset = *self.offsets.get(name).ok_or_else(|| err(&self.path, format!("{name}: not in this .2bit")))?;
            self.file.seek(SeekFrom::Start(offset)).map_err(self.io())?;
            let len = u64::from(self.read_u32()?);
            let n = self.read_u32()? as usize;
            let starts = self.read_u32s(n)?;
            let sizes = self.read_u32s(n)?;
            let masks = self.read_u32()? as usize;
            // Mask starts and sizes, then a reserved word: skipped.
            let dna = offset + 4 + 4 + 8 * n as u64 + 4 + 8 * masks as u64 + 4;
            let mut n_blocks: Vec<(u64, u64)> = starts.iter().zip(&sizes).map(|(&s, &z)| (u64::from(s), u64::from(s) + u64::from(z))).collect();
            n_blocks.sort_unstable();
            self.layouts.insert(name.to_owned(), Layout { len, n_blocks, dna });
        }
        Ok(&self.layouts[name])
    }

    /// Uppercase bases of `name` over the 0-based half-open `start..end` (clamped to the
    /// sequence).
    pub fn fetch(&mut self, name: &str, start: u64, end: u64) -> crate::Result<Vec<u8>> {
        let (len, dna, n_blocks) = {
            let l = self.layout(name)?;
            (l.len, l.dna, l.n_blocks.clone())
        };
        let (start, end) = (start.min(len), end.min(len));
        if start >= end {
            return Ok(Vec::new());
        }
        let (first_byte, last_byte) = (start / 4, (end - 1) / 4);
        let mut packed = vec![0u8; (last_byte - first_byte + 1) as usize];
        self.file.seek(SeekFrom::Start(dna + first_byte)).map_err(self.io())?;
        self.file.read_exact(&mut packed).map_err(self.io())?;
        const BASES: [u8; 4] = *b"TCAG";
        let mut out: Vec<u8> = (start..end)
            .map(|p| {
                let byte = packed[(p / 4 - first_byte) as usize];
                BASES[usize::from((byte >> (6 - 2 * (p % 4))) & 3)]
            })
            .collect();
        // N blocks overlapping the span.
        let i = n_blocks.partition_point(|&(_, e)| e <= start);
        for &(s, e) in n_blocks[i..].iter().take_while(|&&(s, _)| s < end) {
            for b in &mut out[(s.max(start) - start) as usize..(e.min(end) - start) as usize] {
                *b = b'N';
            }
        }
        Ok(out)
    }
}

/// Writes a `.2bit` (version 0, little-endian) of `seqs`, soft-masking lowercase runs; for tests
/// and tools.
pub fn write(path: &Path, seqs: &[(&str, &[u8])]) -> std::io::Result<()> {
    let runs = |s: &[u8], pred: &dyn Fn(u8) -> bool| -> Vec<(u32, u32)> {
        let mut v = Vec::new();
        let mut i = 0;
        while i < s.len() {
            if pred(s[i]) {
                let j = (i..s.len()).find(|&j| !pred(s[j])).unwrap_or(s.len());
                v.push((i as u32, (j - i) as u32));
                i = j;
            } else {
                i += 1;
            }
        }
        v
    };
    let mut out = Vec::new();
    out.extend(0x1A41_2743u32.to_le_bytes());
    out.extend(0u32.to_le_bytes());
    out.extend((seqs.len() as u32).to_le_bytes());
    out.extend(0u32.to_le_bytes());
    let index_len: usize = seqs.iter().map(|(n, _)| 1 + n.len() + 4).sum();
    let mut offset = (16 + index_len) as u32;
    let mut bodies = Vec::new();
    for (name, seq) in seqs {
        out.push(name.len() as u8);
        out.extend(name.as_bytes());
        out.extend(offset.to_le_bytes());
        let ns = runs(seq, &|b| b.eq_ignore_ascii_case(&b'N'));
        let masks = runs(seq, &|b| b.is_ascii_lowercase());
        let mut body = Vec::new();
        body.extend((seq.len() as u32).to_le_bytes());
        for blocks in [&ns, &masks] {
            body.extend((blocks.len() as u32).to_le_bytes());
            blocks.iter().for_each(|b| body.extend(b.0.to_le_bytes()));
            blocks.iter().for_each(|b| body.extend(b.1.to_le_bytes()));
        }
        body.extend(0u32.to_le_bytes());
        for chunk in seq.chunks(4) {
            let mut byte = 0u8;
            for (k, &b) in chunk.iter().enumerate() {
                let code = match b.to_ascii_uppercase() {
                    b'C' => 1,
                    b'A' => 2,
                    b'G' => 3,
                    _ => 0,
                };
                byte |= code << (6 - 2 * k);
            }
            body.push(byte);
        }
        offset += body.len() as u32;
        bodies.push(body);
    }
    for b in bodies {
        out.extend(b);
    }
    std::fs::write(path, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_what_was_written() {
        let dir = std::env::temp_dir().join(format!("ggv-2bit-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.2bit");
        let a: Vec<u8> = b"NNNNacgtACGTTTGGCCAANNNNNNNNggaatt".to_vec();
        let b: Vec<u8> = (0..10_007).map(|i| b"ACGT"[(i * 7 + i / 13) % 4]).collect();
        write(&path, &[("chrA", &a), ("chrB", &b)]).unwrap();
        let mut t = TwoBit::open(&path).unwrap();
        assert_eq!(t.contigs().iter().map(|c| (c.name.as_str(), c.length)).collect::<Vec<_>>(), [("chrA", 34), ("chrB", 10_007)]);
        let up = |s: &[u8]| s.to_ascii_uppercase();
        assert_eq!(t.fetch("chrA", 0, 34).unwrap(), up(&a));
        // Every start and length near byte boundaries, and spans crossing N blocks.
        for s in 0..34u64 {
            for e in s..=34 {
                assert_eq!(t.fetch("chrA", s, e).unwrap(), up(&a[s as usize..e as usize]), "{s}..{e}");
            }
        }
        for (s, e) in [(0, 10_007), (4_999, 5_013), (10_000, 10_007), (3, 4)] {
            assert_eq!(t.fetch("chrB", s, e).unwrap(), b[s as usize..e as usize], "{s}..{e}");
        }
        assert_eq!(t.fetch("chrB", 10_005, 20_000).unwrap(), b[10_005..]);
        assert!(t.fetch("chrZ", 0, 1).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
