use std::{fs::File, path::Path};

use ggv_model::{Contig, Region};
use pasta::{
    core::{Position, Region as PastaRegion},
    fasta,
};

use crate::{Error, Result, error::IoContext, resolve_locus};

/// A reference: an indexed FASTA (`.fa` + `.fa.fai`, or bgzipped `.fa.gz` + `.fai` + `.gzi`)
/// or a UCSC `.2bit`.
pub struct ReferenceFasta {
    path: std::path::PathBuf,
    reader: Inner,
    contigs: Vec<Contig>,
}

enum Inner {
    Fasta(Box<fasta::io::IndexedReader<fasta::io::BufReader<File>>>),
    TwoBit(crate::twobit::TwoBit),
}

impl ReferenceFasta {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_owned();
        if path.to_string_lossy().to_ascii_lowercase().ends_with(".2bit") {
            let t = crate::twobit::TwoBit::open(&path)?;
            let contigs = t.contigs().to_vec();
            return Ok(Self { path, reader: Inner::TwoBit(t), contigs });
        }
        // A missing index is named as such, not as the FASTA itself being missing (opening the
        // absent .fai/.gzi fails with "file not found", which was reported against the FASTA).
        if path.exists() {
            let gz = crate::indexing::is_gzip(&path)?;
            let needed = std::iter::once("fai").chain(gz.then_some("gzi"));
            let missing: Vec<String> = needed.map(|e| path.with_added_extension(e)).filter(|p| !p.exists()).map(|p| crate::assembly::display_name(&p)).collect();
            if !missing.is_empty() {
                let msg = format!("no FASTA index (looked for {}); create it with `ggv index` (or `samtools faidx`), or open the file in the GGV app, which offers to build it", missing.join(", "));
                return Err(Error::Io { path, cause: std::io::Error::new(std::io::ErrorKind::NotFound, msg) });
            }
        }
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(&path)
            .at(&path)?;
        let contigs = reader
            .index()
            .as_ref()
            .iter()
            .map(|r| Contig { name: r.name().to_string(), length: r.length() })
            .collect();
        Ok(Self { path, reader: Inner::Fasta(Box::new(reader)), contigs })
    }

    pub fn contigs(&self) -> &[Contig] {
        &self.contigs
    }

    /// Fetches the uppercase reference bases for `region`. The region's contig is matched to the
    /// FASTA's own naming (`chr7` ⇄ `7`), so a BAM and reference with different conventions
    /// still line up.
    pub fn fetch(&mut self, region: &Region) -> Result<Vec<u8>> {
        let locus = ggv_model::Locus {
            contig: region.contig.clone(),
            start: Some(region.start + 1),
            end: Some(region.end),
        };
        let region = resolve_locus(&self.contigs, &locus).map_err(|e| e.in_reference(&self.path))?;
        let reader = match &mut self.reader {
            Inner::TwoBit(t) => return t.fetch(&region.contig, region.start, region.end),
            Inner::Fasta(r) => r,
        };

        let start = Position::try_from(region.start as usize + 1).expect("start is >= 1");
        let end = Position::try_from(region.end as usize).expect("end is >= 1");
        let record = reader
            .query(&PastaRegion::new(region.contig.as_str(), start..=end))
            .at(&self.path)?;

        let mut seq = record.sequence().as_ref().to_vec();
        seq.make_ascii_uppercase();
        Ok(seq)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_missing_index_is_named_not_the_fasta() {
        let dir = std::env::temp_dir().join(format!("ggv-ref-noindex-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let plain = dir.join("small.fa");
        std::fs::write(&plain, ">c\nACGT\n").unwrap();
        let e = ReferenceFasta::open(&plain).err().expect("no .fai").to_string();
        assert!(e.contains("no FASTA index (looked for small.fa.fai)") && e.contains("samtools faidx") && e.contains("GGV app"), "{e}");
        crate::build_index(&plain).unwrap();
        assert_eq!(ReferenceFasta::open(&plain).unwrap().contigs().len(), 1);

        // Bgzipped: the .gzi is needed too.
        let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/hg002-revio/GRCh38.chr20.slices.fa.gz");
        if src.exists() {
            let gz = dir.join("ref.fa.gz");
            std::fs::copy(&src, &gz).unwrap();
            let e = ReferenceFasta::open(&gz).err().expect("no .fai or .gzi").to_string();
            assert!(e.contains("looked for ref.fa.gz.fai, ref.fa.gz.gzi"), "{e}");
            std::fs::copy(src.with_added_extension("fai"), gz.with_added_extension("fai")).unwrap();
            let e = ReferenceFasta::open(&gz).err().expect("no .gzi").to_string();
            assert!(e.contains("looked for ref.fa.gz.gzi)"), "{e}");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}
