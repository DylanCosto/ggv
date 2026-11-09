//! MANE (Matched Annotation from NCBI and EMBL-EBI) GRCh38 v1.5: the RefSeq and Ensembl
//! transcript of each human protein-coding gene that clinical reporting uses (MANE Select), and
//! the extra ones it adds for some genes (MANE Plus Clinical). Shipped because UCSC's ncbiRefSeq
//! GTF has no MANE tags, and no simple rule over a gene's transcripts finds them (TTN, DMD,
//! NF1). Public domain (NCBI/EMBL-EBI); see THIRD_PARTY_NOTICES.md.
//!
//! The list is gzipped TSV: a `#` comment line, then `RefSeq id`, `Ensembl id`, `S` or `P`.

use std::{collections::HashMap, io::Read, sync::OnceLock};

static DATA: &[u8] = include_bytes!("../data/mane-grch38-v1.5.tsv.gz");

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Kind {
    PlusClinical,
    Select,
}

/// How a transcript id matches the list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Hit {
    pub kind: Kind,
    /// Same version as listed (else the same accession, another version).
    pub exact: bool,
}

/// Unversioned accession → (listed version, kind).
type Table = HashMap<Box<str>, (Box<str>, Kind)>;

fn table() -> &'static Table {
    static TABLE: OnceLock<Table> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut text = String::new();
        // Embedded and tested: a failure here is a broken build, which leaves the table empty.
        let _ = flate2::read::GzDecoder::new(DATA).read_to_string(&mut text);
        let mut out = HashMap::new();
        for line in text.lines().filter(|l| !l.starts_with('#')) {
            let mut f = line.split('\t');
            let (Some(refseq), Some(ensembl), Some(kind)) = (f.next(), f.next(), f.next()) else { continue };
            let kind = if kind == "P" { Kind::PlusClinical } else { Kind::Select };
            for id in [refseq, ensembl] {
                let (base, version) = id.split_once('.').unwrap_or((id, ""));
                out.insert(base.into(), (version.into(), kind));
            }
        }
        out
    })
}

/// A transcript id without UCSC's copy suffix: `NM_000314.8_2` → `NM_000314.8` (UCSC numbers
/// the copies of a transcript placed on several sequences).
pub fn without_copy_suffix(id: &str) -> &str {
    match id.rsplit_once('_') {
        Some((base, n)) if !n.is_empty() && n.bytes().all(|b| b.is_ascii_digit()) && base.rsplit_once('.').is_some_and(|(_, v)| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit())) => base,
        _ => id,
    }
}

/// Where `id` (RefSeq or Ensembl, versioned or not, with or without a copy suffix) is in MANE.
pub(crate) fn lookup(id: &str) -> Option<Hit> {
    // Only RefSeq and Ensembl transcript ids can be listed (NR_ for the non-coding genes, such
    // as H19 and GAS5): don't decode the list for others.
    if !(id.starts_with("NM_") || id.starts_with("NR_") || id.starts_with("ENST")) {
        return None;
    }
    let id = without_copy_suffix(id);
    let (base, version) = id.split_once('.').unwrap_or((id, ""));
    let (listed, kind) = table().get(base)?;
    Some(Hit { kind: *kind, exact: version.is_empty() || **listed == *version })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_list_decodes_and_matches_ids() {
        assert!(table().len() > 38_000, "{}", table().len());
        let select = |exact| Some(Hit { kind: Kind::Select, exact });
        assert_eq!(lookup("NM_004333.6"), select(true));
        assert_eq!(lookup("NM_004333.5"), select(false));
        assert_eq!(lookup("NM_004333"), select(true));
        assert_eq!(lookup("NM_000314.8_2"), select(true));
        assert_eq!(lookup("ENST00000646891.2"), select(true), "BRAF's Ensembl id");
        assert_eq!(lookup("NM_001374258.1").map(|h| h.kind), Some(Kind::PlusClinical), "BRAF's Plus Clinical");
        assert_eq!(lookup("NR_002196.3"), select(true), "H19, a non-coding MANE Select");
        assert_eq!(lookup("NR_002196.2_2"), select(false), "UCSC's older H19 version");
        assert_eq!(lookup("XM_017012559.2"), None);
        assert_eq!(lookup("NM_001374244.1"), None, "a BRAF isoform not in MANE");
        assert_eq!(without_copy_suffix("NM_000314.8_2"), "NM_000314.8");
        assert_eq!(without_copy_suffix("NM_000314"), "NM_000314");
        assert_eq!(without_copy_suffix("my_tx_2"), "my_tx_2");
    }
}
