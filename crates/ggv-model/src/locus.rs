use std::{fmt, str::FromStr};

/// A locus as a person or agent writes it, before it is checked against any file.
///
/// Accepted forms (positions are 1-based, inclusive, and may contain `,` or `_` separators):
///
/// - `chr7` — the whole contig
/// - `chr7:55191822` — one base
/// - `chr7:55,191,822-55,191,900` — a range
///
/// Contig names may themselves contain `:` (e.g. `HLA-A*01:01:01:01`); if the text after the
/// last `:` is not a position or range, the whole string is taken as the contig name. Text alone
/// can't settle every case (`HLA-A*01:01:01:01` ends in what looks like position 1), so
/// `ggv_io::resolve_locus_str` first checks the whole string against the file's contigs.
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Locus {
    pub contig: String,
    /// 1-based inclusive start; `None` means the start of the contig.
    pub start: Option<u64>,
    /// 1-based inclusive end; `None` means the end of the contig.
    pub end: Option<u64>,
}

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum LocusParseError {
    #[error("empty locus")]
    Empty,
    #[error("position must be 1 or greater")]
    ZeroPosition,
    #[error("start {start} is after end {end}")]
    StartAfterEnd { start: u64, end: u64 },
}

impl Locus {
    pub fn whole(contig: impl Into<String>) -> Self {
        Self { contig: contig.into(), start: None, end: None }
    }
}

fn parse_pos(s: &str) -> Option<u64> {
    let digits: String = s.chars().filter(|c| *c != ',' && *c != '_').collect();
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    digits.parse().ok()
}

impl FromStr for Locus {
    type Err = LocusParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err(LocusParseError::Empty);
        }

        let Some((contig, suffix)) = s.rsplit_once(':') else {
            return Ok(Self::whole(s));
        };

        let (start, end) = match suffix.split_once('-') {
            Some((a, b)) => match (parse_pos(a), parse_pos(b)) {
                (Some(a), Some(b)) => (a, b),
                _ => return Ok(Self::whole(s)),
            },
            None => match parse_pos(suffix) {
                Some(p) => (p, p),
                None => return Ok(Self::whole(s)),
            },
        };

        if contig.is_empty() {
            return Err(LocusParseError::Empty);
        }
        if start == 0 || end == 0 {
            return Err(LocusParseError::ZeroPosition);
        }
        if start > end {
            return Err(LocusParseError::StartAfterEnd { start, end });
        }

        Ok(Self { contig: contig.to_owned(), start: Some(start), end: Some(end) })
    }
}

impl fmt::Display for Locus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.start, self.end) {
            (Some(s), Some(e)) if s == e => write!(f, "{}:{s}", self.contig),
            (Some(s), Some(e)) => write!(f, "{}:{s}-{e}", self.contig),
            (Some(s), None) => write!(f, "{}:{s}-", self.contig),
            _ => write!(f, "{}", self.contig),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(s: &str) -> Locus {
        s.parse().unwrap()
    }

    #[test]
    fn whole_contig() {
        assert_eq!(p("chr7"), Locus::whole("chr7"));
    }

    #[test]
    fn single_base_with_commas() {
        assert_eq!(p("chr7:55,191,822"), Locus { contig: "chr7".into(), start: Some(55191822), end: Some(55191822) });
    }

    #[test]
    fn range() {
        assert_eq!(p(" 7:100-200 "), Locus { contig: "7".into(), start: Some(100), end: Some(200) });
    }

    #[test]
    fn colon_in_contig_name() {
        assert_eq!(p("HLA-DRB1*15:01:01:01N"), Locus::whole("HLA-DRB1*15:01:01:01N"));
        assert_eq!(
            p("HLA-A*01:01:01:01:5-10"),
            Locus { contig: "HLA-A*01:01:01:01".into(), start: Some(5), end: Some(10) }
        );
    }

    #[test]
    fn errors() {
        assert_eq!("".parse::<Locus>(), Err(LocusParseError::Empty));
        assert_eq!("chr1:0".parse::<Locus>(), Err(LocusParseError::ZeroPosition));
        assert_eq!("chr1:20-10".parse::<Locus>(), Err(LocusParseError::StartAfterEnd { start: 20, end: 10 }));
    }

    #[test]
    fn display_round_trips() {
        for s in ["chr7", "chr7:5", "chr7:5-9"] {
            assert_eq!(p(s).to_string(), s);
        }
    }
}
