//! Bookmarks: named places with notes, kept in the session and exportable as BED. Shared by
//! every pane of a window.

use std::{cell::RefCell, rc::Rc};

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Bookmark {
    pub contig: String,
    /// 0-based, half-open.
    pub start: u64,
    pub end: u64,
    #[serde(default)]
    pub note: String,
}

impl Bookmark {
    pub fn locus(&self) -> String {
        if self.end == self.start + 1 {
            format!("{}:{}", self.contig, crate::util::commas(self.start + 1))
        } else {
            format!("{}:{}-{}", self.contig, crate::util::commas(self.start + 1), crate::util::commas(self.end))
        }
    }
}

/// The window's bookmarks, shared by its panes.
pub type Shared = Rc<RefCell<Vec<Bookmark>>>;

/// BED: contig, start, end, and the note (spaces kept, tabs and newlines replaced) as the name.
pub fn to_bed(bookmarks: &[Bookmark]) -> String {
    let mut out = String::new();
    for b in bookmarks {
        let name = b.note.replace(['\t', '\n', '\r'], " ");
        out.push_str(&format!("{}\t{}\t{}\t{}\n", b.contig, b.start, b.end, if name.is_empty() { "." } else { &name }));
    }
    out
}

/// Adds `b`, or updates the note of one already at the same place. Returns its index.
pub fn add(list: &mut Vec<Bookmark>, b: Bookmark) -> usize {
    match list.iter().position(|x| x.contig == b.contig && x.start == b.start && x.end == b.end) {
        Some(i) => {
            if !b.note.is_empty() {
                list[i].note = b.note;
            }
            i
        }
        None => {
            list.push(b);
            list.len() - 1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed_export_and_dedup() {
        let mut l = Vec::new();
        add(&mut l, Bookmark { contig: "chr1".into(), start: 99, end: 100, note: String::new() });
        add(&mut l, Bookmark { contig: "chr2".into(), start: 0, end: 500, note: "odd\tdip".into() });
        let i = add(&mut l, Bookmark { contig: "chr1".into(), start: 99, end: 100, note: "het".into() });
        assert_eq!((l.len(), i), (2, 0));
        assert_eq!(to_bed(&l), "chr1\t99\t100\thet\nchr2\t0\t500\todd dip\n");
        assert_eq!(l[0].locus(), "chr1:100");
        assert_eq!(l[1].locus(), "chr2:1-500");
    }
}
