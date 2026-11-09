//! Overlays: regions an agent (or a session) marks on the view, each a tinted band across every
//! track with a label, such as "deletion at the run start, insertion at the site". Shared by a
//! window's panes and saved with the session; the person clears them from the toolbar.

use std::{cell::RefCell, rc::Rc};

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Overlay {
    pub contig: String,
    /// 0-based, half-open.
    pub start: u64,
    pub end: u64,
    #[serde(default)]
    pub label: String,
    /// One of [`COLORS`].
    #[serde(default = "default_color")]
    pub color: String,
}

fn default_color() -> String {
    "yellow".into()
}

/// The colours an overlay can take, by name.
pub const COLORS: [(&str, [u8; 3]); 7] = [
    ("yellow", [235, 190, 30]),
    ("red", [215, 60, 60]),
    ("green", [40, 160, 80]),
    ("blue", [60, 110, 220]),
    ("purple", [150, 80, 200]),
    ("orange", [235, 130, 30]),
    ("grey", [130, 130, 130]),
];

pub fn rgb(name: &str) -> Option<[u8; 3]> {
    COLORS.iter().find(|(n, _)| *n == name).map(|(_, c)| *c)
}

/// The window's overlays, shared by its panes.
pub type Shared = Rc<RefCell<Vec<Overlay>>>;

impl Overlay {
    pub fn locus(&self) -> String {
        if self.end == self.start + 1 {
            format!("{}:{}", self.contig, crate::util::commas(self.start + 1))
        } else {
            format!("{}:{}-{}", self.contig, crate::util::commas(self.start + 1), crate::util::commas(self.end))
        }
    }
}
