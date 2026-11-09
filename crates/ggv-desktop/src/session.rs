//! Sessions: what's open and how it's shown, as a small JSON file. Display options are stored
//! by name (not enum layout) so files stay readable and survive refactors.
//!
//! Files inside the session's folder are stored relative to it, so a project folder (session
//! and data together) can be moved or shared. Saves are atomic with a backup (see `store`).

use std::path::{Path, PathBuf};

use ggv_layout::{GroupBy, SortBy};
use serde::{Deserialize, Serialize};

use crate::tracks::{ColorBy, CoverageScale, ReadFilter, Style};

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Session {
    pub version: u32,
    pub bams: Vec<PathBuf>,
    pub vcfs: Vec<PathBuf>,
    pub genes: Vec<PathBuf>,
    /// Signal tracks (bigWig, bedGraph, SEG).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub signals: Vec<PathBuf>,
    pub reference: Option<PathBuf>,
    /// The view, e.g. `chr17:43,090,000-43,100,000`.
    pub locus: Option<String>,
    /// The outlined base, 1-based, e.g. `chr17:43092919`.
    pub highlight: Option<String>,
    pub color: String,
    pub group: String,
    /// For `group: "base"`: the 1-based position.
    pub group_pos: Option<u64>,
    /// Sort: `{"pos": 43092919, "by": "base"}` (1-based position).
    pub sort: Option<SortState>,
    pub soft_clips: bool,
    /// Split reads linked on one row.
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub link_splits: bool,
    pub squished: bool,
    /// Splice-junction arcs over the coverage.
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub junctions: bool,
    pub all_transcripts: bool,
    /// `view`, `loaded` or `log`.
    pub coverage_scale: String,
    pub shared_scale: bool,
    pub filter: Option<FilterState>,
    pub dark: bool,
    /// Per-file filter expressions that replace the shared one for that track.
    #[serde(skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub track_filters: std::collections::BTreeMap<PathBuf, String>,
    /// Variant tracks coloured other than their default: `af`, or an INFO key (`CLNSIG`).
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub variant_colors: std::collections::BTreeMap<PathBuf, String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub bookmarks: Vec<crate::bookmarks::Bookmark>,
    /// Regions marked on the view, with labels.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub overlays: Vec<crate::overlays::Overlay>,
    /// Track order and heights, when the person arranged them (top to bottom).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub track_layout: Vec<crate::app::arrange::Saved>,
    /// Every pane, left to right, when the window was split (the first is `locus`).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub panes: Vec<PaneView>,
    /// Split panes scroll and zoom together.
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub linked: bool,
    /// The pane last clicked, 1-based.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub active_pane: Option<usize>,
    /// Loaded from the backup because the file itself was damaged.
    #[serde(skip)]
    pub recovered: bool,
}

/// One pane of a split window: where it was and its outlined base.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct PaneView {
    pub locus: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub highlight: Option<String>,
}

/// The window's panes, as the shell last saw them (the shell runs one window per process).
static WORKSPACE: std::sync::Mutex<Option<(Vec<PaneView>, bool, usize)>> = std::sync::Mutex::new(None);

/// Records the window's panes, whether linked, and the active pane (0-based).
pub fn set_workspace(panes: Vec<PaneView>, linked: bool, active: usize) {
    if let Ok(mut w) = WORKSPACE.lock() {
        *w = Some((panes, linked, active));
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FilterState {
    pub hide_duplicates: bool,
    pub hide_secondary: bool,
    pub hide_supplementary: bool,
    pub hide_qc_fail: bool,
    pub min_mapq: u8,
    /// A filter expression, e.g. `mapq >= 20 && hp == 1`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expr: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SortState {
    pub pos: u64,
    pub by: String,
}

pub const VERSION: u32 = 1;

pub fn color_name(c: ColorBy) -> &'static str {
    match c {
        ColorBy::None => "none",
        ColorBy::Strand => "strand",
        ColorBy::Haplotype => "haplotype",
        ColorBy::InsertSize => "insert size",
        ColorBy::BaseMods => "base modifications",
        ColorBy::Tag(_) => "tag",
    }
}

/// The colouring as stored and shown: `tag:CB` for a tag.
pub fn color_text(c: ColorBy) -> String {
    match c {
        ColorBy::Tag(t) => format!("tag:{}", String::from_utf8_lossy(&t)),
        c => color_name(c).to_owned(),
    }
}

/// The grouping as stored: `tag:CB` for a tag.
pub fn group_text(g: GroupBy) -> String {
    match g {
        GroupBy::Tag(t) => format!("tag:{}", String::from_utf8_lossy(&t)),
        g => group_name(g).to_owned(),
    }
}

pub fn group_name(g: GroupBy) -> &'static str {
    match g {
        GroupBy::None => "none",
        GroupBy::Strand => "strand",
        GroupBy::Haplotype => "haplotype",
        GroupBy::ReadGroup => "read group",
        GroupBy::BaseAt(_) => "base",
        GroupBy::Alleles => "alleles",
        GroupBy::Tag(_) => "tag",
    }
}

pub fn sort_name(s: SortBy) -> &'static str {
    match s {
        SortBy::Base => "base",
        SortBy::Strand => "strand",
        SortBy::MappingQuality => "MAPQ",
        SortBy::Start => "start",
        SortBy::InsertSize => "insert size",
    }
}

pub(crate) fn parse_color(s: &str) -> ColorBy {
    if let Some(t) = s.strip_prefix("tag:").and_then(ggv_layout::parse_tag_name) {
        return ColorBy::Tag(t);
    }
    [ColorBy::Strand, ColorBy::Haplotype, ColorBy::InsertSize, ColorBy::BaseMods]
        .into_iter()
        .find(|c| color_name(*c) == s)
        .unwrap_or_default()
}

pub(crate) fn parse_group(s: &str, pos: Option<u64>) -> GroupBy {
    if let Some(t) = s.strip_prefix("tag:").and_then(ggv_layout::parse_tag_name) {
        return GroupBy::Tag(t);
    }
    if s == "base" {
        return pos.and_then(|p| p.checked_sub(1)).map_or(GroupBy::None, GroupBy::BaseAt);
    }
    [GroupBy::Strand, GroupBy::Haplotype, GroupBy::ReadGroup, GroupBy::Alleles].into_iter().find(|g| group_name(*g) == s).unwrap_or_default()
}

pub(crate) fn parse_sort(s: &str) -> Option<SortBy> {
    [SortBy::Base, SortBy::Strand, SortBy::MappingQuality, SortBy::Start, SortBy::InsertSize]
        .into_iter()
        .find(|b| sort_name(*b) == s)
}

impl Session {
    pub fn style(&self) -> Style {
        Style {
            color: parse_color(&self.color),
            group: parse_group(&self.group, self.group_pos),
            sort: self.sort.as_ref().and_then(|s| Some((s.pos.checked_sub(1)?, parse_sort(&s.by)?))),
            soft_clips: self.soft_clips,
            link_splits: self.link_splits,
            squished: self.squished,
            junctions: self.junctions,
            detail: Default::default(),
            coverage_scale: match self.coverage_scale.as_str() {
                "loaded" => CoverageScale::Loaded,
                "log" => CoverageScale::Log,
                _ => CoverageScale::View,
            },
            shared_scale: self.shared_scale,
            dark: self.dark,
            filter: self.filter.as_ref().map_or_else(ReadFilter::default, |f| ReadFilter {
                hide_duplicates: f.hide_duplicates,
                hide_secondary: f.hide_secondary,
                hide_supplementary: f.hide_supplementary,
                hide_qc_fail: f.hide_qc_fail,
                min_mapq: f.min_mapq,
                expr: f.expr.as_deref().and_then(|e| crate::expr::intern(e).ok()),
            }),
        }
    }

    pub fn set_style(&mut self, style: Style) {
        self.color = color_text(style.color);
        self.group = group_text(style.group);
        self.group_pos = match style.group {
            GroupBy::BaseAt(p) => Some(p + 1),
            _ => None,
        };
        self.sort = style.sort.map(|(pos, by)| SortState { pos: pos + 1, by: sort_name(by).to_owned() });
        self.soft_clips = style.soft_clips;
        self.link_splits = style.link_splits;
        self.squished = style.squished;
        self.junctions = style.junctions;
        self.coverage_scale = match style.coverage_scale {
            CoverageScale::View => "view",
            CoverageScale::Loaded => "loaded",
            CoverageScale::Log => "log",
        }
        .to_owned();
        self.shared_scale = style.shared_scale;
        self.dark = style.dark;
        let f = style.filter;
        self.filter = Some(FilterState {
            hide_duplicates: f.hide_duplicates,
            hide_secondary: f.hide_secondary,
            hide_supplementary: f.hide_supplementary,
            hide_qc_fail: f.hide_qc_fail,
            min_mapq: f.min_mapq,
            expr: f.expr.map(crate::expr::text),
        });
    }

    /// Nothing open: no file of any kind, nor a reference.
    pub fn is_empty(&self) -> bool {
        self.bams.is_empty() && self.vcfs.is_empty() && self.genes.is_empty() && self.signals.is_empty() && self.reference.is_none()
    }

    /// Reads a session (or, when the file is damaged, its backup), with relative paths
    /// resolved against the session's folder.
    pub fn load(path: &Path) -> anyhow::Result<Self> {
        // Only unreadable JSON is damage (and opens the backup); a session that reads but is
        // refused (newer version, a filter that doesn't parse) says why instead.
        let parse = |text: &str| -> anyhow::Result<Self> { Ok(serde_json::from_str(text)?) };
        let (mut s, recovered) = crate::store::read_recovering(path, parse).map_err(|e| anyhow::anyhow!("{}: {e:#}", path.display()))?;
        anyhow::ensure!(s.version <= VERSION, "{}: session version {} is newer than this GGV", path.display(), s.version);
        s.validate().map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
        s.recovered = recovered;
        if let Some(dir) = path.parent() {
            // Normalised: saved relative paths use `/`, which joined to a Windows folder gave
            // mixed separators that agents' spellings of the same file didn't match.
            s.map_paths(|p| crate::files::normalize(&if p.is_relative() && !p.as_os_str().is_empty() { dir.join(p) } else { p.to_path_buf() }));
        }
        Ok(s)
    }

    /// Fills in the window's other panes (when split) from what the shell last recorded.
    pub fn with_workspace(mut self) -> Self {
        if let Ok(w) = WORKSPACE.lock()
            && let Some((panes, linked, active)) = w.as_ref()
            && panes.len() > 1
        {
            self.panes = panes.clone();
            self.linked = *linked;
            self.active_pane = Some(active + 1);
        }
        self
    }

    /// Every file path in the session, changed by `f`.
    fn map_paths(&mut self, f: impl Fn(&Path) -> PathBuf) {
        for list in [&mut self.bams, &mut self.vcfs, &mut self.genes, &mut self.signals] {
            for p in list.iter_mut() {
                *p = f(p);
            }
        }
        if let Some(r) = &mut self.reference {
            *r = f(r);
        }
        self.track_filters = std::mem::take(&mut self.track_filters).into_iter().map(|(p, e)| (f(&p), e)).collect();
        self.variant_colors = std::mem::take(&mut self.variant_colors).into_iter().map(|(p, e)| (f(&p), e)).collect();
        for t in &mut self.track_layout {
            t.path = f(&t.path);
        }
    }

    /// Files the session names that aren't there.
    pub fn missing(&self) -> Vec<PathBuf> {
        let mut all: Vec<&PathBuf> = self.bams.iter().chain(&self.vcfs).chain(&self.genes).chain(&self.signals).collect();
        all.extend(&self.reference);
        all.into_iter().filter(|p| !p.exists()).cloned().collect()
    }

    /// Replaces `old` with `new` wherever the session names it.
    pub fn relocate(&mut self, old: &Path, new: &Path) {
        self.map_paths(|p| if p == old { new.to_path_buf() } else { p.to_path_buf() });
    }

    /// Every filter expression parses. A session whose filter doesn't is refused whole rather
    /// than opened unfiltered: dropping the restriction silently changed counts and exported
    /// evidence.
    pub fn validate(&self) -> Result<(), String> {
        if let Some(e) = self.filter.as_ref().and_then(|f| f.expr.as_deref()) {
            crate::expr::intern(e).map_err(|err| format!("filter.expr {e:?}: {err}"))?;
        }
        for (path, e) in &self.track_filters {
            crate::expr::intern(e).map_err(|err| format!("track_filters[{}] {e:?}: {err}", path.display()))?;
        }
        Ok(())
    }

    /// Writes the session atomically (keeping the previous version as a backup), with files
    /// inside the session's folder stored relative to it.
    pub fn save(&self, path: &Path) -> anyhow::Result<()> {
        let mut s = self.clone();
        let dir = path.parent().filter(|d| !d.as_os_str().is_empty()).and_then(|d| std::path::absolute(d).ok());
        if let Some(dir) = dir {
            s.map_paths(|p| {
                let abs = std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf());
                // Relative paths use `/`, which every platform reads, so a project folder moves
                // between Windows, macOS and Linux.
                abs.strip_prefix(&dir).map_or_else(
                    |_| p.to_path_buf(),
                    |r| PathBuf::from(r.components().map(|c| c.as_os_str().to_string_lossy()).collect::<Vec<_>>().join("/")),
                )
            });
        }
        crate::store::write_atomic(path, (serde_json::to_string_pretty(&s)? + "\n").as_bytes())?;
        Ok(())
    }
}

/// Where the last session is kept: `%APPDATA%\ggv\last-session.json` on Windows,
/// `$XDG_CONFIG_HOME/ggv` or `~/.config/ggv` elsewhere.
/// (Settings from before the rename, in the `gv` folder, are copied there first: see
/// [`ggv_io::folders`].)
pub fn last_session_path() -> Option<PathBuf> {
    Some(ggv_io::folders::config_dir()?.join("last-session.json"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A workspace of signal tracks, or just a reference, is a workspace (saved as the last
    /// session).
    #[test]
    fn signals_or_a_reference_alone_are_not_empty() {
        assert!(Session::default().is_empty());
        assert!(!Session { signals: vec!["cov.bw".into()], ..Default::default() }.is_empty());
        assert!(!Session { reference: Some("ref.fa".into()), ..Default::default() }.is_empty());
    }

    #[test]
    fn sessions_with_unparseable_filters_are_refused() {
        let dir = std::env::temp_dir().join(format!("ggv-session-filter-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let ok = Session { version: VERSION, bams: vec!["a.bam".into()], ..Default::default() };
        let mut shared = ok.clone();
        shared.filter = Some(FilterState {
            hide_duplicates: true,
            hide_secondary: true,
            hide_supplementary: false,
            hide_qc_fail: true,
            min_mapq: 0,
            expr: Some("mapq >=".into()),
        });
        let mut per_file = ok.clone();
        per_file.track_filters.insert("a.bam".into(), "mapq >=".into());
        let mut good = ok.clone();
        good.track_filters.insert("a.bam".into(), "hp == 2".into());
        for (name, s, fine) in [("shared", shared, false), ("per-file", per_file, false), ("good", good, true)] {
            let path = dir.join(format!("{name}.ggv.json"));
            s.save(&path).unwrap();
            let loaded = Session::load(&path);
            assert_eq!(loaded.is_ok(), fine, "{name}: {:?}", loaded.err());
            if let Err(e) = loaded {
                assert!(e.to_string().contains("mapq >="), "{name}: names the expression: {e}");
            }
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn style_round_trips_through_json() {
        let style = Style {
            color: ColorBy::BaseMods,
            group: GroupBy::Haplotype,
            sort: Some((43_092_918, SortBy::MappingQuality)),
            soft_clips: true,
            link_splits: true,
            squished: false,
            junctions: true,
            detail: Default::default(),
            coverage_scale: CoverageScale::Log,
            shared_scale: true,
            dark: true,
            filter: ReadFilter { hide_duplicates: false, min_mapq: 20, ..Default::default() },
        };
        let mut s = Session { version: VERSION, bams: vec!["a.bam".into()], locus: Some("chr17:1-100".into()), ..Default::default() };
        s.set_style(style);
        let text = serde_json::to_string(&s).unwrap();
        assert!(text.contains("\"by\":\"MAPQ\"") && text.contains("\"pos\":43092919"), "{text}");
        let back: Session = serde_json::from_str(&text).unwrap();
        assert_eq!(back, s);
        assert_eq!(back.style(), style);
    }

    #[test]
    fn tag_settings_round_trip() {
        let style = Style { color: ColorBy::Tag(*b"CB"), group: GroupBy::Tag(*b"PS"), ..Default::default() };
        let mut s = Session::default();
        s.set_style(style);
        assert_eq!((s.color.as_str(), s.group.as_str()), ("tag:CB", "tag:PS"));
        let back: Session = serde_json::from_str(&serde_json::to_string(&s).unwrap()).unwrap();
        assert_eq!(back.style(), style);
        assert_eq!(style.tag(), Some(*b"PS"));
        // A malformed tag falls back to the default rather than failing the session.
        assert_eq!(parse_color("tag:C"), ColorBy::None);
    }

    #[test]
    fn files_in_the_session_folder_move_with_it() {
        let root = std::env::temp_dir().join(format!("ggv-session-move-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let a = root.join("a");
        std::fs::create_dir_all(a.join("data")).unwrap();
        std::fs::write(a.join("data").join("x.bam"), "").unwrap();
        let outside = std::env::temp_dir().join("ggv-elsewhere-missing.vcf.gz");
        let mut s = Session { version: VERSION, bams: vec![a.join("data").join("x.bam")], vcfs: vec![outside.clone()], ..Default::default() };
        s.track_filters.insert(a.join("data").join("x.bam"), "hp == 1".into());
        s.save(&a.join("s.ggv.json")).unwrap();
        let text = std::fs::read_to_string(a.join("s.ggv.json")).unwrap();
        assert!(text.contains("\"data/x.bam\"") && !text.contains("ggv-session-move"), "stored relative, with /: {text}");
        // The project folder moves; the session still finds its data.
        let b = root.join("b");
        std::fs::rename(&a, &b).unwrap();
        let back = Session::load(&b.join("s.ggv.json")).unwrap();
        assert_eq!(back.bams, vec![b.join("data").join("x.bam")]);
        assert_eq!(back.vcfs, vec![outside.clone()]);
        assert!(back.track_filters.contains_key(&b.join("data").join("x.bam")));
        assert_eq!(back.missing(), vec![outside.clone()]);
        let mut moved = back.clone();
        moved.relocate(&outside, &b.join("data").join("x.bam"));
        assert!(moved.missing().is_empty());
        let _ = std::fs::remove_dir_all(root);
    }

    #[test]
    fn a_damaged_session_loads_from_its_backup() {
        let dir = std::env::temp_dir().join(format!("ggv-session-bak-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("s.ggv.json");
        let s = Session { version: VERSION, locus: Some("chr1:1-100".into()), ..Default::default() };
        s.save(&path).unwrap();
        s.save(&path).unwrap();
        // A crash part-way through some other writer.
        std::fs::write(&path, "{\"version\": 1, \"loc").unwrap();
        let back = Session::load(&path).unwrap();
        assert!(back.recovered);
        assert_eq!(back.locus.as_deref(), Some("chr1:1-100"));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn a_refused_session_is_not_swapped_for_its_backup() {
        let dir = std::env::temp_dir().join(format!("ggv-session-refuse-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("s.ggv.json");
        let s = Session { version: VERSION, locus: Some("chr1:1-100".into()), ..Default::default() };
        s.save(&path).unwrap();
        s.save(&path).unwrap();
        let mut bad = s.clone();
        bad.track_filters.insert("a.bam".into(), "mapq >=".into());
        std::fs::write(&path, serde_json::to_string(&bad).unwrap()).unwrap();
        let e = Session::load(&path).unwrap_err().to_string();
        assert!(e.contains("mapq >="), "refused, naming the filter: {e}");
        let newer = Session { version: VERSION + 1, ..s.clone() };
        std::fs::write(&path, serde_json::to_string(&newer).unwrap()).unwrap();
        assert!(Session::load(&path).unwrap_err().to_string().contains("newer"));
        // A missing file isn't replaced by its backup either.
        std::fs::remove_file(&path).unwrap();
        assert!(Session::load(&path).is_err());
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn missing_and_unknown_fields_default() {
        let s: Session = serde_json::from_str(r#"{"bams": ["x.bam"], "color": "plaid", "future_field": 1}"#).unwrap();
        assert_eq!(s.bams, vec![PathBuf::from("x.bam")]);
        assert_eq!(s.style(), Style::default());
    }
}
