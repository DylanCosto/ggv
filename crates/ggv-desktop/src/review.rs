//! Variant triage: step through a VCF's sites and mark each real, artifact or unsure, with an
//! optional note. Decisions are saved on every change to `<vcf>.ggv-review.tsv` next to the VCF
//! (or, when that folder can't be written, elsewhere: see `store`), so a review can be stopped
//! and resumed.

use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Label {
    Real,
    Artifact,
    Unsure,
}

impl Label {
    pub fn name(self) -> &'static str {
        match self {
            Self::Real => "real",
            Self::Artifact => "artifact",
            Self::Unsure => "unsure",
        }
    }

    fn parse(s: &str) -> Option<Self> {
        [Self::Real, Self::Artifact, Self::Unsure].into_iter().find(|l| l.name() == s)
    }

    pub fn rgb(self) -> [u8; 3] {
        match self {
            Self::Real => [30, 160, 60],
            Self::Artifact => [210, 50, 50],
            Self::Unsure => [230, 170, 20],
        }
    }
}

/// A site: contig, 1-based position, REF, ALT (comma-joined).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SiteKey {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Decision {
    pub label: Label,
    pub note: String,
}

pub struct Review {
    /// Id of the variant track being reviewed (its index changes when tracks move).
    pub track: u64,
    pub path: PathBuf,
    /// The VCF being reviewed.
    pub source: PathBuf,
    pub decisions: BTreeMap<SiteKey, Decision>,
    /// Sites decided here not saved yet (the save failed): kept over what is on disk.
    unsaved: BTreeSet<SiteKey>,
    /// The saves of `path` in this process the decisions include (see `refresh`).
    seen: u64,
    /// The note being typed for the current site.
    pub note: String,
    /// The site `note` was loaded for: arriving at another site loads its saved note.
    note_for: Option<SiteKey>,
    pub error: Option<String>,
    /// Where the decisions went when they couldn't be saved where they were.
    pub notice: Option<String>,
    /// Which of the records starting at the outlined base the review is on (several can):
    /// contig index, 0-based base, and its place among them (`usize::MAX`: the last).
    pub record: Option<(usize, u64, usize)>,
    /// A Prev under way, which arrives at the last record at its base.
    pub back: bool,
}

pub const EXT: &str = "ggv-review.tsv";

const HEADER: &str = "#contig\tpos\tref\talt\tdecision\tnote";

impl Review {
    /// Starts (or resumes, if the file exists) a review of variant track `track` (its id) at `vcf`.
    pub fn open(track: u64, vcf: &Path) -> Self {
        let path = crate::store::sidecar(vcf, EXT);
        let seen = crate::store::saves(&path);
        let mut review = Self {
            track,
            path,
            source: vcf.to_path_buf(),
            decisions: BTreeMap::new(),
            unsaved: BTreeSet::new(),
            seen,
            note: String::new(),
            note_for: None,
            error: None,
            notice: None,
            record: None,
            back: false,
        };
        if review.path.exists() {
            match crate::store::read_recovering(&review.path, parse) {
                Ok((d, recovered)) => {
                    review.decisions = d;
                    if recovered {
                        review.notice = Some(format!("{} was damaged; restored the previous save", review.path.display()));
                    }
                }
                Err(e) => review.error = Some(format!("{} can't be read, nor its backup: {e}", review.path.display())),
            }
        }
        review
    }

    pub fn counts(&self) -> [usize; 3] {
        let mut c = [0; 3];
        for d in self.decisions.values() {
            c[d.label as usize] += 1;
        }
        c
    }

    /// Puts `key`'s saved note (or nothing) in the editor when the review arrives at a new site,
    /// so re-marking a reviewed site keeps its note rather than writing "" over it, and a note
    /// typed for one site doesn't follow onto the next (as for queues).
    pub fn arrive(&mut self, key: &SiteKey) {
        if self.note_for.as_ref() != Some(key) {
            self.note = self.decisions.get(key).map(|d| d.note.clone()).unwrap_or_default();
            self.note_for = Some(key.clone());
        }
    }

    pub fn set(&mut self, key: SiteKey, label: Label) {
        // Kept in the editor: marking the same site again (changing your mind) keeps the note.
        let note = self.note.clone();
        self.decisions.insert(key.clone(), Decision { label, note });
        self.unsaved.insert(key);
        // Another pane reviewing this VCF may have saved since: the save takes its decisions
        // in, and this one, the latest, goes over them.
        self.relocate();
        self.save(self.path.clone());
    }

    /// Takes in decisions another pane reviewing the same VCF saved since this one last read
    /// or saved the file (also when its saves had to move elsewhere).
    pub fn refresh(&mut self) {
        if self.relocate() || crate::store::saves(&self.path) != self.seen {
            self.merge_saved();
        }
    }

    /// Follows the VCF's decisions to where they were last saved (by another pane, after a save
    /// had to go elsewhere); whether they moved.
    fn relocate(&mut self) -> bool {
        match crate::store::located(&self.source, EXT) {
            Some(p) if p != self.path => {
                self.path = p;
                true
            }
            _ => false,
        }
    }

    /// The saved decisions taken in (see `take_in`).
    fn merge_saved(&mut self) {
        self.seen = crate::store::saves(&self.path);
        let path = self.path.clone();
        self.take_in(&path, false);
    }

    /// Takes in the decisions saved at `path` ([`crate::store::merge_decisions`]; a file
    /// that's gone or can't be read keeps what this review has). What it has wins over a
    /// damaged file's older backup, and with `moving` (`path` isn't where they were read).
    /// Whether the file is damaged.
    fn take_in(&mut self, path: &Path, moving: bool) -> bool {
        if !path.exists() {
            return false;
        }
        match crate::store::read_recovering(path, parse) {
            Ok((saved, recovered)) => {
                self.decisions = crate::store::merge_decisions(&self.decisions, &self.unsaved, saved, recovered || moving);
                recovered
            }
            Err(_) => false,
        }
    }

    /// Saves the decisions, taking in what is saved where they go first (`from`: where they
    /// were read).
    fn save(&mut self, from: PathBuf) {
        let mut path = std::mem::take(&mut self.path);
        let source = self.source.clone();
        let result = crate::store::save_sidecar(&mut path, &source, EXT, |p| {
            let damaged = self.take_in(p, p != from);
            (format(&self.decisions).into_bytes(), damaged)
        });
        self.path = path;
        match result {
            Ok(moved) => {
                self.error = None;
                self.unsaved.clear();
                self.seen = crate::store::saved(&self.path);
                if moved.is_some() {
                    self.notice = moved;
                }
            }
            Err(e) => self.error = Some(e),
        }
    }

    /// Moves the decisions to where new ones now go (after a folder for reviews was chosen).
    pub fn move_to_preferred(&mut self) {
        let old = std::mem::replace(&mut self.path, crate::store::preferred(&self.source, EXT));
        if old != self.path {
            self.save(old);
            if self.error.is_none() {
                self.notice = Some(format!("Decisions now saved to {}", self.path.display()));
            }
        }
    }
}

fn clean(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

fn format(decisions: &BTreeMap<SiteKey, Decision>) -> String {
    let mut out = String::from(HEADER);
    out.push('\n');
    for (k, d) in decisions {
        out.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\n", k.contig, k.pos, k.reference, k.alt, d.label.name(), clean(&d.note)));
    }
    out
}

/// A saved review: empty, or the header then one decision per line (blank lines and `#`
/// comments are fine). Anything else is damage, so the backup is tried.
fn parse(text: &str) -> anyhow::Result<BTreeMap<SiteKey, Decision>> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let mut out = BTreeMap::new();
    let Some(first) = lines.next() else { return Ok(out) };
    anyhow::ensure!(first.trim_end() == HEADER, "not a saved variant review");
    for l in lines.filter(|l| !l.starts_with('#')) {
        let row = || {
            let f: Vec<&str> = l.trim_end_matches('\r').split('\t').collect();
            let (&[contig, pos, reference, alt, label] | &[contig, pos, reference, alt, label, _]) = &f[..] else { return None };
            if contig.is_empty() || reference.is_empty() {
                return None;
            }
            let key = SiteKey { contig: contig.to_owned(), pos: pos.parse().ok()?, reference: reference.to_owned(), alt: alt.to_owned() };
            Some((key, Decision { label: Label::parse(label)?, note: f.get(5).unwrap_or(&"").to_string() }))
        };
        let (k, d) = row().ok_or_else(|| anyhow::anyhow!("damaged line: {l}"))?;
        out.insert(k, d);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(pos: u64) -> SiteKey {
        SiteKey { contig: "chr20".into(), pos, reference: "G".into(), alt: "C".into() }
    }

    #[test]
    fn decisions_round_trip_through_tsv() {
        let mut d = BTreeMap::new();
        d.insert(key(1_895_892), Decision { label: Label::Real, note: "clean het\tboth strands".into() });
        d.insert(key(100), Decision { label: Label::Artifact, note: String::new() });
        let text = format(&d);
        assert!(text.starts_with(HEADER));
        assert!(text.contains("chr20\t1895892\tG\tC\treal\tclean het both strands"), "tabs in notes are cleaned: {text}");
        let back = parse(&text).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back[&key(100)].label, Label::Artifact);
        assert_eq!(back[&key(1_895_892)].note, "clean het both strands");
    }

    #[test]
    fn resumes_from_file_and_saves_on_each_decision() {
        let dir = std::env::temp_dir().join(format!("ggv-review-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("calls.vcf.gz");
        let mut r = Review::open(0, &vcf);
        assert!(r.decisions.is_empty());
        r.arrive(&key(5));
        r.note = "strand bias".into();
        r.set(key(5), Label::Unsure);
        r.arrive(&key(6));
        assert!(r.note.is_empty(), "a note stays with its site");
        r.set(key(6), Label::Real);
        let mut again = Review::open(0, &vcf);
        assert_eq!(again.counts(), [1, 0, 1]);
        assert_eq!(again.decisions[&key(5)].note, "strand bias");
        // Reopened, arriving at a reviewed site and marking it again keeps its note.
        again.arrive(&key(5));
        assert_eq!(again.note, "strand bias");
        again.set(key(5), Label::Real);
        again.set(key(5), Label::Artifact);
        assert_eq!(Review::open(0, &vcf).decisions[&key(5)].note, "strand bias");
        assert_eq!(again.path, dir.join("calls.vcf.gz.ggv-review.tsv"));
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Decisions saved before the rename to GGV (`calls.vcf.gz.gv-review.tsv`) are read, and
    /// kept being saved there, while there is no file under the new name; one under the new
    /// name wins.
    #[test]
    fn decisions_from_before_the_rename_are_kept() {
        let dir = std::env::temp_dir().join(format!("ggv-review-legacy-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("calls.vcf.gz");
        let old = dir.join("calls.vcf.gz.gv-review.tsv");
        std::fs::write(&old, format!("{HEADER}
chr20	5	G	C	artifact	strand bias
")).unwrap();
        let mut r = Review::open(0, &vcf);
        assert_eq!(r.path, old);
        assert_eq!(r.counts(), [0, 1, 0]);
        assert_eq!(r.decisions[&key(5)].note, "strand bias");
        r.set(key(6), Label::Real);
        assert!(!dir.join("calls.vcf.gz.ggv-review.tsv").exists());
        assert_eq!(Review::open(0, &vcf).counts(), [1, 1, 0]);
        // Exports must not overwrite it either.
        assert!(crate::store::sidecars(&vcf, EXT).contains(&old));
        // With both, the file under the new name is the one used.
        let both = dir.join("both.vcf.gz");
        std::fs::write(dir.join("both.vcf.gz.gv-review.tsv"), format!("{HEADER}
chr20	5	G	C	artifact	
")).unwrap();
        std::fs::write(dir.join("both.vcf.gz.ggv-review.tsv"), format!("{HEADER}
chr20	7	G	C	unsure	
")).unwrap();
        let r = Review::open(0, &both);
        assert_eq!(r.path, dir.join("both.vcf.gz.ggv-review.tsv"));
        assert_eq!(r.counts(), [0, 0, 1]);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// One VCF reviewed in two panes: each pane's save keeps the other's decisions (the latest
    /// for a site wins), and each sees the other's once refreshed.
    #[test]
    fn two_reviews_of_one_vcf_keep_each_others_decisions() {
        let dir = std::env::temp_dir().join(format!("ggv-review-two-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("calls.vcf.gz");
        let (mut one, mut two) = (Review::open(0, &vcf), Review::open(1, &vcf));
        one.set(key(1), Label::Real);
        two.set(key(2), Label::Artifact);
        one.set(key(2), Label::Unsure);
        assert_eq!(one.counts(), [1, 0, 1]);
        two.refresh();
        assert_eq!(two.counts(), [1, 0, 1], "the other pane's decisions, the latest for each site");
        two.set(key(3), Label::Artifact);
        assert_eq!(Review::open(0, &vcf).counts(), [1, 1, 1]);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A VCF whose folder can't be written (its "folder" is a file), reviewed in two panes
    /// opened before any decision: the saves go to GGV's own folder, and both panes follow them
    /// there, so neither wipes the other's decisions.
    #[test]
    fn two_reviews_of_a_vcf_in_an_unwritable_folder_keep_each_others_decisions() {
        let dir = std::env::temp_dir().join(format!("ggv-review-readonly-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        crate::store::TEST_CONFIG.set(Some(dir.join("config")));
        std::fs::write(dir.join("share"), "not a folder").unwrap();
        let vcf = dir.join("share").join("calls.vcf.gz");
        let (mut one, mut two) = (Review::open(0, &vcf), Review::open(1, &vcf));
        assert_eq!(one.path, two.path);
        one.set(key(1), Label::Real);
        one.set(key(2), Label::Real);
        assert!(one.path.starts_with(dir.join("config")), "{}", one.path.display());
        assert!(one.notice.as_deref().is_some_and(|n| n.contains("now saved to")), "{:?}", one.notice);
        two.refresh();
        assert_eq!((two.counts(), &two.path), ([2, 0, 0], &one.path), "the other pane follows the decisions");
        two.set(key(3), Label::Artifact);
        one.refresh();
        assert_eq!(one.counts(), [2, 1, 0]);
        assert_eq!(Review::open(0, &vcf).counts(), [2, 1, 0]);
        // A pane saving where decisions already are (another process, an earlier session)
        // takes them in first.
        std::fs::write(dir.join("other"), "not a folder").unwrap();
        let other = dir.join("other").join("calls.vcf.gz");
        let mut three = Review::open(2, &other);
        // The last place under the current name (the ones after it are for decisions from before the rename).
        let fallback = crate::store::sidecars(&other, EXT).into_iter().rfind(|p| p.to_string_lossy().ends_with(EXT)).unwrap();
        std::fs::create_dir_all(fallback.parent().unwrap()).unwrap();
        std::fs::write(&fallback, format!("{HEADER}\nchr20\t9\tG\tC\tunsure\tearlier\n")).unwrap();
        three.set(key(1), Label::Real);
        assert_eq!(three.path, fallback);
        assert_eq!(three.counts(), [1, 0, 1]);
        assert_eq!(Review::open(2, &other).decisions[&key(9)].note, "earlier");
        crate::store::TEST_CONFIG.set(None);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// The decisions file damaged while reviewing: the next decision keeps every decision made
    /// (not the backup's older set), and the good backup isn't replaced by the damaged file.
    #[test]
    fn a_damaged_save_keeps_the_latest_decisions() {
        let dir = std::env::temp_dir().join(format!("ggv-review-damaged-live-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("calls.vcf.gz");
        let mut r = Review::open(0, &vcf);
        r.set(key(1), Label::Real);
        r.set(key(2), Label::Real);
        // Changing a decision: the backup has the old one.
        r.note = "important".into();
        r.set(key(2), Label::Artifact);
        let bak = crate::store::backup_of(&r.path);
        std::fs::write(&r.path, "garbage\n").unwrap();
        r.note.clear();
        r.set(key(4), Label::Unsure);
        assert_eq!(r.counts(), [1, 1, 1]);
        let again = Review::open(0, &vcf);
        assert_eq!((again.counts(), again.notice.clone()), ([1, 1, 1], None));
        assert_eq!(again.decisions[&key(2)].note, "important");
        let backup = std::fs::read_to_string(&bak).unwrap();
        assert!(!backup.contains("garbage") && parse(&backup).is_ok_and(|d| d.len() == 2), "{backup}");
        // Refreshing from a damaged file keeps them too.
        std::fs::write(&r.path, "garbage\n").unwrap();
        let mut two = Review::open(1, &vcf);
        r.set(key(5), Label::Real);
        two.refresh();
        assert_eq!(two.counts(), [2, 1, 1]);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn damaged_saves_fall_back_to_the_backup() {
        let dir = std::env::temp_dir().join(format!("ggv-review-damaged-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("calls.bcf");
        let mut r = Review::open(7, &vcf);
        r.note = "first".into();
        r.set(key(5), Label::Real);
        r.note = String::new();
        r.set(key(6), Label::Artifact);
        let bak = crate::store::backup_of(&r.path);
        for damaged in [
            "not a valid saved variant review".to_owned(),
            HEADER[..9].to_owned(),
            format!("{HEADER}\nchr20\t5\tG\tC\tre"),
            format!("{HEADER}\nchr20\t5\tG"),
            format!("{HEADER}\nchr20\tfive\tG\tC\treal\t"),
        ] {
            std::fs::write(&r.path, &damaged).unwrap();
            let again = Review::open(7, &vcf);
            assert_eq!(again.decisions.get(&key(5)), Some(&Decision { label: Label::Real, note: "first".into() }), "{damaged:?}");
            assert!(again.notice.as_deref().is_some_and(|n| n.contains("damaged")), "{damaged:?}");
            assert_eq!(again.error, None);
        }
        // Blank lines, comments, CRLF and an empty file are all fine.
        std::fs::write(&r.path, format!("\r\n{HEADER}\r\n# note\r\nchr20\t5\tG\tC\tunsure\t\r\n\r\n")).unwrap();
        let again = Review::open(7, &vcf);
        assert_eq!((again.counts(), again.notice.clone()), ([0, 0, 1], None));
        std::fs::write(&r.path, "").unwrap();
        assert_eq!(Review::open(7, &vcf).counts(), [0, 0, 0]);
        // Neither readable: an error, not a fresh review.
        std::fs::write(&r.path, "junk").unwrap();
        std::fs::write(&bak, "junk").unwrap();
        let again = Review::open(7, &vcf);
        assert!(again.error.as_deref().is_some_and(|e| e.contains("can't be read")), "{:?}", again.error);
        let _ = std::fs::remove_dir_all(dir);
    }
}
