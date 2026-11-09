//! Which transcripts to name at a change, and which of their consequences sums it up. Used by
//! `ggv batch`; the MCP tools keep their own copy for now and could switch to this.

use ggv_annot::{TAG_CANONICAL, TAG_MANE_PLUS_CLINICAL, TAG_MANE_SELECT, Transcript, coding};

/// Transcripts named at one change unless more genes are there.
pub const MAX_TRANSCRIPTS: usize = 4;

/// Whether `[start, end)` takes coding bases of `t` (not just its span: an intron or UTR of a
/// coding transcript doesn't count).
fn coding_at(t: &Transcript, start: u64, end: u64) -> bool {
    let Some((cs, ce)) = t.cds else { return false };
    let end = end.max(start + 1);
    t.exons.iter().any(|&(s, e)| s.max(cs).max(start) < e.min(ce).min(end))
}

/// How well `t` speaks for its gene at a change over `[start, end)`, best highest: the gene's
/// main transcript; MANE Select, MANE Plus Clinical, canonical (as ggv-annot ranks them for the
/// main one); curated over predicted (`XM_`); one the change is coding in over one where it's
/// intronic or UTR; then the longest coding sequence.
fn rank(t: &Transcript, start: u64, end: u64) -> (bool, u8, bool, bool, bool, u64) {
    let tag = if t.tags & TAG_MANE_SELECT != 0 {
        3
    } else if t.tags & TAG_MANE_PLUS_CLINICAL != 0 {
        2
    } else {
        u8::from(t.tags & TAG_CANONICAL != 0)
    };
    (t.representative, tag, !ggv_annot::is_predicted(&t.id), coding_at(t, start, end), t.cds.is_some(), t.cds_len())
}

/// The transcripts to name at a change over `[start, end)`, from those overlapping it: every
/// gene there once first (its best by [`rank`]), then the others best first up to
/// [`MAX_TRANSCRIPTS`]. Isoforms of one gene, predicted models first, used to crowd out an
/// antisense or nested gene, and a gene's first isoform in the file to stand for it where the
/// change is intronic in it but coding in another.
pub fn pick(ts: Vec<&Transcript>, start: u64, end: u64) -> Vec<&Transcript> {
    let mut ts = ts;
    // Stable: equal ranks keep the annotation's order.
    ts.sort_by_key(|t| std::cmp::Reverse(rank(t, start, end)));
    let mut picks: Vec<&Transcript> = Vec::new();
    for &t in &ts {
        if !picks.iter().any(|p| p.gene_id == t.gene_id) {
            picks.push(t);
        }
    }
    // Coding genes' main transcripts before the rest.
    picks.sort_by_key(|t| (t.cds.is_none(), !t.representative));
    let n = picks.len().max(MAX_TRANSCRIPTS);
    for t in ts {
        if picks.len() >= n {
            break;
        }
        if !picks.iter().any(|p| std::ptr::eq(*p, t)) {
            picks.push(t);
        }
    }
    picks
}

/// How severe an effect is (as [`coding::consequence`] names it), highest worst, so a place is
/// summarised by what its change does most, not by whichever overlapping gene comes first:
/// 6 truncating, start-, stop- or splice-site-changing; 5 missense or in-frame; 4 coding with
/// the protein change not worked out; 3 splice region; 2 synonymous; 1 UTR; 0 intron, non-coding
/// exon and anything else (ggv-annot names effects as text, with no ranked type).
pub fn severity(effect: &str) -> u8 {
    match effect {
        "stop gained" | "frameshift" | "start lost" | "stop lost" | "splice donor" | "splice acceptor" => 6,
        "missense" => 5,
        e if e.starts_with("inframe") => 5,
        e if e == coding::UNCERTAIN || e.starts_with("coding") || e.starts_with("transcript boundary") => 4,
        "splice region" => 3,
        "synonymous" => 2,
        "5' UTR" | "3' UTR" => 1,
        _ => 0,
    }
}

/// Which transcript's consequence wins a tie in [`severity`]: the gene's main one, then MANE.
pub fn main_rank(t: &Transcript) -> (bool, bool) {
    (t.representative, t.tags & (TAG_MANE_SELECT | TAG_MANE_PLUS_CLINICAL) != 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(id: &str, gene: &str, exons: &[(u64, u64)], cds: Option<(u64, u64)>) -> Transcript {
        Transcript {
            contig: 0,
            start: exons.first().unwrap().0,
            end: exons.last().unwrap().1,
            strand: ggv_annot::Strand::Forward,
            id: id.into(),
            gene_id: gene.into(),
            gene_name: gene.into(),
            biotype: "".into(),
            exons: exons.to_vec(),
            cds,
            tags: 0,
            representative: false,
            translation: Default::default(),
            exons_from_span: false,
        }
    }

    #[test]
    fn a_gene_is_named_by_an_isoform_the_change_is_coding_in() {
        // bughunt4 (GNAS XLas): the main transcript doesn't reach; NESP55 (first in the file,
        // intronic here) stood for the gene, and XM_ models filled the rest.
        let nesp = tx("NM_016592.5", "G", &[(0, 100), (900, 1000)], Some((50, 950)));
        let xm = tx("XM_1.1", "G", &[(0, 100), (900, 1000)], Some((50, 950)));
        let short = tx("NM_2.1", "G", &[(0, 100), (400, 600), (900, 1000)], Some((450, 950)));
        let xl = tx("NM_080425.4", "G", &[(0, 100), (300, 700), (900, 1000)], Some((350, 950)));
        let picks = pick(vec![&nesp, &xm, &short, &xl], 500, 501);
        let ids: Vec<&str> = picks.iter().map(|t| &*t.id).collect();
        assert_eq!(ids, ["NM_080425.4", "NM_2.1", "NM_016592.5", "XM_1.1"]);
        // MANE Plus Clinical outranks a longer curated one.
        let mut plus = short.clone();
        plus.tags = TAG_MANE_PLUS_CLINICAL;
        assert_eq!(&*pick(vec![&nesp, &xl, &plus], 500, 501)[0].id, "NM_2.1");
        // The main transcript, wherever it is.
        let mut main = nesp.clone();
        main.representative = true;
        assert_eq!(&*pick(vec![&xl, &main], 500, 501)[0].id, "NM_016592.5");
    }

    #[test]
    fn severity_orders_effects() {
        let order = ["intron", "3' UTR", "synonymous", "splice region", coding::UNCERTAIN, "missense", "stop gained"];
        assert!(order.windows(2).all(|w| severity(w[0]) < severity(w[1])), "{order:?}");
        assert_eq!(severity("inframe deletion"), severity("missense"));
        assert_eq!(severity("frameshift"), severity("splice donor"));
        assert_eq!(severity("non-coding exon"), severity("intron"));
    }
}
