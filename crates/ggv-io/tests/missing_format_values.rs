//! Missing FORMAT values come back as missing, not as htslib's i32::MIN sentinel,
//! from VCF and BCF alike.

use std::path::PathBuf;

use ggv_io::VariantFile;
use ggv_model::Region;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/data/fixtures").join(name)
}

#[test]
fn missing_allele_depths_stay_missing() {
    let path = fixture("missing-fields.vcf.gz");
    assert!(path.exists(), "{} is missing", path.display());
    let mut f = VariantFile::open(&path).unwrap();
    let v = f.query(&Region { contig: "chrQA".into(), start: 140, end: 160 }, false).unwrap();
    assert_eq!(v.len(), 1);
    let [missing, partial] = &v[0].samples[..] else { panic!("two samples") };
    // AD=.
    assert!(missing.ad.is_empty(), "{:?}", missing.ad);
    assert_eq!(missing.ad_text(), None);
    assert_eq!((missing.dp, missing.gq), (None, None));
    // AD=5,.
    assert_eq!(partial.ad, vec![Some(5), None]);
    assert_eq!(partial.ad_text().as_deref(), Some("5,."));
    assert_eq!(serde_json::to_string(&partial.ad).unwrap(), "[5,null]");
}
