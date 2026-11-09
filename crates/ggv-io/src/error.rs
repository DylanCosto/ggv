use std::{io, path::PathBuf};

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{path}: {cause}")]
    Io { path: PathBuf, cause: io::Error },

    #[error("{path}: no index found (looked for {tried}); {how}")]
    MissingIndex { path: PathBuf, tried: String, how: String },

    #[error("invalid locus `{text}`: {cause}")]
    InvalidLocus { text: String, cause: ggv_model::LocusParseError },

    #[error("contig `{name}` is not in this file (it has e.g. {})", examples.join(", "))]
    UnknownContig { name: String, examples: Vec<String> },

    #[error("{locus} starts past the end of {contig} (length {length})")]
    OutOfBounds { locus: String, contig: String, length: u64 },
}

pub(crate) trait IoContext<T> {
    fn at(self, path: &std::path::Path) -> Result<T>;
}

impl<T> IoContext<T> for io::Result<T> {
    fn at(self, path: &std::path::Path) -> Result<T> {
        self.map_err(|cause| Error::Io { path: path.to_owned(), cause })
    }
}
