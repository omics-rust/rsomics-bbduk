// clap --help text names BBDuk / bbduk.sh params that read as code; backticking them clutters the help
#![allow(clippy::doc_markdown)]

mod filter;
mod kmer_set;
mod ref_kmers;

pub use filter::{Config, KTrim, QTrim, is_contaminant, process};
pub use ref_kmers::{MAX_HDIST, MAX_K, RefKmers};
