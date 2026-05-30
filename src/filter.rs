use crate::ref_kmers::RefKmers;

// BBDuk ktrim: None=f (kfilter, a hit removes the read), Right=r (trim hit→3'), Left=l (trim 5'→hit end)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KTrim {
    None,
    Right,
    Left,
}

// BBDuk qtrim=f|r|l|rl
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QTrim {
    None,
    Right,
    Left,
    Both,
}

/// BBDuk processing configuration.
#[derive(Debug, Clone)]
pub struct Config {
    pub k: usize,
    pub mink: usize,
    pub hdist: usize,
    pub hdist2: usize,
    pub rcomp: bool,
    pub maskmiddle: bool,
    pub min_kmer_hits: usize,
    pub min_kmer_fraction: f64,
    pub ktrim: KTrim,
    pub qtrim: QTrim,
    pub trimq: u8,
    pub min_length: usize,
    pub qual_offset: u8,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            k: 27,
            mink: 0,
            hdist: 0,
            hdist2: 0,
            rcomp: true,
            maskmiddle: true,
            min_kmer_hits: 1,
            min_kmer_fraction: 0.0,
            ktrim: KTrim::None,
            qtrim: QTrim::None,
            trimq: 6,
            min_length: 10,
            qual_offset: 33,
        }
    }
}

/// BBDuk kfilter contaminant test.
#[must_use]
pub fn is_contaminant(seq: &[u8], refs: &RefKmers, cfg: &Config) -> bool {
    let hits = refs.hits(seq);
    if hits == 0 {
        return false;
    }
    if hits >= cfg.min_kmer_hits {
        return true;
    }
    if cfg.min_kmer_fraction > 0.0 {
        let total = seq.len().saturating_sub(cfg.k) + 1;
        if total > 0 {
            #[allow(clippy::cast_precision_loss)]
            let frac = hits as f64 / total as f64;
            return frac >= cfg.min_kmer_fraction;
        }
    }
    false
}

/// BBDuk pipeline order: k-mer filter / k-mer trim → quality trim → min-length.
///
/// Returns `Some((start, end))` for the kept slice, `None` if the read is discarded.
#[must_use]
pub fn process(seq: &[u8], qual: &[u8], refs: &RefKmers, cfg: &Config) -> Option<(usize, usize)> {
    let mut start = 0usize;
    let mut end = seq.len();

    if !refs.is_empty() {
        match cfg.ktrim {
            KTrim::None => {
                if is_contaminant(seq, refs, cfg) {
                    return None;
                }
            }
            KTrim::Right => {
                if let Some(pos) = refs.right_cut(seq) {
                    end = pos;
                }
            }
            KTrim::Left => {
                if let Some(pos) = refs.left_cut(seq) {
                    start = pos.min(end);
                }
            }
        }
    }

    if end > start && cfg.qtrim != QTrim::None {
        let thr = cfg.trimq;
        let phred = |q: u8| q.saturating_sub(cfg.qual_offset);
        if matches!(cfg.qtrim, QTrim::Right | QTrim::Both) {
            while end > start && phred(qual[end - 1]) < thr {
                end -= 1;
            }
        }
        if matches!(cfg.qtrim, QTrim::Left | QTrim::Both) {
            while start < end && phred(qual[start]) < thr {
                start += 1;
            }
        }
    }

    if end.saturating_sub(start) < cfg.min_length {
        return None;
    }
    Some((start, end))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ref_kmers::RefKmers;

    fn refs(seqs: &[&[u8]], cfg: &Config) -> RefKmers {
        RefKmers::build(seqs.iter().copied(), cfg)
    }

    // Exact-match property tests disable maskmiddle so the centre base is not wildcarded.
    fn exact(k: usize) -> Config {
        Config {
            k,
            maskmiddle: false,
            ..Config::default()
        }
    }

    #[test]
    fn contaminant_read_with_shared_kmer_is_filtered() {
        let cfg = exact(11);
        let r = refs(&[b"ACGTTGCAACGTTGCA"], &cfg);
        assert!(is_contaminant(b"TTTTACGTTGCAACGTTGCATTTT", &r, &cfg));
        assert!(!is_contaminant(b"TTTTTTTTTTTTTTTTTTTTTTTT", &r, &cfg));
    }

    #[test]
    fn rcomp_match_is_found_when_enabled() {
        let fwd: &[u8] = b"ACGTTGCAACG";
        let on = exact(11);
        let off = Config {
            rcomp: false,
            ..on.clone()
        };
        let probe: Vec<u8> = fwd
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                _ => b'C',
            })
            .collect();
        assert!(is_contaminant(&probe, &refs(&[fwd], &on), &on));
        assert!(!is_contaminant(&probe, &refs(&[fwd], &off), &off));
    }

    #[test]
    fn hdist_one_matches_single_mismatch_kmer() {
        let cfg = Config {
            hdist: 1,
            rcomp: false,
            ..exact(11)
        };
        let r = refs(&[b"ACGTTGCAACG"], &cfg);
        // index 1 is not the middle (k=11 centre is index 5), so hdist=0 must reject it.
        assert!(is_contaminant(b"ATGTTGCAACG", &r, &cfg));
        let strict = Config {
            hdist: 0,
            ..cfg.clone()
        };
        assert!(!is_contaminant(
            b"ATGTTGCAACG",
            &refs(&[b"ACGTTGCAACG"], &strict),
            &strict
        ));
    }

    #[test]
    fn maskmiddle_default_wildcards_the_centre_base() {
        // k=11 ⇒ middle index 5. A single mismatch *at the middle* is a
        // hit with maskmiddle (the BBDuk default), a miss without it.
        let mm = Config {
            k: 11,
            ..Config::default()
        };
        let r = refs(&[b"ACGTTGCAACG"], &mm);
        assert!(is_contaminant(b"ACGTTACAACG", &r, &mm)); // index-5 G→A
        let no_mm = exact(11);
        assert!(!is_contaminant(
            b"ACGTTACAACG",
            &refs(&[b"ACGTTGCAACG"], &no_mm),
            &no_mm
        ));
    }

    #[test]
    fn mink_disables_maskmiddle() {
        // mink>0 overrides the defaulted-true maskmiddle.
        let cfg = Config {
            k: 11,
            mink: 6,
            ..Config::default()
        };
        let r = refs(&[b"ACGTTGCAACG"], &cfg);
        assert!(!is_contaminant(b"ACGTTACAACG", &r, &cfg));
    }

    #[test]
    fn ktrim_right_trims_from_adapter_kmer() {
        let cfg = Config {
            ktrim: KTrim::Right,
            rcomp: false,
            ..exact(12)
        };
        let adapter = b"AGATCGGAAGAGC";
        let r = refs(&[adapter], &cfg);
        let read = b"ACGTTGCAACGTACGTTGCAAGATCGGAAGAGC";
        let qual = vec![b'I'; read.len()];
        let (s, e) = process(read, &qual, &r, &cfg).unwrap();
        assert_eq!(s, 0);
        assert_eq!(e, 20);
    }

    #[test]
    fn mink_trims_partial_adapter_at_read_tip() {
        let cfg = Config {
            k: 12,
            mink: 6,
            ktrim: KTrim::Right,
            rcomp: false,
            ..Config::default()
        };
        let adapter = b"AGATCGGAAGAGCACAC";
        let r = refs(&[adapter], &cfg);
        // only the first 7 bp of the adapter are present (< k=12, ≥ mink=6).
        let read = b"ACGTTGCAACGTACGTTGCAAGATCGG";
        let qual = vec![b'I'; read.len()];
        let (s, e) = process(read, &qual, &r, &cfg).unwrap();
        assert_eq!(s, 0);
        assert_eq!(e, 20);
    }

    #[test]
    fn short_read_whole_length_tip_is_checked() {
        // The tip-length upper bound is inclusive: a read shorter than k
        // whose entire body is the adapter tip must still clip.
        let cfg = Config {
            k: 12,
            mink: 6,
            ktrim: KTrim::Right,
            rcomp: false,
            ..Config::default()
        };
        let adapter = b"AGATCGGAAGAGCACAC";
        let r = refs(&[adapter], &cfg);
        let read = b"AGATCGG";
        let qual = vec![b'I'; read.len()];
        // entirely adapter ⇒ trimmed to length 0 ⇒ discarded (< minlength).
        assert!(process(read, &qual, &r, &cfg).is_none());
    }

    #[test]
    fn variants_skip_ambiguous_ref_window() {
        // An N in the ref window must not fabricate ACGT variants.
        let cfg = Config {
            hdist: 1,
            rcomp: false,
            ..exact(11)
        };
        let r = refs(&[b"ACGTNGCAACG"], &cfg);
        assert!(r.is_empty(), "N-bearing ref window must index nothing");
    }

    #[test]
    fn qtrim_right_and_minlength_discard() {
        let cfg = Config {
            qtrim: QTrim::Right,
            trimq: 20,
            min_length: 10,
            ..Config::default()
        };
        let empty = RefKmers::build(std::iter::empty(), &cfg);
        let seq = b"ACGTACGTACGTACGTACGT";
        let mut q = vec![b'I'; 8];
        q.extend(std::iter::repeat_n(b'#', 12));
        assert!(process(seq, &q, &empty, &cfg).is_none());
    }

    #[test]
    fn clean_read_survives_unchanged() {
        let cfg = Config::default();
        let empty = RefKmers::build(std::iter::empty(), &cfg);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let q = vec![b'I'; seq.len()];
        assert_eq!(process(seq, &q, &empty, &cfg), Some((0, seq.len())));
    }
}
