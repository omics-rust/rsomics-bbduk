use std::collections::HashMap;

use rsomics_kmer::{RollingKmers, encode, reverse_complement};

use crate::filter::Config;
use crate::kmer_set::KmerSet;

// k ≤ 31: reverse_complement does 1u64 << (2*k), UB at k=32; BBDuk's own default is 27
pub const MAX_K: usize = 31;

// variant index grows (3k)^hdist; past 3 it outgrows memory — fail loud, not OOM
pub const MAX_HDIST: usize = 3;

// zero the middle base's 2-bit field (BBDuk maskmiddle); for odd k the centre is
// revcomp-invariant, so masking at insert and probe is symmetric
#[must_use]
pub(crate) fn mask_mid(code: u64, k: usize) -> u64 {
    let m = k / 2;
    let shift = 2 * (k - 1 - m);
    code & !(0b11u64 << shift)
}

// non-ACGT window → nothing: BBDuk drops ambiguous ref k-mers; expanding one would
// fabricate reference content
pub(crate) fn variants(w: &[u8], hdist: usize) -> Vec<u64> {
    if w.iter().any(|b| !matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return Vec::new();
    }
    let mut buf = w.to_vec();
    let mut out = Vec::new();
    expand(&mut buf, 0, hdist, &mut out);
    out
}

pub(crate) fn expand(buf: &mut Vec<u8>, start: usize, budget: usize, out: &mut Vec<u64>) {
    if let Ok(code) = encode(buf) {
        out.push(code);
    }
    if budget == 0 {
        return;
    }
    for pos in start..buf.len() {
        let orig = buf[pos];
        for &b in b"ACGT" {
            if b != orig {
                buf[pos] = b;
                expand(buf, pos + 1, budget - 1, out);
            }
        }
        buf[pos] = orig;
    }
}

// Exact ref k-mer index (not a Bloom filter — false positives would break kfilter
// byte-compat); mink>0 also indexes shorter tip k-mers (lengths mink..k) governed by hdist2
pub struct RefKmers {
    pub(crate) full: KmerSet,
    pub(crate) tips: HashMap<usize, KmerSet>, // length → ref k-mers of that length (mink..k)
    pub(crate) k: usize,
    pub(crate) mink: usize,
    // effective maskmiddle = requested && mink==0 — BBDuk disables maskmiddle when short
    // k-mers are in play
    pub(crate) mm: bool,
}

impl RefKmers {
    #[must_use]
    pub fn build<'a>(refs: impl IntoIterator<Item = &'a [u8]>, cfg: &Config) -> Self {
        assert!(
            (1..=MAX_K).contains(&cfg.k),
            "k must be in 1..={MAX_K} (2-bit codec / rcomp limit); got {}",
            cfg.k
        );
        assert!(
            cfg.hdist <= MAX_HDIST && cfg.hdist2 <= MAX_HDIST,
            "hdist/hdist2 must be ≤ {MAX_HDIST} (variant index memory bound)"
        );
        assert!(cfg.mink <= cfg.k, "mink must be ≤ k");

        let mm = cfg.maskmiddle && cfg.mink == 0;
        let key = |c: u64| if mm { mask_mid(c, cfg.k) } else { c };

        let mut full_vec: Vec<u64> = Vec::new();
        let mut tip_vecs: HashMap<usize, Vec<u64>> = HashMap::new();
        let tip_lens: Vec<usize> = if cfg.mink > 0 {
            (cfg.mink..cfg.k).collect()
        } else {
            Vec::new()
        };
        for seq in refs {
            for w in seq.windows(cfg.k) {
                for code in variants(w, cfg.hdist) {
                    full_vec.push(key(code));
                    if cfg.rcomp {
                        full_vec.push(key(reverse_complement(code, cfg.k)));
                    }
                }
            }
            for &l in &tip_lens {
                if seq.len() < l {
                    continue;
                }
                let v = tip_vecs.entry(l).or_default();
                for w in seq.windows(l) {
                    for code in variants(w, cfg.hdist2) {
                        v.push(code);
                        if cfg.rcomp {
                            v.push(reverse_complement(code, l));
                        }
                    }
                }
            }
        }
        let full = KmerSet::from_vec(full_vec);
        let tips: HashMap<usize, KmerSet> = tip_vecs
            .into_iter()
            .map(|(l, v)| (l, KmerSet::from_vec(v)))
            .collect();
        Self {
            full,
            tips,
            k: cfg.k,
            mink: cfg.mink,
            mm,
        }
    }

    fn key(&self, c: u64) -> u64 {
        if self.mm { mask_mid(c, self.k) } else { c }
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.full.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.full.is_empty()
    }

    #[must_use]
    pub fn hits(&self, seq: &[u8]) -> usize {
        if seq.len() < self.k {
            return 0;
        }
        RollingKmers::new(seq, self.k)
            .flatten()
            .filter(|&c| self.full.contains(&self.key(c)))
            .count()
    }

    #[must_use]
    pub fn first_full_hit(&self, seq: &[u8]) -> Option<usize> {
        if seq.len() < self.k {
            return None;
        }
        RollingKmers::new(seq, self.k)
            .enumerate()
            .find_map(|(i, opt)| {
                opt.filter(|&c| self.full.contains(&self.key(c)))
                    .map(|_| i + 1 - self.k)
            })
    }

    // ktrim=r cut: full-k-mer hit, else (mink>0) start of the longest ref-matching 3' tip —
    // BBDuk never trims an internal short k-mer
    #[must_use]
    pub fn right_cut(&self, seq: &[u8]) -> Option<usize> {
        if let Some(p) = self.first_full_hit(seq) {
            return Some(p);
        }
        if self.mink == 0 {
            return None;
        }
        // when the read is shorter than k it is itself the longest tip, so the bound is
        // inclusive: min(k-1, seq.len())
        for l in (self.mink..=(self.k - 1).min(seq.len())).rev() {
            let tip = &seq[seq.len() - l..];
            if self
                .tips
                .get(&l)
                .is_some_and(|s| encode(tip).is_ok_and(|c| s.contains(&c)))
            {
                return Some(seq.len() - l);
            }
        }
        None
    }

    // ktrim=l keep-from: end of full-k-mer hit, else (mink>0) end of the longest ref-matching
    // 5' tip
    #[must_use]
    pub fn left_cut(&self, seq: &[u8]) -> Option<usize> {
        if let Some(p) = self.first_full_hit(seq) {
            return Some(p + self.k);
        }
        if self.mink == 0 {
            return None;
        }
        for l in (self.mink..=(self.k - 1).min(seq.len())).rev() {
            let tip = &seq[..l];
            if self
                .tips
                .get(&l)
                .is_some_and(|s| encode(tip).is_ok_and(|c| s.contains(&c)))
            {
                return Some(l);
            }
        }
        None
    }
}
