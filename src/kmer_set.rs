use std::hash::{BuildHasherDefault, Hasher};

// k-mer codes are 2-bit encodings — already well-distributed; SipHash is wasted work
#[derive(Default)]
pub(crate) struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    fn write_u64(&mut self, n: u64) {
        self.0 = n;
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!();
    }
    fn finish(&self) -> u64 {
        self.0
    }
}

pub(crate) type KmerHashSet = std::collections::HashSet<u64, BuildHasherDefault<IdentityHasher>>;

// Sorted-array probe: ≤ SMALL_THRESHOLD unique k-mers fit in a few cache lines;
// binary search beats HashSet bucket chasing at this scale.
const SMALL_THRESHOLD: usize = 64;

// Adaptive k-mer set: sorted array for small refs, HashSet for large ones.
// Both hold deduplicated, fully-expanded (hdist variants + RC) k-mer codes.
pub(crate) enum KmerSet {
    Small(Box<[u64]>),
    Large(KmerHashSet),
}

impl KmerSet {
    pub(crate) fn from_vec(mut v: Vec<u64>) -> Self {
        v.sort_unstable();
        v.dedup();
        if v.len() <= SMALL_THRESHOLD {
            Self::Small(v.into_boxed_slice())
        } else {
            let mut s = KmerHashSet::default();
            s.extend(v);
            Self::Large(s)
        }
    }

    pub(crate) fn contains(&self, k: &u64) -> bool {
        match self {
            Self::Small(a) => a.binary_search(k).is_ok(),
            Self::Large(s) => s.contains(k),
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Small(a) => a.len(),
            Self::Large(s) => s.len(),
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
