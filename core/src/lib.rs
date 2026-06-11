//! te-core library — shared io + canonical k-mer index + seed harvest.
//! Used by the `te-seed` (Phase 1 stats) and `te-discover` (Phase 2 A1) binaries.

use std::fs::File;
use std::io::{BufRead, BufReader};

/// A,C,G,T -> 0,1,2,3 (matches the prototype's encoding); everything else invalid (255).
#[inline]
pub fn base_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 255,
    }
}

/// Decode a canonical k-mer code (base0 in the high 2 bits) to ACGT.
pub fn decode(code: u32, k: u32) -> String {
    let mut s = Vec::with_capacity(k as usize);
    for i in (0..k).rev() {
        s.push(b"ACGT"[((code >> (2 * i)) & 3) as usize]);
    }
    String::from_utf8(s).unwrap()
}

/// Deterministic order-independent mixer (hash-bottom-k occurrence capping).
#[inline]
pub fn mix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^ (x >> 33)
}

/// One genome record: name + base codes (0-3, 255 = N/other) for O(1) random access.
pub struct Record {
    pub name: String,
    pub codes: Vec<u8>,
}

pub fn read_fasta(path: &str) -> Vec<Record> {
    let f = File::open(path).unwrap_or_else(|e| panic!("cannot open {path}: {e}"));
    let mut recs: Vec<Record> = Vec::new();
    let mut name = String::new();
    let mut codes: Vec<u8> = Vec::new();
    let mut push = |name: &mut String, codes: &mut Vec<u8>, recs: &mut Vec<Record>| {
        if !name.is_empty() {
            recs.push(Record { name: std::mem::take(name), codes: std::mem::take(codes) });
        }
    };
    for line in BufReader::new(f).lines() {
        let line = line.unwrap();
        if let Some(stripped) = line.strip_prefix('>') {
            push(&mut name, &mut codes, &mut recs);
            name = stripped.split_whitespace().next().unwrap_or("").to_string();
        } else {
            codes.extend(line.trim_end().bytes().map(base_code));
        }
    }
    push(&mut name, &mut codes, &mut recs);
    recs
}

/// A repetitive seed: canonical code + total genome count + (capped) occurrences.
pub struct Seed {
    pub code: u32,
    pub count: u64,
    /// occurrences: (record_index, offset, is_plus_strand)
    pub occ: Vec<(u32, u32, bool)>,
}

/// Harvest canonical (k<=16) k-mers; return seeds with count >= min_count, occurrences
/// hash-bottom-k capped to `cap`. Packs each k-mer as (canon:u32<<32 | gpos:u32), sorts
/// once, then groups — O(N log N), single linear grouping pass. Also returns global stats.
pub struct HarvestStats {
    pub total: u64,
    pub distinct: u64,
    pub unique: u64,
    pub count2: u64,
    pub r_ge3: u64,
    pub max_count: u64,
}

/// Scan one record -> packed (canon:u32<<32 | gpos:u32) entries.
/// w==1: every canonical k-mer (exact, reproduces jellyfish). w>1: (w,k)-minimizers only
/// (1/w density via a monotonic deque) — the memory/throughput lever for multi-Gb genomes.
fn scan_record(rec: &Record, ri: usize, k: u32, w: u32, mask: u32, hi: u32) -> Vec<u64> {
    let mut out: Vec<u64> = Vec::new();
    let (mut fwd, mut rc, mut valid) = (0u32, 0u32, 0u32);
    // deque of (window_pos, hash, packed) for minimizer mode
    let mut dq: std::collections::VecDeque<(u64, u64, u64)> = std::collections::VecDeque::new();
    let mut wpos: u64 = 0;          // index among valid k-mers since last reset
    let mut last_emit: u64 = u64::MAX;
    for (i, &b) in rec.codes.iter().enumerate() {
        if b == 255 { valid = 0; fwd = 0; rc = 0; dq.clear(); wpos = 0; last_emit = u64::MAX; continue; }
        fwd = ((fwd << 2) | b as u32) & mask;
        rc = (rc >> 2) | (((3 - b) as u32) << hi);
        valid += 1;
        if valid < k { continue; }
        let canon = fwd.min(rc);
        let strand = (fwd <= rc) as u32;
        let off = (i + 1 - k as usize) as u32;
        let packed = ((canon as u64) << 32) | ((strand << 31) | ((ri as u32) << 25) | (off & 0x1FF_FFFF)) as u64;
        if w <= 1 { out.push(packed); continue; }
        let h = mix64(canon as u64);
        while matches!(dq.back(), Some(&(_, bh, _)) if bh >= h) { dq.pop_back(); }
        dq.push_back((wpos, h, packed));
        while matches!(dq.front(), Some(&(p, _, _)) if p + w as u64 <= wpos + 1) { dq.pop_front(); }
        if wpos + 1 >= w as u64 {
            let (mp, _, mpacked) = *dq.front().unwrap();
            if mp != last_emit { out.push(mpacked); last_emit = mp; }
        }
        wpos += 1;
    }
    out
}

pub fn harvest(recs: &[Record], k: u32, w: u32, min_count: u64, cap: usize) -> (Vec<Seed>, HarvestStats) {
    use rayon::prelude::*;
    assert!(k >= 1 && k <= 16, "build supports 1<=k<=16");
    let mask: u32 = if k == 16 { u32::MAX } else { (1u32 << (2 * k)) - 1 };
    let hi = 2 * (k - 1);
    // parallel per-record scan, then concatenate (O(N)) and parallel sort
    let chunks: Vec<Vec<u64>> = recs.par_iter().enumerate()
        .map(|(ri, rec)| scan_record(rec, ri, k, w, mask, hi))
        .collect();
    let mut packed: Vec<u64> = chunks.concat();
    let total = packed.len() as u64;
    packed.par_sort_unstable();

    let mut seeds = Vec::new();
    let mut st = HarvestStats { total, distinct: 0, unique: 0, count2: 0, r_ge3: 0, max_count: 0 };
    let n = packed.len();
    let mut i = 0;
    while i < n {
        let canon = (packed[i] >> 32) as u32;
        let mut j = i + 1;
        while j < n && (packed[j] >> 32) as u32 == canon { j += 1; }
        let count = (j - i) as u64;
        st.distinct += 1;
        match count { 1 => st.unique += 1, 2 => st.count2 += 1, _ => st.r_ge3 += 1 }
        if count > st.max_count { st.max_count = count; }
        if count >= min_count {
            let mut g: Vec<u32> = packed[i..j].iter().map(|&p| p as u32).collect();
            if cap > 0 && g.len() > cap {   // cap==0 => full occurrence list (te-discover masks all copies)
                g.sort_unstable_by_key(|&x| mix64(x as u64));
                g.truncate(cap);
            }
            let occ = g.iter().map(|&x| ((x >> 25) & 0x3F, x & 0x1FF_FFFF, (x >> 31) & 1 == 1)).collect();
            seeds.push(Seed { code: canon, count, occ });
        }
        i = j;
    }
    (seeds, st)
}
