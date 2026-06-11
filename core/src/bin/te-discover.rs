//! te-discover (Phase 2) — A1 extend-to-growing-consensus discovery engine.
//!
//! Replaces cd-hit all-vs-all clustering with RepeatScout/RECON-style O(occ) family
//! formation: for each high-count seed (greedy, highest first), align each occurrence
//! ONCE against a consensus grown left/right by per-column majority fit; occurrences that
//! keep matching are the family, those that diverge are dropped; the family's genomic
//! spans are masked so later seeds don't rebuild it. Output = a rough consensus library
//! (FASTA) per family + a BED of member spans. The high-quality gapped consensus is then
//! produced by spoa on the member spans (kept as a mature tool — see core/README.md).
//!
//! Cost is O(total_occ x extension_length): each occurrence is touched a bounded number
//! of times, never pairwise — this is the multi-Gb scalability fix (V4_DESIGN §A1).

use std::env;
use std::io::{BufWriter, Write};
use std::fs::File;
use te_core::{harvest, mix64, read_fasta, Record};

// --- tuning (boundary = where flanks diverge to background; member = stays aligned) ---
const AGREE: f64 = 0.60;      // per-column majority agreement to accept a consensus column
const STOP_RUN: u32 = 30;     // stop extending after this many consecutive low-agree columns
const MIN_COV: u32 = 3;       // need >= this many members reaching a column
const MAX_MISM: f64 = 0.35;   // drop a member if mismatch rate exceeds this (after MIN_LEN)
const MIN_LEN: u32 = 40;      // min checked length before a member can be dropped
const MAX_EXT: isize = 6000;  // cap extension each direction (runaway guard)
const MIN_MEMBERS: usize = 3;
const MIN_CONS_LEN: usize = 80;
const BUILD_CAP: usize = 200;  // occurrences used to GROW the consensus (O(occ)); all copies are masked

/// §E1/§E6 runaway guard: short/medium-period self-similarity => tandem/satellite, not a
/// dispersed interspersed element. Routes those to a separate track (and keeps the long
/// runaway consensi out of the interspersed library, where they cripple RepeatMasker).
fn is_tandem(cons: &[u8]) -> bool {
    let n = cons.len();
    if n < 100 { return false; }
    let mut p = 2usize;
    while p <= 600 && p < n / 2 {
        let lim = n - p;
        let mut m = 0usize;
        for i in 0..lim { if cons[i] == cons[i + p] { m += 1; } }
        if (m as f64) / (lim as f64) > 0.60 { return true; }
        p += 1;
    }
    false
}

/// Genomic-tandem guard: a LARGE-period satellite escapes `is_tandem` (its consensus is a
/// single unit with no internal periodicity), but its copies sit TANDEM (adjacent) in the
/// genome, whereas a real interspersed TE's copies are dispersed. If most members have a
/// same-record neighbor within ~1 consensus length, it's a satellite -> tandem track.
fn genomic_tandem(spans: &mut [(u32, i64, i64, bool)], cons_len: i64) -> bool {
    if spans.len() < 3 { return false; }
    spans.sort_unstable();
    let mut adj = 0usize;
    for i in 1..spans.len() {
        if spans[i].0 == spans[i - 1].0 && spans[i].1 - spans[i - 1].2 < cons_len { adj += 1; }
    }
    (adj as f64) / (spans.len() as f64) > 0.5
}

/// Oriented base at consensus-relative position `p` (p=0..k is the seed) for one occurrence.
#[inline]
fn base_at(rec: &Record, off: u32, strand: bool, k: u32, p: isize) -> u8 {
    let len = rec.codes.len() as isize;
    let gi = if strand { off as isize + p } else { off as isize + (k as isize - 1) - p };
    if gi < 0 || gi >= len { return 255; }
    let b = rec.codes[gi as usize];
    if b == 255 { 255 } else if strand { b } else { 3 - b }
}

/// Does this occurrence's region align (ungapped) to the consensus? Used to recruit & mask
/// ALL genuine copies of a family (not just the build sample, and not unrelated seed hits).
fn member_match(rec: &Record, off: u32, strand: bool, cons: &[u8], ll: i64, k: u32) -> bool {
    let (mut m, mut t) = (0u32, 0u32);
    for (ci, &cb) in cons.iter().enumerate() {
        let b = base_at(rec, off, strand, k, (ci as i64 - ll) as isize);
        if b < 4 { t += 1; if b == cb { m += 1; } }
    }
    t > 30 && (m as f64) / (t as f64) >= 0.70
}

/// Extend a consensus in one direction; returns the accepted base codes (in extension order).
/// Updates `active`/`mism`/`elen` (member recruitment) in place.
fn extend(
    recs: &[Record], occ: &[(u32, u32, bool)], k: u32,
    start_p: isize, step: isize,
    active: &mut [bool], mism: &mut [u32], elen: &mut [u32],
) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    let mut last_good = 0usize;
    let mut low_run = 0u32;
    let mut p = start_p;
    while (p - start_p).abs() < MAX_EXT {
        let mut tally = [0u32; 4];
        let mut cov = 0u32;
        for (m, &(ri, off, st)) in occ.iter().enumerate() {
            if !active[m] { continue; }
            let b = base_at(&recs[ri as usize], off, st, k, p);
            if b < 4 { tally[b as usize] += 1; cov += 1; }
        }
        if cov < MIN_COV { break; }
        let mb = (0..4).max_by_key(|&x| tally[x]).unwrap() as u8;
        let agree = tally[mb as usize] as f64 / cov as f64;
        out.push(mb);
        if agree >= AGREE { last_good = out.len(); low_run = 0; }
        else { low_run += 1; if low_run >= STOP_RUN { break; } }
        // member recruitment: penalize mismatch, drop persistent non-members
        for (m, &(ri, off, st)) in occ.iter().enumerate() {
            if !active[m] { continue; }
            let b = base_at(&recs[ri as usize], off, st, k, p);
            if b < 4 {
                if b != mb { mism[m] += 1; }
                elen[m] += 1;
                if elen[m] > MIN_LEN && (mism[m] as f64 / elen[m] as f64) > MAX_MISM {
                    active[m] = false;
                }
            }
        }
        p += step;
    }
    out.truncate(last_good);
    out
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: te-discover <genome.fa> [k=16] [min_count=200] [cap=200] [out_prefix=te_disc]");
        std::process::exit(1);
    }
    let fasta = &args[1];
    let k: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(16);
    let min_count: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(200);
    let cap: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(200);
    let prefix = args.get(5).map(|s| s.as_str()).unwrap_or("te_disc");

    let t0 = std::time::Instant::now();
    let recs = read_fasta(fasta);
    let w: u32 = args.iter().position(|a| a == "--w").and_then(|i| args.get(i + 1)).and_then(|s| s.parse().ok()).unwrap_or(1);
    let _ = cap; // discover keeps FULL occurrences (cap=0) so it can mask ALL copies
    let (mut seeds, _st) = harvest(&recs, k, w, min_count, 0);
    seeds.sort_unstable_by(|a, b| b.count.cmp(&a.count)); // greedy: highest count first
    eprintln!("[seed] {} seeds count>={}, {:.1}s", seeds.len(), min_count, t0.elapsed().as_secs_f64());

    let mut masked: Vec<Vec<bool>> = recs.iter().map(|r| vec![false; r.codes.len()]).collect();
    let mut fa = BufWriter::new(File::create(format!("{prefix}.consensi.fa")).unwrap());
    let mut bed = BufWriter::new(File::create(format!("{prefix}.members.bed")).unwrap());
    let mut tfa = BufWriter::new(File::create(format!("{prefix}.tandem.fa")).unwrap());

    let mut nfam = 0usize;
    let mut n_tandem = 0usize;
    let mut total_align: u64 = 0;
    for s in &seeds {
        // drop occurrences whose seed start is already masked (family likely found)
        let live: Vec<(u32, u32, bool)> = s.occ.iter().cloned()
            .filter(|&(ri, off, _)| !masked[ri as usize][off as usize]).collect();
        if live.len() < MIN_MEMBERS { continue; }

        // grow the consensus from a bounded sample (O(occ)); ALL copies are masked below
        let mut build = live.clone();
        if build.len() > BUILD_CAP {
            build.sort_unstable_by_key(|&(ri, off, _)| mix64(((ri as u64) << 26) | off as u64));
            build.truncate(BUILD_CAP);
        }
        let mut active = vec![true; build.len()];
        let mut mism = vec![0u32; build.len()];
        let mut elen = vec![0u32; build.len()];
        let right = extend(&recs, &build, k, k as isize, 1, &mut active, &mut mism, &mut elen);
        let mut left = extend(&recs, &build, k, -1, -1, &mut active, &mut mism, &mut elen);
        total_align += build.len() as u64; // each sampled member aligned once (O(occ), not O(occ^2))
        if active.iter().filter(|&&a| a).count() < MIN_MEMBERS { continue; }

        // consensus = reverse(left) + seed + right
        left.reverse();
        let mut cons: Vec<u8> = Vec::with_capacity(left.len() + k as usize + right.len());
        cons.extend(&left);
        for i in (0..k).rev() { cons.push(((s.code >> (2 * i)) & 3) as u8); }
        cons.extend(&right);
        if cons.len() < MIN_CONS_LEN { continue; }

        // recruit ALL copies among live that align to the consensus (>=70%), and mask them:
        // masks every genuine copy (so high-copy families aren't rebuilt) without touching
        // unrelated seed hits (which over-masked and suppressed real families).
        let (ll, rl) = (left.len() as i64, right.len() as i64);
        let mut spans: Vec<(u32, i64, i64, bool)> = Vec::new();
        for &(ri, off, st) in &live {
            if !member_match(&recs[ri as usize], off, st, &cons, ll, k) { continue; }
            let (gs, ge) = if st { (off as i64 - ll, off as i64 + k as i64 + rl) }
                           else { (off as i64 - rl, off as i64 + k as i64 + ll) };
            let rlen = recs[ri as usize].codes.len() as i64;
            spans.push((ri, gs.max(0), ge.min(rlen), st));
        }
        if spans.len() < MIN_MEMBERS { continue; }
        let n_mem = spans.len();
        // §E6 routing: internal-period OR genomic-tandem (large-period satellite) -> tandem track
        let tandem = is_tandem(&cons) || genomic_tandem(&mut spans, cons.len() as i64);
        let id = if tandem { n_tandem += 1; format!("tand_{n_tandem}") } else { nfam += 1; format!("fam_{nfam}") };
        let seq: String = cons.iter().map(|&b| b"ACGT"[b as usize] as char).collect();
        {
            let w: &mut dyn Write = if tandem { &mut tfa } else { &mut fa };
            writeln!(w, ">{id} members={n_mem} len={}", cons.len()).unwrap();
            for chunk in seq.as_bytes().chunks(80) { w.write_all(chunk).unwrap(); w.write_all(b"\n").unwrap(); }
        }
        // always mask member spans (so the family isn't rebuilt); BED only for interspersed
        for &(ri, gs, ge, st) in &spans {
            for x in gs..ge { masked[ri as usize][x as usize] = true; }
            if !tandem {
                writeln!(bed, "{}\t{}\t{}\t{}\t.\t{}",
                    recs[ri as usize].name, gs, ge, id, if st { '+' } else { '-' }).unwrap();
            }
        }
    }
    fa.flush().unwrap(); bed.flush().unwrap(); tfa.flush().unwrap();
    eprintln!("[A1] {} interspersed families + {} tandem; {} member alignments (O(occ)); {:.1}s",
        nfam, n_tandem, total_align, t0.elapsed().as_secs_f64());
    println!("interspersed_families\t{nfam}");
    println!("tandem_families\t{n_tandem}");
    println!("member_alignments\t{total_align}");
}
