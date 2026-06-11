//! te-seed (Phase 1) — io + canonical k-mer index + seed-stat / seed-dump.
//! Validates the index by reproducing jellyfish's canonical k-mer statistics.

use std::env;
use std::io::{BufWriter, Write};
use std::fs::File;
use te_core::{decode, harvest, read_fasta};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: te-seed <genome.fa> [k=16] [min_count=3] [cap=200] [--seeds out.tsv]");
        std::process::exit(1);
    }
    let fasta = &args[1];
    let k: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(16);
    let min_count: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3);
    let cap: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(200);
    let seeds_out: Option<String> = args.iter().position(|a| a == "--seeds")
        .and_then(|i| args.get(i + 1)).cloned();
    let w: u32 = args.iter().position(|a| a == "--w").and_then(|i| args.get(i + 1)).and_then(|s| s.parse().ok()).unwrap_or(1);

    let t0 = std::time::Instant::now();
    let recs = read_fasta(fasta);
    eprintln!("[io] {} records, {:.1} Mb ({:.1}s)",
        recs.len(), recs.iter().map(|r| r.codes.len()).sum::<usize>() as f64 / 1e6, t0.elapsed().as_secs_f64());

    let (seeds, st) = harvest(&recs, k, w, min_count, cap);
    eprintln!("[harvest+sort] w={w} (1=exact, >1=minimizer), {:.1}s", t0.elapsed().as_secs_f64());

    if let Some(p) = seeds_out {
        let mut w = BufWriter::new(File::create(&p).unwrap());
        for s in &seeds {
            write!(w, "{}\t{}", decode(s.code, k), s.count).unwrap();
            for &(ri, off, strand) in &s.occ {
                write!(w, "\t{}:{}:{}", recs[ri as usize].name, off, if strand { '+' } else { '-' }).unwrap();
            }
            writeln!(w).unwrap();
        }
        w.flush().unwrap();
    }

    println!("k\t{k}");
    println!("Total\t{}", st.total);
    println!("Distinct\t{}", st.distinct);
    println!("Unique(count=1)\t{}", st.unique);
    println!("count=2\t{}", st.count2);
    println!("R(count>=3)\t{}", st.r_ge3);
    println!("Max_count\t{}", st.max_count);
    println!("seeds(count>={min_count})\t{}", seeds.len());
    eprintln!("[done] {:.1}s", t0.elapsed().as_secs_f64());
}
