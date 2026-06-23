//! te-refine (Phase 2b) — gapped POA refinement of A1 families via spoa.
//!
//! "Build where tools fail (A1 clustering), reuse where they excel (POA)": the A1 engine
//! (te-discover) does the scalable O(occ) family formation and emits member spans
//! (members.bed); this step hands each family's member sequences to spoa for the
//! high-quality gapped consensus that the ungapped A1 consensus can't match.

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use te_core::read_fasta;

/// Resolve the `spoa` binary: honor `$SPOA_BIN` / `$SPOA`, else find `spoa` on `PATH`
/// (the conda package puts it there). No host-specific path is baked in.
fn spoa_bin() -> String {
    env::var("SPOA_BIN")
        .or_else(|_| env::var("SPOA"))
        .unwrap_or_else(|_| "spoa".to_string())
}

const MAX_MEMBERS: usize = 30;
const MIN_MEMBERS: usize = 3;

#[inline]
fn code_to_base(c: u8) -> u8 {
    if c < 4 {
        b"ACGT"[c as usize]
    } else {
        b'N'
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: te-refine <genome.fa> <members.bed> <out.fa>");
        std::process::exit(1);
    }
    let t0 = std::time::Instant::now();
    let recs = read_fasta(&args[1]);
    let idx: HashMap<&str, usize> = recs
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.as_str(), i))
        .collect();

    // group members.bed by family: name -> Vec<(rec, start, end, plus)>
    let mut fams: HashMap<String, Vec<(usize, usize, usize, bool)>> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    for line in BufReader::new(File::open(&args[2]).unwrap()).lines() {
        let l = line.unwrap();
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 6 {
            continue;
        }
        let ri = match idx.get(f[0]) {
            Some(&i) => i,
            None => continue,
        };
        let (s, e) = (
            f[1].parse::<usize>().unwrap(),
            f[2].parse::<usize>().unwrap(),
        );
        let fam = f[3].to_string();
        let plus = f[5] == "+";
        fams.entry(fam.clone())
            .or_insert_with(|| {
                order.push(fam.clone());
                Vec::new()
            })
            .push((ri, s, e, plus));
    }

    let mut out = BufWriter::new(File::create(&args[3]).unwrap());
    let mut nf = 0usize;
    let tmp = format!("{}.spoa_in.fa", args[3]);
    let spoa = spoa_bin();
    for fam in &order {
        let mut mem = fams[fam].clone();
        if mem.len() < MIN_MEMBERS {
            continue;
        }
        if mem.len() > MAX_MEMBERS {
            mem.truncate(MAX_MEMBERS);
        } // bounded per-family POA
          // write oriented member sequences for spoa
        {
            let mut w = BufWriter::new(File::create(&tmp).unwrap());
            for (i, &(ri, s, e, plus)) in mem.iter().enumerate() {
                let seq: Vec<u8> = recs[ri].codes[s..e.min(recs[ri].codes.len())]
                    .iter()
                    .map(|&c| code_to_base(c))
                    .collect();
                let seq = if plus { seq } else { revcomp(&seq) };
                writeln!(w, ">{i}").unwrap();
                w.write_all(&seq).unwrap();
                w.write_all(b"\n").unwrap();
            }
            w.flush().unwrap();
        }
        // spoa -r 0 = consensus only
        let o = Command::new(&spoa)
            .args(["-r", "0", &tmp])
            .output()
            .unwrap_or_else(|e| {
                eprintln!(
                    "[refine] failed to run spoa ('{spoa}'): {e}\n  \
                     install spoa (the conda package provides it) or set SPOA_BIN=/path/to/spoa"
                );
                std::process::exit(3);
            });
        let cons: String = String::from_utf8_lossy(&o.stdout)
            .lines()
            .filter(|l| !l.starts_with('>'))
            .collect();
        if cons.len() < 80 {
            continue;
        }
        nf += 1;
        writeln!(out, ">{fam}_ref members={} len={}", mem.len(), cons.len()).unwrap();
        for chunk in cons.as_bytes().chunks(80) {
            out.write_all(chunk).unwrap();
            out.write_all(b"\n").unwrap();
        }
    }
    out.flush().unwrap();
    let _ = std::fs::remove_file(&tmp);
    eprintln!(
        "[refine] {} families refined via spoa; {:.1}s",
        nf,
        t0.elapsed().as_secs_f64()
    );
    println!("refined\t{nf}");
}
