//! dtr — te-looker orchestrator CLI. Drop-in for the Pan_TE "te-looker" slot
//! (replaces the legacy repgraph / dtr placeholder).
//!
//! Pan_TE invokes:  dtr run --genome <fa> --out <dir> --threads <n>
//!                      [--dfam-hmm <h> --enable-known-hmm-discovery] [--window-stride <n>] [extra]
//! and requires:    <dir>/families.fasta  (a TE consensus library FASTA for its Combine step).
//!
//! This chains the validated te-looker components in-process by subprocess:
//!   te-discover (A1 extend-to-consensus, replaces cd-hit; tandem split) -> consensi + members.bed
//!   te-refine   (spoa gapped POA refinement of the member spans)        -> families.fasta
//! HMM Track-1 seeding (--dfam-hmm) and graph window-stride are accepted but not yet
//! implemented; they are ignored with a notice (de novo discovery still runs).

use std::env;
use std::fs;
use std::path::PathBuf;
use std::process::{exit, Command};

fn flag_val(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}
fn has_flag(args: &[String], name: &str) -> bool {
    args.iter().any(|a| a == name)
}

/// resolve a sibling binary next to this executable, else fall back to PATH name
fn sibling(name: &str) -> String {
    if let Ok(exe) = env::current_exe() {
        if let Some(dir) = exe.parent() {
            let p = dir.join(name);
            if p.exists() {
                return p.to_string_lossy().into_owned();
            }
        }
    }
    name.to_string()
}

fn run(cmd: &mut Command, what: &str) {
    eprintln!("[dtr] {what}: {cmd:?}");
    let st = cmd.status().unwrap_or_else(|e| {
        eprintln!("[dtr] failed to spawn {what}: {e}");
        exit(1);
    });
    if !st.success() {
        eprintln!("[dtr] {what} failed ({st})");
        exit(1);
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] != "run" {
        eprintln!("usage: dtr run --genome <fa> --out <dir> [--threads N] [--min-count C] [--window-stride N] [...]");
        exit(2);
    }
    let genome = flag_val(&args, "--genome").unwrap_or_else(|| {
        eprintln!("[dtr] --genome required");
        exit(2);
    });
    let out = flag_val(&args, "--out").unwrap_or_else(|| {
        eprintln!("[dtr] --out required");
        exit(2);
    });
    let threads = flag_val(&args, "--threads").unwrap_or_else(|| "1".into());
    // de novo knobs (tunable via --te-looker-extra-args); sensible defaults for the slot
    let min_count = flag_val(&args, "--min-count").unwrap_or_else(|| "20".into());
    let w = flag_val(&args, "--w").unwrap_or_else(|| "1".into()); // minimizer window; 1=exact
                                                                  // accepted-but-not-yet-implemented Pan_TE flags
    if has_flag(&args, "--dfam-hmm") || has_flag(&args, "--enable-known-hmm-discovery") {
        eprintln!("[dtr] note: HMM Track-1 seeding (--dfam-hmm) not yet implemented; running de novo discovery only");
    }
    if let Some(s) = flag_val(&args, "--window-stride") {
        eprintln!("[dtr] note: --window-stride {s} accepted (graph-mode knob); using minimizer w={w} for scaling instead");
    }

    fs::create_dir_all(&out).ok();
    let disc_prefix = PathBuf::from(&out).join("disc");
    let disc = disc_prefix.to_string_lossy().into_owned();
    let members_bed = format!("{disc}.members.bed");
    let families = PathBuf::from(&out).join("families.fasta");

    // Stage A: A1 discovery (te-discover) — interspersed families + member spans (+ tandem track)
    let mut c = Command::new(sibling("te-discover"));
    c.args([
        genome.as_str(),
        "16",
        min_count.as_str(),
        "200",
        disc.as_str(),
        "--w",
        w.as_str(),
    ])
    .env("RAYON_NUM_THREADS", &threads);
    run(&mut c, "discover (A1)");

    // Stage B: spoa gapped POA refinement (te-refine) -> families.fasta
    let families_str = families.to_string_lossy().into_owned();
    if fs::metadata(&members_bed)
        .map(|m| m.len() > 0)
        .unwrap_or(false)
    {
        let mut r = Command::new(sibling("te-refine"));
        r.args([genome.as_str(), members_bed.as_str(), families_str.as_str()]);
        run(&mut r, "refine (spoa)");
    } else {
        eprintln!("[dtr] no interspersed members discovered; writing empty families.fasta");
    }

    // Pan_TE requires families.fasta to exist (empty is tolerated)
    if !families.exists() {
        fs::File::create(&families).ok();
    }
    let n = fs::read_to_string(&families)
        .map(|s| s.matches('>').count())
        .unwrap_or(0);
    eprintln!("[dtr] done: {} families -> {}", n, families.display());
    println!("families\t{n}");
}
