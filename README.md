# te-looker — de novo detection of divergent transposable elements

`te-looker` is a de novo (reference-library-free) transposable-element / repeat detector. Its
command-line tool, **`dtr`** (Dark TE Rebuilder), discovers repeat families directly from a
genome FASTA and writes a consensus library — targeting the older, divergent elements that
exact-match library searches miss. It is the Step-4 discovery track of the
[Pan_TE](https://github.com/unavailable-2374/Pan_TE) pan-genome TE-annotation pipeline.

The method was developed and validated against the *Arabidopsis thaliana* T2T Col-CEN genome,
where it raised RepeatMasker masking coverage from the `mdl-repeat` baseline of **22.15%** to
**~27.0%** (**+~4.85 pp**, every increment copy-number-validated — see [Validated results](#validated-results)).

## Install

**From source** (needs a Rust toolchain and [`spoa`](https://github.com/rvaser/spoa) on `PATH`):

```bash
cd core
cargo build --release
# binaries: core/target/release/{dtr,te-discover,te-refine,te-seed}
```

`dtr` locates its `te-discover` / `te-refine` helpers next to its own executable, so keep the
four binaries in the same directory (or install them all onto `PATH`).

`spoa` is resolved from `PATH` by default; override with `SPOA_BIN=/path/to/spoa` if needed.

Within Pan_TE, `dtr` and `spoa` are provided by the conda packaging (`te-looker` + `spoa`).

## Usage

```bash
dtr run --genome genome.fa --out out_dir [--threads N] [--min-count C] [--w W]
```

- `--genome`  input genome FASTA (required)
- `--out`     output directory (required); the library is written to `out_dir/families.fasta`
- `--threads` worker threads (default 1)
- `--min-count` minimum occurrences of a canonical 16-mer to seed a family (default 20; lower
  for small genomes / test data)
- `--w`       minimizer window (default 1 = exact k-mers; raise to scale to large genomes)

Pipeline:

```
dtr run
  └─ te-discover : A1 extend-to-consensus seeding (O(occ) linear kernel) → consensi + member spans + tandem track
  └─ te-refine   : spoa gapped-POA refinement of each family's members   → out_dir/families.fasta
```

## Validated results

Validation uses **RepeatMasker masking-coverage increment** as a non-circular metric (rather
than matching against a possibly noisy library). **Every increment passes a real-genome copy-number
completeness gate (≥3 copies).**

| Step | Script | Masked | Δ |
|---|---|---|---|
| `mdl-repeat` baseline | — | 22.15% | — |
| k-mer pass 1 (count ≥200, residual-targeted) | `proto_scaled.py` | 24.32% | +2.17 |
| k-mer iterated (count ≥20 + low-complexity filter + copy validation) | `proto_scaled2.py` | 25.52% | +1.20 |
| protein-anchored ancient TE — consensus (diamond × RepeatPeps) | `proto_protein.py` | 25.99% | +0.47 |
| protein-anchored ancient TE — element-span footprint | (bed workflow) | 26.22% | +0.23 |
| TRF tandem (chr3+5 extrapolated) | (TRF + bedtools) | ~27.0% | +0.80 |
| RepeatMasker `-s` / naive boundary extension | `proto_refine.py` | — | ~0 / −5.17 (falsified) |

Two complementary discovery channels — **exact k-mers for high-copy families (+3.37 pp)** and
**protein anchoring for 70–90%-divergent ancient elements (+0.70 pp)** — attack orthogonal element
classes. This reaches the estimated true-repeat ceiling of ~27.3–28%; the remaining gap sits behind
the completeness wall (1–2-copy domesticated genes / <70% identity), where pushing further trades
away precision.

## Repository layout

- **`core/`** — the Rust crate (`te-core`) and the usable CLI tools (`dtr`, `te-discover`,
  `te-refine`, `te-seed`). This is the production code.
- **`src/`** — the Python research prototypes used to develop and validate the method, in order
  of method evolution. These are experimental, demo-bound reproduction scripts (they hardcode
  *A. thaliana* demo paths and intermediate products), **not** the tool. See the table below.
- **`docs/`** — design and review records.

### Research prototypes (`src/`)

| Script | Role |
|---|---|
| `proto_a1.py` | Stage 1–2 seeding + A1 extend-to-consensus (O(occ) linear kernel) |
| `proto_stage4.py` | + Stage-4 cd-hit cross-seed clustering + spoa POA consensus + boundaries |
| `proto_e3e4.py` | + E4 adaptive window + E3 SINE attempt (failed) + E1 runaway guard |
| `proto_scaled.py` | residual-targeted whole-genome discovery on `mdl-repeat`-missed regions (count ≥200) |
| `proto_scaled2.py` | iterated threshold (count ≥20) + homopolymer/low-complexity seed filter + min-copy gate |
| `proto_refine.py` | naive boundary extension — a backfire case (chimeras, lower masking) |
| `proto_refgate.py` | Refiner-style per-family acceptance gate (accept only if refinement covers more) |
| `proto_protein.py` | diamond blastx × RepeatPeps anchoring of ancient TEs → cluster → POA |

External tools used by the prototypes: jellyfish, cd-hit-est, spoa, minimap2, blastn, diamond,
RepeatMasker, TRF, RepeatProteinMask.

### Documentation (`docs/`)

- `NEXTGEN_DESIGN.md` — the original design proposal that was reviewed.
- `V4_DESIGN.md` — the revised v4 design (A1 kernel, organelle guard, Helitron, adaptive window),
  each item annotated with its empirical-validation status.
- `DESIGN_REVIEW.md` — the main record: dual-perspective review, gating experiments, three-tier
  algorithmic prototyping, and the full masking-improvement trail (non-circular metric, copy-number
  completeness gate, discovery directions #1–#5).

## Reusable methodology

1. **Non-circular metric** — RepeatMasker masking-coverage increment instead of matching a noisy library.
2. **Real copy-number completeness gate** at every step — sample newly masked spans, count genomic
   copies by blastn, require **≥3 copies** to credit a gain. This gate caught a naive iteration's
   spurious +6 pp (60% single-copy), a `-s` false gain, and the −5.17 pp boundary-extension backfire.
3. **Two orthogonal channels** — exact k-mer (high-copy) and protein-anchored (divergent ancient),
   with non-overlapping footprints.
4. **Residual targeting** — discover on `mdl-repeat`-unmasked regions to aim directly at the gap.
