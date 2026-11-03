# soil-metagenome-assembly

Bash pipeline for **hybrid metagenomic assembly of soil samples** with Illumina paired-end reads and PacBio long reads, built around **metaSPAdes** and designed to fit high-complexity soil data on a single server with ~250 GB RAM.

```
Illumina R1/R2 ──► fastp ──► tadpole ──► bbnorm ──► metaSPAdes (hybrid) ──► contigs / scaffolds
                  (optional)  (error     (digital        ▲
                              correction) normalization)  │
PacBio reads ────────────────────────────────────────────┘
```

## Why these steps

Soil metagenomes contain billions of distinct k-mers, and metaSPAdes needs memory roughly proportional to that number. To make the assembly fit in memory, the pipeline:

1. **Trims** reads with fastp (skip with `-S` if reads are already trimmed).
2. **Corrects sequencing errors** with BBTools `tadpole.sh`, which needs far less memory than SPAdes' built-in BayesHammer.
3. **Normalizes coverage** with BBTools `bbnorm.sh`: caps highly abundant organisms at a target depth and removes reads with very low apparent depth.
4. **Assembles** with `metaspades.py --only-assembler`, using the PacBio reads for gap closing and repeat resolution.

## Installation

```bash
mamba env create -f environment.yml     # or: conda env create -f environment.yml
conda activate soil_asm
```

Tested with fastp 1.4.0, BBTools 40.05 and SPAdes 4.3.0. GNU `time` (`/usr/bin/time`) is optional; if present, the pipeline records peak RAM for every step.

## Usage

### One sample

```bash
bash scripts/soil_assembly.sh \
  -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
  -l sample.hifi_reads.fastq \
  -t 32 -m 220 -n 50 \
  -o sample_assembly
```

| Option | Description | Default |
|---|---|---|
| `-1`, `-2` | Illumina paired-end reads (required) | |
| `-l` | Long reads (FASTQ or FASTQ.gz) | none |
| `-L` | Long-read type: `clr`, `hifi` or `nanopore`. `hifi` also runs a HiFi-only assembly with metaMDBG | `clr` |
| `-S` | Skip fastp (reads already trimmed) | off |
| `-o` | Output directory | `soil_assembly` |
| `-t` | Threads | 32 |
| `-m` | Memory limit for metaSPAdes, in GB | 220 |
| `-T` | Temporary directory (use a fast disk) | `OUTDIR/tmp` |
| `-k` | k-mer sizes for metaSPAdes | `21,33,55,77` |
| `-n` | bbnorm target depth | 50 |

### Many samples

`scripts/run_all.sh` runs the pipeline sample by sample, pairing every *N* Illumina samples with one PacBio sample (by default 3 Illumina samples per PacBio sample, e.g. `Met_19`, `Met_20`, `Met_21` → `Met_7`). Each Illumina sample is assembled separately with its group's PacBio reads.

Edit the settings at the top of the script (input directories, output directory, sample numbering and the `r1_path` / `r2_path` / `pacbio_path` functions for your file names), then:

```bash
bash scripts/run_all.sh -n    # dry run: shows the pairs and checks that all files exist
bash scripts/run_all.sh       # run
```

The batch script waits for any metaSPAdes run that is already in progress, continues with the next sample if one fails, writes a summary to `OUT_BASE/batch_summary.tsv` and, with `CLEANUP=1`, deletes large intermediate files after each successful assembly.

### Comparing assemblies

```bash
bash scripts/assembly_stats.sh -m 1000 \
  sampleA/05_metaspades/contigs.fasta \
  other_assembler/final.contigs.fa.gz
```

Prints contig count, total length, largest contig, N50, L50, N90, contigs ≥ 1/5/10 kb and GC% side by side. Use `-m` to apply the same minimum contig length to every assembly (assemblers differ in the shortest contig they report), and `-t` for TSV output.

## Output

```
OUTDIR/
├── 01_fastp/          fastp reports (or links to the input reads with -S)
├── 02_tadpole/        error-corrected reads
├── 03_bbnorm/         normalized reads and k-mer histograms (khist_before/after.txt)
├── 04_metamdbg/       HiFi-only assembly (only with -L hifi)
├── 05_metaspades/     contigs.fasta, scaffolds.fasta, assembly graph, spades.log
├── logs/              one log per step, plus run time and peak RAM (*.time.txt)
└── assembly_stats.txt
```

For binning, map the **trimmed, non-normalized** reads to the contigs: normalization distorts abundance.

## Resuming

Every completed step is recorded in `OUTDIR/.checkpoints/`. If a run stops (crash, reboot, closed terminal), run the same command again: completed steps are skipped and an unfinished metaSPAdes run continues from its last checkpoint (`--restart-from last`). To redo a step, delete its `.done` file and its outputs.

## Memory notes

These lessons come from a ~39 Gbp soil sample (after trimming) on a 251 GB RAM server:

- **Keep `-m` well below physical RAM.** With `-m 245` on a 251 GB machine and only 1 GB of swap, the server froze and rebooted while metaSPAdes built the K21 graph. With `-m 220`, SPAdes stops with an error instead if it runs out of memory, and the system stays responsive.
- **The K21 graph construction is the memory peak.** In this sample it reached ~187 GB (~7 billion distinct 21-mers, ~1.9 billion graph edges) during `Connecting the graph`. Watch the third column of `spades.log` (peak memory) during this stage.
- **Java uses 15–25 GB more than its heap.** The BBTools heap is set to 70% of `-m`.
- **Adding swap** (e.g. a 64 GB swap file) is a cheap safety net against the out-of-memory killer.
- If memory is still insufficient, the most effective lever is a higher bbnorm `min` depth (`NORM_MIN` in the script; default 2), at the cost of discarding reads from very rare organisms. Alternatively, use MEGAHIT (`--presets meta-large`), which needs much less memory.

Typical run time for the sample above with 32 threads: fastp ~5 min, tadpole ~20 min, bbnorm ~1 h, metaSPAdes ~1–1.5 days (K21 ~5 h).

## Citation

If you use this pipeline, please cite the underlying tools: fastp, BBTools (BBMap), SPAdes / metaSPAdes and, if used, metaMDBG.
