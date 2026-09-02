# Bacterial 16S rRNA Sequencing Pipeline

A standalone Python script that takes a folder of paired Sanger `.ab1`
chromatograms (27F forward / 1492R reverse reads) and produces cleaned
sequences, a multiple sequence alignment, identity-based clustering, and a
phylogenetic tree.

## Setup

```bash
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

```bash
python bacterial_16s_pipeline.py --input /path/to/ab1_folder --output /path/to/results_folder
```

`--input` should be a folder containing `.ab1` files named like
`<sample>_27-F.ab1` and `<sample>_149-R.ab1` (adjust the direction-detection
logic in `process_ab1_file()` if your filenames use a different pattern,
e.g. `1492-R` instead of `149-R`).

## Outputs (written to `--output`)

| File | Description |
|---|---|
| `01_trim_qc_results.json` | Per-read quality-trim + pass/fail QC log |
| `02_consensus.json` | Per-sample consensus assembly log |
| `cleaned_sequences.fasta` | One cleaned sequence per sample that passed QC |
| `aligned.fasta` | Multiple sequence alignment of the cleaned sequences |
| `identity_clusters.json` | Pairwise %identity clustering (97% threshold) |
| `phylo_tree_nj.nwk` | Neighbor-Joining tree, Newick format |
| `phylo_tree.png` | Rendered tree image, colored by cluster |

## Pipeline steps

1. **Quality trimming** — parses each `.ab1` with Biopython, extracts
   per-base Phred quality scores, and trims each read to its best
   contiguous window using a modified Mott algorithm (Kadane's
   maximum-subarray algorithm on `quality - 20`).
2. **QC filtering** — a trimmed read must be ≥350bp and average ≥Q20 to pass;
   otherwise it's excluded. A sample is dropped entirely if both its reads fail.
3. **Consensus assembly** — where both reads pass, the reverse read is
   reverse-complemented and locally aligned against the forward read; if they
   overlap at ≥97% identity over ≥50bp, they're merged into one extended
   consensus sequence. Otherwise the better single read is used.
4. **Multiple sequence alignment** — all cleaned sequences aligned with FAMSA.
5. **Identity clustering** — pairwise %identity computed directly from the
   alignment; single-linkage clustering at 97% groups likely-same-species
   samples together, cross-checked against GC content.
6. **Phylogenetic tree** — Neighbor-Joining tree built from an
   identity-based distance matrix, midpoint-rooted, rendered as a colored PNG.

## Notes / things you may want to change

- **Species/genus identification is NOT included.** This pipeline clusters
  sequences by similarity but does not call BLAST or any reference database.
  To identify species, take representative sequences from
  `cleaned_sequences.fasta` (or the whole file) and BLAST them at
  https://blast.ncbi.nlm.nih.gov/Blast.cgi against the "16S ribosomal RNA
  sequences (Bacteria and Archaea)" database. This can also be automated with
  Biopython's `Bio.Blast.NCBIWWW.qblast()` if you have normal internet access.
- **Thresholds are adjustable** — QC length/quality cutoffs, overlap
  identity/length cutoffs, and the clustering identity threshold are all
  named constants near the top of the script.
- **Primer sequences** (`PRIMER_27F`, `PRIMER_1492R`) assume standard
  universal bacterial 16S primers; change these if you used different primers.
