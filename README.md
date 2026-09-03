# Bacterial 16S rRNA Sequencing Pipeline
[![Python 3.14.7](https://img.shields.io/badge/Python-3.14.7-blue.svg)](https://www.python.org)
[![VS Code](https://img.shields.io/badge/VS--Code-green.svg)](https://code.visualstudio.com)

A Python script that takes paired Sanger `.ab1` chromatograms and produces cleaned sequences, a multiple sequence alignment, identity-based clustering, and a phylogenetic tree.

## Setup on Mac

```bash
python3 -m venv venv
source venv/bin/activate       
pip install -r requirements.txt
```

## Usage

```bash
python bacterial_16s_pipeline.py --input /path/to/ab1_folder --output /path/to/results_folder
```

`--input` should be a folder containing `.ab1` files.
adjust the direction-detection logic in `process_ab1_file()` if your filenames use a different pattern.


## Pipeline steps

1. **Quality trimming**: parses each `.ab1` with Biopython, extracts
   per-base Phred quality scores, and trims each read to its best
   contiguous window using a modified Mott algorithm (Kadane's
   maximum-subarray algorithm on `quality - 20`).
2. **QC filtering**: a trimmed read must be ≥350bp and average ≥Q20 to pass;
   otherwise it's excluded. A sample is dropped entirely if both its reads fail.
3. **Consensus assembly**: where both reads pass, the reverse read is
   reverse-complemented and locally aligned against the forward read; if they
   overlap at ≥97% identity over ≥50bp, they're merged into one extended
   consensus sequence. Otherwise the better single read is used.
4. **Multiple sequence alignment**: all cleaned sequences aligned with FAMSA.
5. **Identity clustering**: pairwise %identity computed directly from the
   alignment; single-linkage clustering at 97% groups likely-same-species
   samples together, cross-checked against GC content.
6. **Phylogenetic tree**: Neighbor-Joining tree built from an
   identity-based distance matrix, midpoint-rooted, rendered as a colored PNG.

## Notes

- **Species/genus identification is NOT included**: This pipeline clusters
  sequences by similarity but does not call BLAST or any reference database.
  To identify species, take representative sequences from
  `cleaned_sequences.fasta` (or the whole file) and BLAST them at
  https://blast.ncbi.nlm.nih.gov/Blast.cgi against the "16S ribosomal RNA
  sequences (Bacteria and Archaea)" database. This can also be automated with
  Biopython's `Bio.Blast.NCBIWWW.qblast()`.
- **Thresholds are adjustable**: QC length/quality cutoffs, overlap
  identity/length cutoffs, and the clustering identity threshold are all
  named constants near the top of the script.
- **Primer sequences**: assume standard
  universal bacterial 16S primers; change these if you used different primers.
