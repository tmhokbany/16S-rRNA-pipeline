"""
Bacterial 16S rRNA gene sequencing pipeline

Requirements (pip install):
  biopython
  pyfamsa
  matplotlib

"""

import argparse
import glob
import itertools
import json
import os
import re

from Bio import SeqIO, Align, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.BaseTree import BranchColor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyfamsa import Aligner, Sequence

# Configuration constants

PRIMER_27F = "AGAGTTTGATCMTGGCTCAG"
PRIMER_1492R = "TACGGYTACCTTGTTACGACTT"

# IUPAC ambiguity code -> regex character class, used to detect primers
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "M": "[AC]", "R": "[AG]", "W": "[AT]", "S": "[CG]",
    "Y": "[CT]", "K": "[GT]", "V": "[ACG]", "H": "[ACT]",
    "D": "[AGT]", "B": "[CGT]", "N": "[ACGT]",
}

QC_MIN_LENGTH = 350      
QC_MIN_MEAN_Q = 20        
OVERLAP_MIN_IDENTITY = 0.97   
OVERLAP_MIN_LENGTH = 50      
CLUSTER_IDENTITY_THRESHOLD = 97.0  

def iupac_to_regex(primer: str) -> str:
    return "".join(IUPAC.get(base, base) for base in primer)

# 1 quality trimming
def mott_trim(quals, q_threshold=20):
    vals = [q - q_threshold for q in quals]
    max_sum = 0
    cur_sum = 0
    cur_start = 0
    best_start, best_end = 0, 0
    for i, v in enumerate(vals):
        if cur_sum <= 0:
            cur_start = i
            cur_sum = v
        else:
            cur_sum += v
        if cur_sum > max_sum:
            max_sum = cur_sum
            best_start = cur_start
            best_end = i + 1
    return best_start, best_end


def process_ab1_file(path: str) -> dict:
    rec = SeqIO.read(path, "abi")
    seq = str(rec.seq)
    qual = rec.letter_annotations["phred_quality"]

    start, end = mott_trim(qual, q_threshold=QC_MIN_MEAN_Q)
    trimmed_seq = seq[start:end]
    trimmed_qual = qual[start:end]
    mean_q = sum(trimmed_qual) / len(trimmed_qual) if trimmed_qual else 0

    name = os.path.basename(path).replace(".ab1", "")
    sample, read_label = name.split("_", 1)
    direction = "F" if "27-F" in read_label else "R"

    status = "PASS" if (len(trimmed_seq) >= QC_MIN_LENGTH and mean_q >= QC_MIN_MEAN_Q) else "FAIL"

    # Screen the first 30bp of the trimmed read for residual primer sequence.
    primer_trimmed = trimmed_seq
    primer_found = None
    search_region = trimmed_seq[:30]
    primer = PRIMER_27F if direction == "F" else PRIMER_1492R
    match = re.search(iupac_to_regex(primer), search_region)
    if match:
        primer_trimmed = trimmed_seq[match.end():]
        primer_found = primer

    return {
        "file": os.path.basename(path),
        "sample": sample,
        "direction": direction,
        "raw_length": len(seq),
        "raw_mean_q": round(sum(qual) / len(qual), 1),
        "trimmed_length": len(trimmed_seq),
        "trimmed_mean_q": round(mean_q, 1),
        "primer_found": primer_found,
        "final_length": len(primer_trimmed),
        "status": status,
        "sequence": primer_trimmed,
    }

def run_trim_qc(input_dir: str, output_dir: str) -> list:
    files = sorted(glob.glob(os.path.join(input_dir, "*.ab1")))
    if not files:
        raise FileNotFoundError(f"No .ab1 files found in {input_dir}")

    results = [process_ab1_file(f) for f in files]
    results.sort(key=lambda r: (int(r["sample"]), r["direction"]))

    out_path = os.path.join(output_dir, "01_trim_qc_results.json")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)

    n_pass = sum(1 for r in results if r["status"] == "PASS")
    print(f"[trim/QC] {len(results)} reads processed, {n_pass} passed QC -> {out_path}")
    return results

# 2: per-sample consensus assembly (merge F + R reads)
def build_aligner():
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -2
    return aligner

def try_merge(aligner, f_seq: str, r_seq_rc: str):
    """
    Attempt to merge a forward read with a reverse-complemented reverse
    read by finding their overlap via local alignment. Returns the merged
    sequence if a confident overlap (>=97% identity, >=50bp) is found,
    otherwise None.
    """
    alignment = aligner.align(f_seq, r_seq_rc)[0]
    f_blocks = alignment.aligned[0]
    r_blocks = alignment.aligned[1]
    if len(f_blocks) == 0:
        return None

    f_start, f_end = f_blocks[0][0], f_blocks[-1][1]
    r_end = r_blocks[-1][1]

    s1, s2 = str(alignment[0]), str(alignment[1])
    matches = sum(1 for a, b in zip(s1, s2) if a == b and a != "-")
    aligned_cols = sum(1 for a, b in zip(s1, s2) if a != "-" and b != "-")
    if aligned_cols == 0:
        return None

    identity = matches / aligned_cols
    if identity < OVERLAP_MIN_IDENTITY or aligned_cols < OVERLAP_MIN_LENGTH:
        return None

    merged_seq = f_seq[:f_start] + f_seq[f_start:f_end] + r_seq_rc[r_end:]
    return {
        "identity": round(identity, 4),
        "overlap_len": aligned_cols,
        "merged_seq": merged_seq,
        "merged_len": len(merged_seq),
    }

def run_consensus_assembly(qc_results: list, output_dir: str) -> list:
    by_sample = {}
    for r in qc_results:
        by_sample.setdefault(r["sample"], {})[r["direction"]] = r

    aligner = build_aligner()
    consensus_records = []

    for sample in sorted(by_sample.keys(), key=int):
        entry = by_sample[sample]
        f = entry.get("F")
        r = entry.get("R")
        f_pass = bool(f) and f["status"] == "PASS"
        r_pass = bool(r) and r["status"] == "PASS"

        if f_pass and r_pass:
            f_seq = f["sequence"]
            r_seq_rc = str(Seq(r["sequence"]).reverse_complement())
            merge = try_merge(aligner, f_seq, r_seq_rc)
            if merge:
                consensus_records.append({
                    "sample": sample,
                    "method": "F+R merged consensus (overlap assembly)",
                    "sequence": merge["merged_seq"],
                    "length": merge["merged_len"],
                    "overlap_identity": merge["identity"],
                    "overlap_len": merge["overlap_len"],
                })
            else:
                best = f if f["final_length"] >= r["final_length"] else r
                seq = best["sequence"] if best is f else str(Seq(r["sequence"]).reverse_complement())
                consensus_records.append({
                    "sample": sample,
                    "method": f"No reliable F/R overlap found; used {best['direction']} read alone",
                    "sequence": seq,
                    "length": len(seq),
                })
        elif f_pass or r_pass:
            best = f if f_pass else r
            seq = best["sequence"] if best is f else str(Seq(r["sequence"]).reverse_complement())
            consensus_records.append({
                "sample": sample,
                "method": f"Only {best['direction']} read passed QC; used alone",
                "sequence": seq,
                "length": len(seq),
            })
        else:
            print(f"[assembly] Sample {sample}: EXCLUDED - both reads failed QC")

    out_path = os.path.join(output_dir, "02_consensus.json")
    with open(out_path, "w") as fh:
        json.dump(consensus_records, fh, indent=2)

    print(f"[assembly] {len(consensus_records)}/{len(by_sample)} samples retained -> {out_path}")
    return consensus_records

# 3 cleaned up FASTA
def write_fasta(records: list, path: str, seq_key: str = "sequence",
                 header_fn=None, wrap: int = 70):
    with open(path, "w") as fh:
        for rec in records:
            header = header_fn(rec) if header_fn else rec.get("id", "seq")
            fh.write(f">{header}\n")
            seq = rec[seq_key]
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")


def run_write_cleaned_fasta(consensus_records: list, output_dir: str) -> str:
    out_path = os.path.join(output_dir, "cleaned_sequences.fasta")
    write_fasta(
        consensus_records,
        out_path,
        header_fn=lambda r: f"Sample_{r['sample']} length={r['length']} method={r['method']}",
    )
    print(f"[fasta] wrote {len(consensus_records)} sequences -> {out_path}")
    return out_path

# 4 multiple sequence alignment (FAMSA)

def run_msa(consensus_records: list, output_dir: str) -> str:
    seqs = [
        Sequence(f"Sample_{r['sample']}".encode(), r["sequence"].encode())
        for r in consensus_records
    ]
    aligner = Aligner(guide_tree="upgma")
    msa = aligner.align(seqs)

    out_path = os.path.join(output_dir, "aligned.fasta")
    with open(out_path, "w") as fh:
        for s in msa:
            fh.write(f">{s.id.decode()}\n")
            seq = s.sequence.decode()
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    print(f"[msa] aligned {len(msa)} sequences, length {len(msa[0].sequence)} -> {out_path}")
    return out_path

# 5 pairwise identity + clustering + GC content
def pct_identity(a: str, b: str) -> float:
    matches = valid = 0
    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue
        valid += 1
        if x == y:
            matches += 1
    return 100 * matches / valid if valid else 0.0


def cluster_by_identity(names: list, identity_matrix: dict, threshold: float) -> dict:
    parent = {name: name for name in names}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for (a, b), v in identity_matrix.items():
        if v >= threshold:
            union(a, b)

    clusters = {}
    for name in names:
        clusters.setdefault(find(name), []).append(name)
    return clusters

def run_identity_clustering(aligned_path: str, consensus_records: list, output_dir: str) -> dict:
    aln = AlignIO.read(aligned_path, "fasta")
    records = {r.id: str(r.seq) for r in aln}
    names = sorted(records.keys(), key=lambda x: int(x.split("_")[1]))

    identity_matrix = {
        (a, b): pct_identity(records[a], records[b])
        for a, b in itertools.combinations(names, 2)
    }

    clusters = cluster_by_identity(names, identity_matrix, CLUSTER_IDENTITY_THRESHOLD)

    # GC content per sample, as an independent cross-check
    gc_by_sample = {
        f"Sample_{r['sample']}": round(gc_fraction(r["sequence"]) * 100, 1)
        for r in consensus_records
    }

    cluster_summary = []
    for i, (root, members) in enumerate(sorted(clusters.items(), key=lambda x: -len(x[1])), 1):
        members_sorted = sorted(members, key=lambda x: int(x.split("_")[1]))
        gcs = [gc_by_sample[m] for m in members_sorted]
        cluster_summary.append({
            "cluster": f"C{i}",
            "members": members_sorted,
            "n": len(members_sorted),
            "gc_percent": gcs,
            "mean_gc": round(sum(gcs) / len(gcs), 1),
        })
        print(f"[cluster] C{i} (n={len(members_sorted)}): {members_sorted}  mean GC={sum(gcs)/len(gcs):.1f}%")

    out_path = os.path.join(output_dir, "identity_clusters.json")
    with open(out_path, "w") as fh:
        json.dump(cluster_summary, fh, indent=2)
    print(f"[cluster] wrote cluster summary -> {out_path}")

    # sample -> cluster label lookup, used later for tree coloring
    sample_to_cluster = {}
    for c in cluster_summary:
        for m in c["members"]:
            sample_to_cluster[m] = c["cluster"]
    return sample_to_cluster

# 6 phylogenetic tree
def run_tree(aligned_path: str, output_dir: str) -> str:
    aln = AlignIO.read(aligned_path, "fasta")
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(aln)

    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)
    nj_tree.root_at_midpoint()
    nj_tree.ladderize()

    out_path = os.path.join(output_dir, "phylo_tree_nj.nwk")
    Phylo.write(nj_tree, out_path, "newick")
    print(f"[tree] Neighbor-Joining tree written -> {out_path}")
    return out_path


# fixed palette
CLUSTER_COLORS = [
    (178, 34, 34), (255, 140, 0), (46, 139, 87), (65, 105, 225),
    (128, 0, 128), (0, 128, 128), (105, 105, 105), (0, 0, 0),
    (218, 112, 214), (139, 69, 19),
]

def species_from_hitdef(hit_def: str) -> str:
    """Pull just 'Genus species' out of a full BLAST hit definition string."""
    parts = hit_def.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else hit_def


def load_cluster_species(output_dir: str) -> dict:
    """Map cluster label (e.g. 'C1') -> short BLAST species name, read from
    blast_results.json if blast_identify.py has already been run on this
    output folder. Returns {} if it hasn't (tree falls back to cluster
    labels only). This is one BLAST call per cluster REPRESENTATIVE, so it
    assumes every member of a cluster is the same species."""
    path = os.path.join(output_dir, "blast_results.json")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        blast_results = json.load(fh)
    species = {}
    for cluster_id, r in blast_results.items():
        hits = r.get("top_hits") or []
        if hits:
            species[cluster_id] = species_from_hitdef(hits[0]["hit_def"])
    return species


def load_sample_species(output_dir: str) -> dict:
    """Map sample label (e.g. 'Sample_21') -> confirmed species name, read
    from confirmed_species.json if present. This is per-SAMPLE (not just
    per-cluster-representative), so it can distinguish species within a
    cluster that a single representative BLAST call would miss."""
    path = os.path.join(output_dir, "confirmed_species.json")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def run_tree_plot(tree_path: str, sample_to_cluster: dict, output_dir: str) -> str:
    tree = Phylo.read(tree_path, "newick")

    distinct_clusters = sorted(set(sample_to_cluster.values()))
    color_for_cluster = {
        c: CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
        for i, c in enumerate(distinct_clusters)
    }

    cluster_species = load_cluster_species(output_dir)
    sample_species = load_sample_species(output_dir)

    for clade in tree.get_terminals():
        cluster = sample_to_cluster.get(clade.name, "?")
        rgb = color_for_cluster.get(cluster, (0, 0, 0))
        clade.color = BranchColor(*rgb)
        label = sample_species.get(clade.name) or cluster_species.get(cluster, cluster)
        clade.name = f"{clade.name}  [{label}]"

    def label_only_terminals(clade):
        return clade.name if clade.is_terminal() else None

    n_leaves = tree.count_terminals()
    fig_height = max(9, n_leaves * 0.35)
    fig = plt.figure(figsize=(14, fig_height))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        branch_labels=None,
        label_func=label_only_terminals,
    )
    label_note = "species (BLAST-confirmed)" if (sample_species or cluster_species) else "cluster (BLAST not yet run)"
    ax.set_title("Phylogenetic tree of 16S rRNA gene consensus sequences\n"
                 f"(leaf color = sequence-identity cluster, leaf label = {label_note})", fontsize=11)

    xmin, xmax = ax.get_xlim()
    ax.set_xlim(xmin, xmax + 0.25 * (xmax - xmin))
    plt.tight_layout()

    out_path = os.path.join(output_dir, "phylo_tree.png")
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[tree] rendered image -> {out_path}")
    return out_path

# Main

def main():
    parser = argparse.ArgumentParser(description="Bacterial 16S rRNA sequencing pipeline")
    parser.add_argument("--input", required=True, help="Folder containing .ab1 chromatogram files")
    parser.add_argument("--output", required=True, help="Folder to write results into")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    qc_results = run_trim_qc(args.input, args.output)
    consensus_records = run_consensus_assembly(qc_results, args.output)
    run_write_cleaned_fasta(consensus_records, args.output)
    aligned_path = run_msa(consensus_records, args.output)
    sample_to_cluster = run_identity_clustering(aligned_path, consensus_records, args.output)
    tree_path = run_tree(aligned_path, args.output)
    run_tree_plot(tree_path, sample_to_cluster, args.output)

    print("\nDone. All outputs written to:", args.output)


if __name__ == "__main__":
    main()
