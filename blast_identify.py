"""
BLAST identification of cluster 16s sequences using biopython's NCBIWWW.qblast() function.

Requirements (pip install):
  biopython
"""

# NCBI database identifier for the curated 16S ribosomal RNA reference set (Bacteria and Archaea)
import argparse
import json
import os
import time
 
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

BLAST_DATABASE = "rRNA_typestrains/16S_ribosomal_RNA"
BLAST_PROGRAM = "blastn"
TOP_N_HITS = 3            
SECONDS_BETWEEN_REQUESTS = 5  


def load_cluster_representatives(results_dir: str) -> dict:
    """Pick the first member of each cluster as its representative sequence,
    and pull that sequence out of cleaned_sequences.fasta."""
    clusters_path = os.path.join(results_dir, "identity_clusters.json")
    fasta_path = os.path.join(results_dir, "cleaned_sequences.fasta")

    with open(clusters_path) as fh:
        clusters = json.load(fh)

    # cleaned_sequences.fasta headers look like ">Sample_26 length=... method=..."
    # so we index by the first whitespace-delimited token of each record's id.
    seq_index = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}

    representatives = {}
    for c in clusters:
        rep_name = c["members"][0]           
        if rep_name not in seq_index:
            print(f"[warn] {rep_name} not found in {fasta_path}, skipping cluster {c['cluster']}")
            continue
        representatives[c["cluster"]] = {
            "representative": rep_name,
            "cluster_members": c["members"],
            "sequence": seq_index[rep_name],
        }
    return representatives


def blast_one_sequence(name: str, sequence: str):
    """Submit one sequence to NCBI BLAST and return the top N hits."""
    print(f"[blast] submitting {name} ({len(sequence)}bp) ... this can take a few minutes")
    result_handle = NCBIWWW.qblast(
        program=BLAST_PROGRAM,
        database=BLAST_DATABASE,
        sequence=sequence,
        megablast=True,
        hitlist_size=TOP_N_HITS,
    )

    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    hits = []
    for alignment in blast_record.alignments[:TOP_N_HITS]:
        hsp = alignment.hsps[0]   # best local alignment segment for this hit
        percent_identity = 100 * hsp.identities / hsp.align_length
        query_coverage = 100 * (hsp.query_end - hsp.query_start + 1) / len(sequence)
        hits.append({
            "hit_id": alignment.hit_id,
            "hit_def": alignment.hit_def,
            "percent_identity": round(percent_identity, 2),
            "query_coverage_percent": round(query_coverage, 1),
            "e_value": hsp.expect,
            "align_length": hsp.align_length,
        })

    print(f"[blast] {name}: top hit -> {hits[0]['hit_def'] if hits else 'NO HITS FOUND'}")
    return hits


def run_blast_all(representatives: dict) -> dict:
    results = {}
    cluster_ids = list(representatives.keys())
    for i, cluster_id in enumerate(cluster_ids):
        info = representatives[cluster_id]
        hits = blast_one_sequence(info["representative"], info["sequence"])
        results[cluster_id] = {
            "representative": info["representative"],
            "cluster_members": info["cluster_members"],
            "top_hits": hits,
        }
        # Be polite between submissions (skip the wait after the last one)
        if i < len(cluster_ids) - 1:
            time.sleep(SECONDS_BETWEEN_REQUESTS)
    return results


def write_outputs(results: dict, results_dir: str):
    json_path = os.path.join(results_dir, "blast_results.json")
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\n[output] full results -> {json_path}")

    txt_path = os.path.join(results_dir, "blast_summary.txt")
    with open(txt_path, "w") as fh:
        fh.write("BLAST identification summary\n")
        fh.write("=" * 60 + "\n\n")
        for cluster_id, r in results.items():
            fh.write(f"{cluster_id}  (representative: {r['representative']}, "
                      f"{len(r['cluster_members'])} members: {', '.join(r['cluster_members'])})\n")
            if not r["top_hits"]:
                fh.write("  NO HITS FOUND\n\n")
                continue
            for rank, hit in enumerate(r["top_hits"], 1):
                fh.write(f"  {rank}. {hit['hit_def']}\n")
                fh.write(f"     identity={hit['percent_identity']}%  "
                          f"coverage={hit['query_coverage_percent']}%  "
                          f"e-value={hit['e_value']}\n")
            fh.write("\n")
    print(f"[output] readable summary -> {txt_path}")


def main():
    parser = argparse.ArgumentParser(description="BLAST-identify 16S cluster representatives")
    parser.add_argument("--results", required=True,
                         help="Folder containing identity_clusters.json and cleaned_sequences.fasta "
                              "(the --output folder from bacterial_16s_pipeline.py)")
    args = parser.parse_args()

    representatives = load_cluster_representatives(args.results)
    print(f"[setup] {len(representatives)} cluster representatives to BLAST: "
          f"{[v['representative'] for v in representatives.values()]}")

    results = run_blast_all(representatives)
    write_outputs(results, args.results)

    print("\nDone.")


if __name__ == "__main__":
    main()