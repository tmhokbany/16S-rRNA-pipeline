"""
Confirm species identification by BLASTing a SECOND member of every
multi-sample cluster and checking it agrees with the first.
This is a double-check step: identity_clusters.json groups samples at 97%
sequence identity
"""

import argparse
import json
import os
import time

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

BLAST_DATABASE = "rRNA_typestrains/16S_ribosomal_RNA"
BLAST_PROGRAM = "blastn"
TOP_N_HITS = 3
RETRIES = 3
SECONDS_BETWEEN_RETRIES = 15
SECONDS_BETWEEN_REQUESTS = 5


def load_second_members(results_dir: str) -> dict:
    """For each cluster with >=2 members, pick a SECOND sample (not the one
    already BLASTed in blast_identify.py) to cross-check against."""
    with open(os.path.join(results_dir, "identity_clusters.json")) as fh:
        clusters = json.load(fh)
    seq_index = {rec.id: str(rec.seq)
                 for rec in SeqIO.parse(os.path.join(results_dir, "cleaned_sequences.fasta"), "fasta")}

    # what the first pass already called, so we can check agreement
    first_pass_path = os.path.join(results_dir, "blast_results.json")
    first_pass = {}
    if os.path.exists(first_pass_path):
        with open(first_pass_path) as fh:
            first_pass = json.load(fh)

    to_check = {}
    for c in clusters:
        if len(c["members"]) < 2:
            continue  # nothing to cross-check for a singleton
        second_member = c["members"][1]
        if second_member not in seq_index:
            continue
        first_hit_def = None
        if c["cluster"] in first_pass and first_pass[c["cluster"]]["top_hits"]:
            first_hit_def = first_pass[c["cluster"]]["top_hits"][0]["hit_def"]
        to_check[c["cluster"]] = {
            "sample": second_member,
            "sequence": seq_index[second_member],
            "first_pass_top_hit": first_hit_def,
        }
    return to_check


def blast_with_retries(name: str, sequence: str):
    last_error = None
    for attempt in range(1, RETRIES + 2):
        try:
            print(f"[confirm] submitting {name} ({len(sequence)}bp), attempt {attempt}...")
            handle = NCBIWWW.qblast(BLAST_PROGRAM, BLAST_DATABASE, sequence,
                                     megablast=True, hitlist_size=TOP_N_HITS)
            record = NCBIXML.read(handle)
            handle.close()
            hits = []
            for alignment in record.alignments[:TOP_N_HITS]:
                hsp = alignment.hsps[0]
                hits.append({
                    "hit_def": alignment.hit_def,
                    "percent_identity": round(100 * hsp.identities / hsp.align_length, 2),
                    "e_value": hsp.expect,
                })
            return hits
        except Exception as e:
            last_error = e
            print(f"[confirm] {name}: attempt {attempt} failed ({e!r})")
            if attempt <= RETRIES:
                time.sleep(SECONDS_BETWEEN_RETRIES)
    print(f"[confirm] {name}: gave up after {RETRIES + 1} attempts ({last_error!r})")
    return None


def species_from_hitdef(hit_def: str) -> str:
    """Pull just 'Genus species' out of a full hit definition string."""
    parts = hit_def.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else hit_def


def main():
    parser = argparse.ArgumentParser(description="Cross-check cluster identification with a 2nd BLAST per cluster")
    parser.add_argument("--results", required=True)
    args = parser.parse_args()

    to_check = load_second_members(args.results)
    print(f"[setup] {len(to_check)} clusters to cross-check: "
          f"{[v['sample'] for v in to_check.values()]}")

    out_json = os.path.join(args.results, "blast_confirmation.json")
    existing = {}
    if os.path.exists(out_json):
        with open(out_json) as fh:
            existing = json.load(fh)

    cluster_ids = list(to_check.keys())
    for i, cluster_id in enumerate(cluster_ids):
        if cluster_id in existing:
            print(f"[skip] {cluster_id} already confirmed, skipping")
            continue

        info = to_check[cluster_id]
        hits = blast_with_retries(info["sample"], info["sequence"])
        if hits is None:
            continue  # leave it out; re-running the script will retry it

        second_species = species_from_hitdef(hits[0]["hit_def"])
        first_species = species_from_hitdef(info["first_pass_top_hit"] or "")
        agrees = (second_species.lower() == first_species.lower())

        existing[cluster_id] = {
            "first_pass_sample_top_hit": info["first_pass_top_hit"],
            "second_sample": info["sample"],
            "second_sample_top_hits": hits,
            "agrees_with_first_pass": agrees,
        }
        with open(out_json, "w") as fh:
            json.dump(existing, fh, indent=2)

        status = "AGREES" if agrees else "*** DISAGREES - REVIEW ***"
        print(f"[confirm] {cluster_id}: {status}  "
              f"(1st: {first_species!r}  2nd: {second_species!r})")

        if i < len(cluster_ids) - 1:
            time.sleep(SECONDS_BETWEEN_REQUESTS)

    # Write readable summary
    txt_path = os.path.join(args.results, "blast_confirmation.txt")
    with open(txt_path, "w") as fh:
        fh.write("BLAST cross-check summary (2nd member per cluster)\n")
        fh.write("=" * 60 + "\n\n")
        for cluster_id, r in existing.items():
            status = "AGREES" if r["agrees_with_first_pass"] else "*** DISAGREES ***"
            fh.write(f"{cluster_id}: {status}\n")
            fh.write(f"  1st-pass sample top hit: {r['first_pass_sample_top_hit']}\n")
            fh.write(f"  2nd sample ({r['second_sample']}) top hit: "
                      f"{r['second_sample_top_hits'][0]['hit_def']}\n\n")
    print(f"\n[output] -> {out_json}\n[output] -> {txt_path}")

    missing = set(to_check.keys()) - set(existing.keys())
    if missing:
        print(f"\n[note] {len(missing)} cluster(s) still need cross-checking: {sorted(missing)}. "
              f"Re-run the same command to retry just those.")
    else:
        disagreements = [c for c, r in existing.items() if not r["agrees_with_first_pass"]]
        if disagreements:
            print(f"\n[REVIEW NEEDED] Disagreement found in: {disagreements}")
        else:
            print("\nAll cross-checked clusters agree with the first-pass identification.")


if __name__ == "__main__":
    main()
