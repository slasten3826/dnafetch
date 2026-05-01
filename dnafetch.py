#!/usr/bin/env python3
# dnafetch.py — DNA ProcessLang Translator
# Translates DNA sequences into ProcessLang operator descriptions (.dpl)
# For machine-to-machine communication. Humans need not apply.
#
# Usage:
#   python dnafetch.py single gene.fasta
#   python dnafetch.py batch gencode.v47.pc_transcripts.fa.gz -o output/
#   python dnafetch.py compare output/BRCA1.dpl output/TP53.dpl

import sys
import os
import time
import argparse

from fetch import fetch_single, fetch_batch
from translate import translate_gene
from describe import generate_dpl, write_dpl

LOGO = r"""
 ▄▄                                                                             ▄▄
▐██▌  ██████╗ ███╗   ██╗ █████╗ ███████╗███████╗████████╗ ██████╗ ██╗  ██╗  ▐██▌
 ▐█▌  ██╔══██╗████╗  ██║██╔══██╗██╔════╝██╔════╝╚══██╔══╝██╔════╝ ██║  ██║   ▐█▌
 ▐█▌  ██║  ██║██╔██╗ ██║███████║█████╗  █████╗     ██║   ██║      ███████║   ▐█▌
 ▐█▌  ██║  ██║██║╚██╗██║██╔══██║██╔══╝  ██╔══╝     ██║   ██║      ██╔══██║   ▐█▌
▐██▌  ██████╔╝██║ ╚████║██║  ██║██║     ███████╗   ██║   ╚██████╗ ██║  ██║  ▐██▌
 ▀▀   ╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝     ╚══════╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝   ▀▀
 ▐█▌                                                                           ▐█▌
▐██▌          >>  MAPPING THE SOURCE CODE OF LIFE  <<                         ▐██▌
 ▐█▌                                                                           ▐█▌
  ▀▀ ▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄ ▀▀
"""


def cmd_single(args):
    """Process a single FASTA file."""
    print(f"[FETCH] {args.input}")
    gene = fetch_single(args.input)
    print(f"  Gene: {gene.name} | {gene.organism}")
    print(f"  CDS: {gene.cds_length} bp ({gene.codon_count} codons)")

    print(f"[TRANSLATE] Mapping to ProcessLang...")
    result = translate_gene(gene)

    print(f"[DESCRIBE] Generating .dpl...")
    dpl_content = generate_dpl(result)

    # Output
    if args.output:
        outpath = args.output
    else:
        outpath = f"{gene.name}.dpl"

    if args.stdout:
        print(dpl_content)
    else:
        write_dpl(result, outpath)
        print(f"[DONE] {outpath}")


def cmd_batch(args):
    """Process a batch FASTA file (e.g., GENCODE .gz)."""
    outdir = args.output or "output/profiles"
    os.makedirs(outdir, exist_ok=True)

    print(f"[BATCH] Processing {args.input}")
    print(f"[OUTPUT] {outdir}/")

    limit = args.limit or 0
    start = time.time()
    count = 0
    errors = 0

    for gene in fetch_batch(args.input, limit=limit):
        try:
            result = translate_gene(gene)
            outpath = os.path.join(outdir, f"{gene.name}.dpl")
            write_dpl(result, outpath)
            count += 1

            if count % 500 == 0:
                elapsed = time.time() - start
                rate = count / elapsed
                print(f"  [{count}] {gene.name} | {rate:.0f} genes/sec")

        except Exception as e:
            errors += 1
            if args.verbose:
                print(f"  [ERROR] {gene.name}: {e}", file=sys.stderr)

    elapsed = time.time() - start
    print(f"\n[DONE] {count} genes processed in {elapsed:.1f}s")
    if errors:
        print(f"  {errors} errors skipped")
    print(f"  Output: {outdir}/")


def cmd_compare(args):
    """Compare two .dpl files."""
    print(f"[COMPARE] {args.file1} vs {args.file2}")

    with open(args.file1) as f:
        dpl1 = f.read()
    with open(args.file2) as f:
        dpl2 = f.read()

    # Parse key metrics from both
    def parse_dpl_values(content):
        values = {}
        for line in content.split("\n"):
            line = line.strip()
            if "=" in line and not line.startswith("#") and not line.startswith("["):
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.strip()
                # strip inline comments
                if "#" in val:
                    val = val[:val.index("#")].strip()
                try:
                    values[key] = float(val)
                except ValueError:
                    values[key] = val
        return values

    v1 = parse_dpl_values(dpl1)
    v2 = parse_dpl_values(dpl2)

    name1 = v1.get("gene", "file1")
    name2 = v2.get("gene", "file2")

    print(f"\n{'METRIC':<35} {'<'+str(name1)+'>':<20} {'<'+str(name2)+'>':<20} DELTA")
    print("─" * 95)

    compare_keys = [
        "dominant_pct", "FLOW", "OBSERVE", "LOGIC", "CHOOSE",
        "gc_content", "wobble_gc", "wobble_entropy",
        "complementary_balance", "top2_gap",
        "hydrophobic_ratio", "aa_std",
        "codon_diversity", "emergence_potential"
    ]

    for key in compare_keys:
        val1 = v1.get(key, "—")
        val2 = v2.get(key, "—")

        if isinstance(val1, float) and isinstance(val2, float):
            delta = val2 - val1
            sign = "+" if delta > 0 else ""
            print(f"  {key:<33} {val1:<20.4f} {val2:<20.4f} {sign}{delta:.4f}")
        else:
            print(f"  {key:<33} {str(val1):<20} {str(val2):<20}")

    # Compare string classifications
    str_keys = [
        "dominant_operator", "flow_state", "encoding_efficiency",
        "observation_level", "cycle_stability", "dissolution_potential",
        "connection_type", "manifestation_class", "fingerprint"
    ]

    print(f"\n{'CLASSIFICATION':<35} {'<'+str(name1)+'>':<25} {'<'+str(name2)+'>':<25}")
    print("─" * 85)
    for key in str_keys:
        val1 = v1.get(key, "—")
        val2 = v2.get(key, "—")
        match = "=" if val1 == val2 else "≠"
        print(f"  {key:<33} {str(val1):<25} {str(val2):<25} {match}")


def main():
    print(LOGO)
    parser = argparse.ArgumentParser(
        prog="dnafetch",
        description="DNAfetch v3 — DNA ProcessLang Translator"
    )

    subparsers = parser.add_subparsers(dest="command")

    # single
    p_single = subparsers.add_parser("single", help="Process a single FASTA file")
    p_single.add_argument("input", help="FASTA file path")
    p_single.add_argument("-o", "--output", help="Output .dpl file path")
    p_single.add_argument("--stdout", action="store_true", help="Print to stdout")

    # batch
    p_batch = subparsers.add_parser("batch", help="Process batch FASTA (supports .gz)")
    p_batch.add_argument("input", help="FASTA file path (.fa, .fa.gz)")
    p_batch.add_argument("-o", "--output", help="Output directory")
    p_batch.add_argument("-l", "--limit", type=int, help="Limit number of genes")
    p_batch.add_argument("-v", "--verbose", action="store_true")

    # compare
    p_compare = subparsers.add_parser("compare", help="Compare two .dpl files")
    p_compare.add_argument("file1", help="First .dpl file")
    p_compare.add_argument("file2", help="Second .dpl file")

    args = parser.parse_args()

    if args.command == "single":
        cmd_single(args)
    elif args.command == "batch":
        cmd_batch(args)
    elif args.command == "compare":
        cmd_compare(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
