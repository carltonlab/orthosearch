#!/usr/bin/env python3
import argparse, os
from Bio import SeqIO

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--id-file", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    chosen = ""
    if os.path.exists(args.id_file):
        with open(args.id_file) as fh:
            chosen = fh.readline().strip()

    if not chosen:
        open(args.out, "w").close()
        return

    found = False
    with open(args.out, "w") as out_f:
        for rec in SeqIO.parse(args.fasta, "fasta"):
            if rec.id == chosen:
                SeqIO.write(rec, out_f, "fasta")
                found = True
                break

    if not found:
        msg = f"[extract_fasta_record] ID '{chosen}' not found in {args.fasta}"
        if args.strict:
            raise SystemExit(msg)
        open(args.out, "w").close()

if __name__ == "__main__":
    main()
