#!/usr/bin/env python3
import argparse, os

def nonempty(p: str) -> bool:
    return os.path.exists(p) and os.path.getsize(p) > 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-species", required=True)
    ap.add_argument("--orthologs-dir", required=True)
    ap.add_argument("--suffix", required=True)  # protein.faa or cds.fna
    ap.add_argument("--out", required=True)
    ap.add_argument("--require-at-least-one", action="store_true")
    args = ap.parse_args()

    with open(args.keep_species) as fh:
        keep = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]

    inputs = [os.path.join(args.orthologs_dir, f"{sp}.{args.suffix}") for sp in keep]
    good = [p for p in inputs if nonempty(p)]
    if args.require_at_least_one and not good:
        raise SystemExit(f"[filter_combined] No non-empty inputs for {args.out}")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as out:
        for p in good:
            with open(p) as fh:
                s = fh.read()
            if s and not s.endswith("\n"):
                s += "\n"
            out.write(s)

if __name__ == "__main__":
    main()
