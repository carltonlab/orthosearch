#!/usr/bin/env python3
import argparse, csv, os

def f(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")

def isnan(x: float) -> bool:
    return x != x

def nonempty(p: str) -> bool:
    return os.path.exists(p) and os.path.getsize(p) > 0

def fails_len_cov(row, min_qcov, min_tcov, min_lr, max_lr) -> bool:
    qcov = f(row.get("qcov", "nan"))
    tcov = f(row.get("tcov", "nan"))
    lr = f(row.get("len_ratio", "nan"))
    if isnan(qcov) or qcov < min_qcov:
        return True
    if isnan(tcov) or tcov < min_tcov:
        return True
    if isnan(lr) or lr < min_lr or lr > max_lr:
        return True
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--orthologs-dir", required=True)
    ap.add_argument("--suffix", default="protein.faa")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-qcov", type=float, required=True)
    ap.add_argument("--min-tcov", type=float, required=True)
    ap.add_argument("--min-len-ratio", type=float, required=True)
    ap.add_argument("--max-len-ratio", type=float, required=True)
    ap.add_argument("--require-at-least-one", action="store_true")
    args = ap.parse_args()

    rejects = []
    with open(args.manifest, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if not row.get("species"):
                continue
            keep = row.get("keep", "")
            if keep != "1":
                rejects.append(row["species"])
                continue
            # backward compatibility: also treat failures in case keep not set
            if fails_len_cov(row, args.min_qcov, args.min_tcov, args.min_len_ratio, args.max_len_ratio):
                rejects.append(row["species"])

    paths = [os.path.join(args.orthologs_dir, f"{sp}.{args.suffix}") for sp in rejects]
    good = [p for p in paths if nonempty(p)]
    if args.require_at_least_one and not good:
        raise SystemExit(f"[collect_rejected_by_manifest] No rejected sequences found for {args.out}")

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
