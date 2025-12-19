#!/usr/bin/env python3
import argparse, csv, os

def f(x):
    try: return float(x)
    except Exception: return float("nan")

def isnan(x): return x != x

def nonempty(p: str) -> bool:
    return os.path.exists(p) and os.path.getsize(p) > 0

def contains_unknown_x(path: str) -> bool:
    if not nonempty(path):
        return False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            if "X" in line.upper():
                return True
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--topk-tsvs", nargs="+", required=True)
    ap.add_argument("--qc-tsvs", nargs="+", required=True)
    ap.add_argument("--out-manifest", required=True)
    ap.add_argument("--out-keep", required=True)
    ap.add_argument("--min-qcov", type=float, default=0.6)
    ap.add_argument("--min-tcov", type=float, default=0.2)
    ap.add_argument("--min-len-ratio", type=float, default=0.6)
    ap.add_argument("--max-len-ratio", type=float, default=1.8)
    ap.add_argument("--min-bits", type=float, default=0.0)
    ap.add_argument("--orthologs-dir", default="orthologs")
    ap.add_argument("--suffix", default="protein.faa")
    args = ap.parse_args()

    qc_flags = {}
    for qc_path in args.qc_tsvs:
        if not os.path.exists(qc_path) or os.path.getsize(qc_path) == 0:
            continue
        with open(qc_path, newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                qc_flags[row.get("species","")] = row.get("flags","")

    rows = []
    for tsv in args.topk_tsvs:
        if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
            continue
        with open(tsv, newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                if row.get("rank") != "1":
                    continue
                sp = row.get("species","")
                target = row.get("target","")
                qcov = f(row.get("qcov","nan"))
                tcov = f(row.get("tcov","nan"))
                bits = f(row.get("bits","nan"))
                lr = f(row.get("len_ratio","nan"))
                keep = True
                if (not isnan(qcov)) and qcov < args.min_qcov: keep = False
                if (not isnan(tcov)) and tcov < args.min_tcov: keep = False
                if (not isnan(bits)) and bits < args.min_bits: keep = False
                if (not isnan(lr)) and (lr < args.min_len_ratio or lr > args.max_len_ratio): keep = False

                flags = ",".join([x for x in [row.get("flags",""), qc_flags.get(sp,"")] if x])
                prot_path = os.path.join(args.orthologs_dir, f"{sp}.{args.suffix}")
                has_x = contains_unknown_x(prot_path)
                if has_x:
                    keep = False
                    flags = ",".join([x for x in [flags, "contains_X"] if x])
                rows.append({
                    "species": sp,
                    "query": row.get("query",""),
                    "target": target,
                    "bits": bits,
                    "evalue": row.get("evalue",""),
                    "qcov": qcov,
                    "tcov": tcov,
                    "fident": f(row.get("fident","nan")),
                    "qlen": row.get("qlen",""),
                    "tlen": row.get("tlen",""),
                    "len_ratio": lr,
                    "flags": flags,
                    "keep": "1" if keep else "0",
                })
                break

    os.makedirs(os.path.dirname(args.out_manifest), exist_ok=True)
    with open(args.out_manifest, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["species","query","target","bits","evalue","qcov","tcov","fident","qlen","tlen","len_ratio","flags","keep"], delimiter="\t")
        w.writeheader()
        for row in sorted(rows, key=lambda x: x["species"]):
            w.writerow(row)

    kept = sorted({r["species"] for r in rows if r["keep"] == "1" and r["target"]})
    os.makedirs(os.path.dirname(args.out_keep), exist_ok=True)
    with open(args.out_keep, "w") as out:
        for sp in kept:
            out.write(sp + "\n")

if __name__ == "__main__":
    main()
