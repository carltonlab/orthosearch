#!/usr/bin/env python3
import argparse, csv, os, sys
from dataclasses import dataclass
from typing import List

@dataclass
class Hit:
    query: str
    target: str
    fident: float
    alnlen: int
    evalue: float
    bits: float
    qcov: float
    tcov: float
    qlen: int
    tlen: int

def f(x: str) -> float:
    try: return float(x)
    except Exception: return float("nan")

def i(x: str) -> int:
    try: return int(float(x))
    except Exception: return 0

def isnan(x: float) -> bool:
    return x != x

def load_hits(path: str) -> List[Hit]:
    if (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        req = {"query","target","fident","alnlen","evalue","bits","qcov","tcov","qlen","tlen"}
        miss = req - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"[pick_best_hit] Missing columns in {path}: {sorted(miss)}")
        hits = []
        for row in r:
            if not row["target"]:
                    continue
            hits.append(Hit(
                query=row["query"],
                target=row["target"],
                fident=f(row["fident"]),
                alnlen=i(row["alnlen"]),
                evalue=f(row["evalue"]),
                bits=f(row["bits"]),
                qcov=f(row["qcov"]),
                tcov=f(row["tcov"]),
                qlen=i(row["qlen"]),
                tlen=i(row["tlen"]),
            ))
        return hits

def score(h: Hit) -> float:
    # Primary: bits; Secondary: min coverage; Tertiary: identity
    if isnan(h.bits): return -1e18
    cov = 0.0
    if (not isnan(h.qcov)) and (not isnan(h.tcov)):
        cov = min(h.qcov, h.tcov)
    ident = 0.0 if isnan(h.fident) else h.fident
    return h.bits + 10.0*cov + 0.1*ident

def read_override(overrides_path: str, species: str):
    if not overrides_path: return None
    if not os.path.exists(overrides_path): return None
    with open(overrides_path) as fh:
        for line in fh:
            if (not line.strip()) or line.startswith("#"): 
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] == species:
                return parts[1]
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--min-qcov", type=float, default=0.6)
    ap.add_argument("--min-tcov", type=float, default=0.6)
    ap.add_argument("--min-fident", type=float, default=0.0)
    ap.add_argument("--min-bits", type=float, default=0.0)
    ap.add_argument("--min-len-ratio", type=float, default=0.6)
    ap.add_argument("--max-len-ratio", type=float, default=1.5)
    ap.add_argument("--overrides", default="")
    ap.add_argument("--out-best-id", required=True)
    ap.add_argument("--out-topk", required=True)
    ap.add_argument("--out-qc", required=True)
    args = ap.parse_args()

    override = read_override(args.overrides, args.species)
    hits = load_hits(args.hits)
    hits_sorted = sorted(hits, key=score, reverse=True)

    os.makedirs(os.path.dirname(args.out_topk), exist_ok=True)
    with open(args.out_topk, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species","rank","query","target","bits","evalue","qcov","tcov","fident","qlen","tlen","len_ratio","flags"])
        for rank, h in enumerate(hits_sorted[:args.top_k], start=1):
            lr = (h.tlen / h.qlen) if h.qlen else float("nan")
            flags = []
            if (not isnan(lr)) and (lr < args.min_len_ratio or lr > args.max_len_ratio):
                flags.append("len_ratio_outside")
            if (not isnan(h.qcov)) and h.qcov < args.min_qcov:
                flags.append("low_qcov")
            if (not isnan(h.tcov)) and h.tcov < args.min_tcov:
                flags.append("low_tcov")
            if (not isnan(h.fident)) and h.fident < args.min_fident:
                flags.append("low_fident")
            if (not isnan(h.bits)) and h.bits < args.min_bits:
                flags.append("low_bits")
            w.writerow([args.species, rank, h.query, h.target, h.bits, h.evalue, h.qcov, h.tcov, h.fident, h.qlen, h.tlen, lr, ",".join(flags)])

    chosen = None
    reason = ""
    if override:
        chosen = override
        reason = "override"
    else:
        for h in hits_sorted:
            if isnan(h.bits) or h.bits < args.min_bits: 
                continue
            if (not isnan(h.qcov)) and h.qcov < args.min_qcov:
                continue
            if (not isnan(h.tcov)) and h.tcov < args.min_tcov:
                continue
            if (not isnan(h.fident)) and h.fident < args.min_fident:
                continue
            chosen = h.target
            reason = "best_passing"
            break

    flags = []
    if len(hits_sorted) >= 2 and (not isnan(hits_sorted[0].bits)) and (not isnan(hits_sorted[1].bits)):
        if hits_sorted[0].bits > 0 and (hits_sorted[1].bits / hits_sorted[0].bits) >= 0.97:
            flags.append("top2_close_duplication_candidate")

    status = "chosen" if chosen else "no_passing_hit"

    os.makedirs(os.path.dirname(args.out_best_id), exist_ok=True)
    with open(args.out_best_id, "w") as fh:
        fh.write("" if not chosen else chosen)
        fh.write("\n")

    os.makedirs(os.path.dirname(args.out_qc), exist_ok=True)
    with open(args.out_qc, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species","status","chosen_id","reason","n_hits","flags"])
        w.writerow([args.species, status, "" if not chosen else chosen, reason, len(hits), ",".join(flags)])

if __name__ == "__main__":
    main()
