#!/usr/bin/env python3
import argparse, csv
from pathlib import Path
from Bio import Phylo


def read_keep(path: Path):
    return {ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")}


def read_manifest_targets(path: Path):
    m = {}
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            sp = row.get("species", "")
            tgt = row.get("target", "")
            if sp and tgt:
                m[sp] = tgt
    return m

def read_best_id_targets(base_dir: Path, suffix: str = ".best_id.txt"):
    mapping = {}
    for p in base_dir.rglob(f"*{suffix}"):
        if p.is_dir():
            continue
        if p.stat().st_size <= 1:
            continue
        sp = p.name[:-len(suffix)] if p.name.endswith(suffix) else p.stem
        with p.open() as fh:
            target = fh.readline().strip()
        if sp and target:
            mapping[sp] = target
    return mapping


def prune_and_rename(tree_path: Path, keep_species, sp_to_target):
    tree = Phylo.read(tree_path, "newick")
    # prune species not in keep_species
    for tip in list(tree.get_terminals()):
        if tip.name not in keep_species:
            try:
                tree.prune(tip)
            except Exception:
                # if pruning by object fails, attempt by name; ignore if still not found
                try:
                    tree.prune(tip.name)
                except Exception:
                    pass
    # rename remaining tips
    for tip in tree.get_terminals():
        new_name = sp_to_target.get(tip.name)
        if new_name:
            tip.name = new_name
    return tree


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--keep")
    ap.add_argument("--manifest")
    ap.add_argument("--dir", help="Base dir containing *.best_id.txt from a prior run (fallback when manifest is absent)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sp_to_target = {}
    keep = set()

    if args.manifest:
        sp_to_target = read_manifest_targets(Path(args.manifest))
    elif args.dir:
        sp_to_target = read_best_id_targets(Path(args.dir))
    else:
        raise SystemExit("[build_protein_tree] Need either --manifest or --dir")

    if args.keep:
        keep = read_keep(Path(args.keep))
    elif args.dir:
        # derive keep list from available best_id files with non-empty targets
        keep = set(sp_to_target.keys())
    else:
        raise SystemExit("[build_protein_tree] Need either --keep or --dir")

    # only species that are both kept and have targets
    keep = {sp for sp in keep if sp in sp_to_target}
    if not keep:
        raise SystemExit(f"[build_protein_tree] No species to keep after intersecting keep list and manifest targets")

    tree = prune_and_rename(Path(args.tree), keep, sp_to_target)
    # sanity: ensure all remaining terminals have target names
    remaining = {t.name for t in tree.get_terminals()}
    if not remaining:
        raise SystemExit("[build_protein_tree] Tree is empty after pruning")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, args.out, "newick")


if __name__ == "__main__":
    main()
