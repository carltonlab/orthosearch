# OrthoSearch Snakemake Pipeline v2 (MMseqs2 + PAL2NAL + downstream filtering)

This workflow:
1) Takes a **query protein FASTA** (`query/query.protein.faa`, single sequence).
2) Searches it against each species proteome in `data/protein/` using **MMseqs2**.
3) Chooses the **best hit per species** (by bitscore/score) and extracts:
   - per-species protein FASTA
   - per-species CDS FASTA (matching record ID)
4) Builds combined FASTAs and a codon-aware CDS alignment using **MAFFT** + **PAL2NAL**.
5) Supports **cheap, downstream filtering** of final outputs (without re-running MMseqs) based on
   qcov/tcov/len_ratio/bits, etc.

## Input layout

- `data/protein/*.protein.fa`
- `data/cds/*.cds.fa`
- `query/*.fa|*.fasta|*.faa` (one or many files; each may contain one or many query proteins)

Filename grammar supported:

- `genus_species.protein.fa` OR `genus_species.<info>.protein.fa`
- `genus_species.cds.fa` OR `genus_species.<info>.cds.fa`

Species key is the substring **before the first dot** in each protein filename.

**Fail-fast**: for each species, the workflow requires exactly one matching protein file and exactly one matching CDS file.
Matching is *dot-anchored* to avoid sp2/sp24 collisions:
- `{species}.protein.fa` or `{species}.*.protein.fa` (analogous for CDS).

Query:
- Place any number of query sequences under `query/` (FASTA). Headers must be unique across all query files. The pipeline creates per-query outputs using a sanitized version of each FASTA header as the directory name. If two headers sanitize to the same name, the workflow fails fast.

Optional:
- `overrides/overrides.tsv` with columns: `species<TAB>chosen_protein_id`

## Run

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip snakemake

snakemake --use-conda --cores 16
```

## Running a new query
Drop additional FASTA files or sequences into `query/`. Keep headers unique. Then rerun:
```bash
snakemake --use-conda --cores 16 --rerun-incomplete
```
Outputs are written under `results/<query_id>/` and `orthologs/<query_id>/`, where `<query_id>` is the sanitized FASTA header (non-alphanumeric characters replaced with `_`). A collision (two headers mapping to the same sanitized name) will fail fast.

## Outputs

Unfiltered (all best hits):
- `results/<query_id>/combined.proteins.all.faa`
- `results/<query_id>/combined.cds.all.fna`
- `results/<query_id>/protein.all.aln.faa`
- `results/<query_id>/cds.codon_aware.all.aln.fasta`

Filtered (based on `config.yaml: final_filter` thresholds):
- `results/<query_id>/combined.proteins.filtered.faa`
- `results/<query_id>/combined.cds.filtered.fna`
- `results/<query_id>/protein.filtered.aln.faa`
- `results/<query_id>/cds.codon_aware.filtered.aln.fasta`

Rejected candidates (fail coverage/length-ratio thresholds; for manual rescue):
- `results/<query_id>/combined.proteins.rejected.faa`

Manifests:
- `results/<query_id>/manifest.all.tsv` (per-species metrics for chosen hit)
- `results/<query_id>/keep_species.txt` (species passing final_filter)

Per species:
- `orthologs/<query_id>/{species}.mmseqs_topk.tsv`
- `orthologs/<query_id>/{species}.qc.tsv`
- `orthologs/<query_id>/{species}.best_id.txt`
- `orthologs/<query_id>/{species}.protein.faa`
- `orthologs/<query_id>/{species}.cds.fna`

## Iterating thresholds without re-running MMseqs

Edit `config.yaml: final_filter` and rerun only the cheap downstream steps for one query:
```bash
TARGET_QID=<query_id>
snakemake --use-conda --cores 16 results/${TARGET_QID}/cds.codon_aware.filtered.aln.fasta \
  --forcerun build_manifest_and_keep \
            combined_proteins_filtered \
            combined_cds_filtered \
            mafft_protein_alignment_filtered \
            pal2nal_codon_alignment_filtered
```

This keeps all MMseqs searches intact and refreshes the manifest, filtered FASTAs, and filtered alignments with the new thresholds. To refresh everything for all queries, drop the target path (or use a glob) and keep the same `--forcerun` list.
