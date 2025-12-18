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

Filename grammar supported:

- `genus_species.protein.fa` OR `genus_species.<info>.protein.fa`
- `genus_species.cds.fa` OR `genus_species.<info>.cds.fa`

Species key is the substring **before the first dot** in each protein filename.

**Fail-fast**: for each species, the workflow requires exactly one matching protein file and exactly one matching CDS file.
Matching is *dot-anchored* to avoid sp2/sp24 collisions:
- `{species}.protein.fa` or `{species}.*.protein.fa` (analogous for CDS).

Query:
- `query/query.protein.faa`

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

1) Put your new query sequence (single FASTA record) at `query/query.protein.faa` (or point `paths.query_protein` in `config.yaml` to a different file).  
2) Rerun:
   ```bash
   snakemake --use-conda --cores 16 --rerun-incomplete
   ```
   Snakemake will detect the changed query FASTA and rebuild only the MMseqs-dependent steps and everything downstream.

## Outputs

Unfiltered (all best hits):
- `results/combined.proteins.all.faa`
- `results/combined.cds.all.fna`
- `results/protein.all.aln.faa`
- `results/cds.codon_aware.all.aln.fasta`

Filtered (based on `config.yaml: final_filter` thresholds):
- `results/combined.proteins.filtered.faa`
- `results/combined.cds.filtered.fna`
- `results/protein.filtered.aln.faa`
- `results/cds.codon_aware.filtered.aln.fasta`

Manifests:
- `results/manifest.all.tsv` (per-species metrics for chosen hit)
- `results/keep_species.txt` (species passing final_filter)

Per species:
- `orthologs/{species}.mmseqs_topk.tsv`
- `orthologs/{species}.qc.tsv`
- `orthologs/{species}.best_id.txt`
- `orthologs/{species}.protein.faa`
- `orthologs/{species}.cds.fna`

## Iterating thresholds without re-running MMseqs

Edit `config.yaml: final_filter` and rerun only the cheap downstream steps:

```bash
snakemake --use-conda --cores 16 results/cds.codon_aware.filtered.aln.fasta \
  --forcerun build_manifest_and_keep \
            combined_proteins_filtered \
            combined_cds_filtered \
            mafft_protein_alignment_filtered \
            pal2nal_codon_alignment_filtered
```

This keeps all MMseqs searches intact and refreshes the manifest, filtered FASTAs, and filtered alignments with the new thresholds.
