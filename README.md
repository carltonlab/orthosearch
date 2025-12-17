# OrthoSearch Snakemake Pipeline (MMseqs2 + PAL2NAL)

This Snakemake workflow:
1) Takes a **query protein FASTA** (`query/query.protein.faa`) as input (single sequence).
2) Searches it against each species proteome in `data/protein/` using **MMseqs2**.
3) Chooses the best candidate per species (with QC + top-K table), then extracts:
   - per-species protein FASTA
   - per-species CDS FASTA (matching record ID)
4) Produces combined FASTAs and a **codon-aware CDS alignment** using **MAFFT** (protein MSA) + **PAL2NAL** (back-translation).

## Input layout

Place your files:

- `data/protein/*.protein.fa`
- `data/cds/*.cds.fa`

Naming convention supported:

- `genus_species.info.protein.fa` or `genus_species.protein.fa`
- `genus_species.info.cds.fa` or `genus_species.cds.fa`

The species key is the substring **before the first dot** in each protein filename.
For each species, the workflow requires **exactly one** matching protein file and **exactly one** matching CDS file.
Missing or non-unique matches cause a **hard error** (fail-fast).

Query:
- `query/query.protein.faa` (single-sequence FASTA)

Optional overrides:
- `overrides/overrides.tsv` with columns: `species<TAB>chosen_protein_id`

## Run

Create a small venv for Snakemake (optional):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip snakemake
```

Then run:

```bash
snakemake --use-conda --cores 16 codon_alignment
```

Outputs:
- `results/combined.proteins.faa`
- `results/combined.cds.fna`
- `results/protein.aln.faa`
- `results/cds.codon_aware.aln.fasta`
- per-species: `orthologs/{species}.qc.tsv`, `orthologs/{species}.mmseqs_topk.tsv`, plus extracted FASTAs.

## Notes

- Species with **no passing hit** are skipped (their extracted FASTAs are empty).
- If *no* species produces a passing hit overall, the workflow fails.
