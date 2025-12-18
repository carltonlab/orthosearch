import os, glob

configfile: "config.yaml"

PROTEIN_DIR = config["paths"]["protein_dir"]
CDS_DIR = config["paths"]["cds_dir"]
QUERY_PROTEIN = config["paths"]["query_protein"]
OVERRIDES = "overrides/overrides.tsv"

def species_from_filename(fn: str) -> str:
    return os.path.basename(fn).split(".", 1)[0]

def require_unique_species_file(dirpath: str, species: str, type_token: str) -> str:
    patterns = [
        os.path.join(dirpath, f"{species}.{type_token}.fa"),
        os.path.join(dirpath, f"{species}.*.{type_token}.fa"),
    ]
    matches = []
    for pat in patterns:
        matches.extend(glob.glob(pat))
    matches = sorted(set(matches))
    if len(matches) == 0:
        raise ValueError(f"Missing {type_token} file for {species}: tried {patterns}")
    if len(matches) > 1:
        raise ValueError(f"Non-unique {type_token} file for {species}: {matches}")
    return matches[0]

protein_files = sorted(glob.glob(os.path.join(PROTEIN_DIR, "*.protein.fa")))
if not protein_files:
    raise ValueError(f"No protein files found in {PROTEIN_DIR} matching *.protein.fa")

SPECIES = sorted({species_from_filename(p) for p in protein_files})

PROT_BY_SPECIES = {sp: require_unique_species_file(PROTEIN_DIR, sp, "protein") for sp in SPECIES}
CDS_BY_SPECIES  = {sp: require_unique_species_file(CDS_DIR, sp, "cds") for sp in SPECIES}

rule all:
    input:
        "results/cds.codon_aware.filtered.aln.fasta"

rule codon_alignment_all:
    input:
        "results/cds.codon_aware.all.aln.fasta"

rule codon_alignment_filtered:
    input:
        "results/cds.codon_aware.filtered.aln.fasta"

rule mmseqs_createdb_query:
    input:
        QUERY_PROTEIN
    output:
        "work/mmseqs/queryDB"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "mkdir -p work/mmseqs && mmseqs createdb {input} {output}"

rule mmseqs_createdb_species:
    input:
        lambda wc: PROT_BY_SPECIES[wc.sp]
    output:
        "work/mmseqs/{sp}.protDB"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "mkdir -p work/mmseqs && mmseqs createdb {input} {output}"

rule mmseqs_search:
    input:
        qdb="work/mmseqs/queryDB",
        tdb="work/mmseqs/{sp}.protDB"
    output:
        "work/mmseqs/{sp}.alnDB"
    threads:
        config["threads"]["mmseqs"]
    conda:
        "workflow/envs/search.yaml"
    shell:
        "mkdir -p work/mmseqs/tmp_{wildcards.sp} && mmseqs search {input.qdb} {input.tdb} {output} work/mmseqs/tmp_{wildcards.sp} --threads {threads}"

rule mmseqs_tsv:
    input:
        qdb="work/mmseqs/queryDB",
        tdb="work/mmseqs/{sp}.protDB",
        adb="work/mmseqs/{sp}.alnDB"
    output:
        "work/mmseqs/{sp}.hits.tsv"
    conda:
        "workflow/envs/search.yaml"
    shell:
        r'''
        set -euo pipefail
        tmp="{output}.body"
        rm -f "$tmp"
        mmseqs convertalis {input.qdb} {input.tdb} {input.adb} "$tmp"           --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qlen,tlen"           2>/dev/null || true
        printf "query	target	fident	alnlen	mismatch	gapopen	qstart	qend	tstart	tend	evalue	bits	qcov	tcov	qlen	tlen
" > {output}
        if [ -s "$tmp" ]; then
          cat "$tmp" >> {output}
        fi
        '''

rule pick_best_hit:
    input:
        hits="work/mmseqs/{sp}.hits.tsv"
    output:
        best="orthologs/{sp}.best_id.txt",
        topk="orthologs/{sp}.mmseqs_topk.tsv",
        qc="orthologs/{sp}.qc.tsv"
    params:
        top_k=config["search"]["top_k"],
        qc_min_qcov=config["search"]["qc_min_query_cov"],
        qc_min_tcov=config["search"]["qc_min_target_cov"],
        qc_min_len_ratio=config["search"]["qc_min_len_ratio"],
        qc_max_len_ratio=config["search"]["qc_max_len_ratio"],
        overrides=OVERRIDES
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/pick_best_hit.py --hits {input.hits} --species {wildcards.sp} --top-k {params.top_k} --qc-min-qcov {params.qc_min_qcov} --qc-min-tcov {params.qc_min_tcov} --qc-min-len-ratio {params.qc_min_len_ratio} --qc-max-len-ratio {params.qc_max_len_ratio} --overrides {params.overrides} --out-best-id {output.best} --out-topk {output.topk} --out-qc {output.qc}"

rule extract_protein:
    input:
        fasta=lambda wc: PROT_BY_SPECIES[wc.sp],
        best="orthologs/{sp}.best_id.txt"
    output:
        "orthologs/{sp}.protein.faa"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/extract_fasta_record.py --fasta {input.fasta} --id-file {input.best} --out {output} --strict"

rule extract_cds:
    input:
        fasta=lambda wc: CDS_BY_SPECIES[wc.sp],
        best="orthologs/{sp}.best_id.txt"
    output:
        "orthologs/{sp}.cds.fna"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/extract_fasta_record.py --fasta {input.fasta} --id-file {input.best} --out {output} --strict"

rule combined_proteins_all:
    input:
        expand("orthologs/{sp}.protein.faa", sp=SPECIES)
    output:
        "results/combined.proteins.all.faa"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/combine_fasta.py --inputs {input} --out {output} --require-at-least-one"

rule combined_cds_all:
    input:
        expand("orthologs/{sp}.cds.fna", sp=SPECIES)
    output:
        "results/combined.cds.all.fna"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/combine_fasta.py --inputs {input} --out {output} --require-at-least-one"

rule mafft_protein_alignment_all:
    input:
        "results/combined.proteins.all.faa"
    output:
        "results/protein.all.aln.faa"
    threads:
        config["threads"]["mafft"]
    params:
        mafft_args=config["alignment"]["mafft_args"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "mafft {params.mafft_args} --thread {threads} {input} > {output}"

rule pal2nal_codon_alignment_all:
    input:
        prot_aln="results/protein.all.aln.faa",
        cds="results/combined.cds.all.fna"
    output:
        "results/cds.codon_aware.all.aln.fasta"
    params:
        pal2nal_output=config["alignment"]["pal2nal_output"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "pal2nal.pl {input.prot_aln} {input.cds} -output {params.pal2nal_output} > {output}"

rule build_manifest_and_keep:
    input:
        topk=expand("orthologs/{sp}.mmseqs_topk.tsv", sp=SPECIES),
        qc=expand("orthologs/{sp}.qc.tsv", sp=SPECIES)
    output:
        manifest="results/manifest.all.tsv",
        keep="results/keep_species.txt"
    params:
        min_qcov=config["final_filter"]["min_qcov"],
        min_tcov=config["final_filter"]["min_tcov"],
        min_len_ratio=config["final_filter"]["min_len_ratio"],
        max_len_ratio=config["final_filter"]["max_len_ratio"],
        min_bits=config["final_filter"]["min_bits"]
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/build_manifest_and_keep.py --topk-tsvs {input.topk} --qc-tsvs {input.qc} --out-manifest {output.manifest} --out-keep {output.keep} --min-qcov {params.min_qcov} --min-tcov {params.min_tcov} --min-len-ratio {params.min_len_ratio} --max-len-ratio {params.max_len_ratio} --min-bits {params.min_bits}"

rule combined_proteins_filtered:
    input:
        keep="results/keep_species.txt",
        manifest="results/manifest.all.tsv",
        proteins=expand("orthologs/{sp}.protein.faa", sp=SPECIES)
    output:
        "results/combined.proteins.filtered.faa"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/filter_combined_by_keep.py --keep-species {input.keep} --orthologs-dir orthologs --suffix protein.faa --out {output} --require-at-least-one"

rule combined_cds_filtered:
    input:
        keep="results/keep_species.txt",
        manifest="results/manifest.all.tsv",
        cds=expand("orthologs/{sp}.cds.fna", sp=SPECIES)
    output:
        "results/combined.cds.filtered.fna"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/filter_combined_by_keep.py --keep-species {input.keep} --orthologs-dir orthologs --suffix cds.fna --out {output} --require-at-least-one"

rule mafft_protein_alignment_filtered:
    input:
        "results/combined.proteins.filtered.faa"
    output:
        "results/protein.filtered.aln.faa"
    threads:
        config["threads"]["mafft"]
    params:
        mafft_args=config["alignment"]["mafft_args"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "mafft {params.mafft_args} --thread {threads} {input} > {output}"

rule pal2nal_codon_alignment_filtered:
    input:
        prot_aln="results/protein.filtered.aln.faa",
        cds="results/combined.cds.filtered.fna"
    output:
        "results/cds.codon_aware.filtered.aln.fasta"
    params:
        pal2nal_output=config["alignment"]["pal2nal_output"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "pal2nal.pl {input.prot_aln} {input.cds} -output {params.pal2nal_output} > {output}"
