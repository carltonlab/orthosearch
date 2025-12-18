import os, glob, re

configfile: "config.yaml"

PROTEIN_DIR = config["paths"]["protein_dir"]
CDS_DIR = config["paths"]["cds_dir"]
QUERY_DIR = config["paths"]["query_dir"]
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

def slug_id(s: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", s.strip())
    if not safe:
        raise ValueError(f"Empty/unsafe query id after sanitizing: {s!r}")
    return safe

def parse_query_headers(query_dir: str):
    files = sorted(sum([glob.glob(os.path.join(query_dir, ext)) for ext in ("*.fa", "*.fasta", "*.faa")], []))
    if not files:
        raise ValueError(f"No query FASTA files found in {query_dir}")
    seen_ids = {}
    for path in files:
        cur_id = ""
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    if not header:
                        raise ValueError(f"Empty FASTA header in {path}")
                    if header in seen_ids:
                        raise ValueError(f"Duplicate query id '{header}' in {path} and {seen_ids[header]}")
                    seen_ids[header] = path
    if not seen_ids:
        raise ValueError(f"No sequences found across query files: {files}")
    # Build mapping from safe dir name -> (orig id, source file)
    by_safe = {}
    for orig, path in seen_ids.items():
        safe = slug_id(orig)
        if safe in by_safe:
            raise ValueError(f"Different query ids collide to same directory name '{safe}': {by_safe[safe]['orig']} vs {orig}")
        by_safe[safe] = {"orig": orig, "file": path}
    return by_safe

QUERY_INFO = parse_query_headers(QUERY_DIR)
QUERY_IDS = sorted(QUERY_INFO.keys())

def query_file_for(wc):
    return QUERY_INFO[wc.qid]["file"]

def query_orig_id_for(wc):
    return QUERY_INFO[wc.qid]["orig"]

rule all:
    input:
        expand("results/{qid}/cds.codon_aware.filtered.aln.fasta", qid=QUERY_IDS),
        expand("results/{qid}/combined.proteins.rejected.faa", qid=QUERY_IDS)

rule codon_alignment_all:
    input:
        expand("results/{qid}/cds.codon_aware.all.aln.fasta", qid=QUERY_IDS)

rule codon_alignment_filtered:
    input:
        expand("results/{qid}/cds.codon_aware.filtered.aln.fasta", qid=QUERY_IDS)

rule mmseqs_createdb_query:
    input:
        lambda wc: query_file_for(wc)
    params:
        orig_id=lambda wc: query_orig_id_for(wc)
    output:
        "work/queries/{qid}.faa",
        "work/mmseqs/{qid}/queryDB"
    conda:
        "workflow/envs/search.yaml"
    shell:
        r'''
        set -euo pipefail
        mkdir -p work/queries work/mmseqs/{wildcards.qid}
        # Extract the specific query sequence by id into a single-seq fasta
        python - "{params.orig_id}" {input} {output[0]} <<'PY'
import sys, os
orig_id = sys.argv[1]
src = sys.argv[2]
out = sys.argv[3]
found = False
seq = []
with open(src) as fh:
    cur = None
    for line in fh:
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            continue
        if cur == orig_id:
            seq.append(line.rstrip())
            found = True
if not found:
    raise SystemExit("[query extract] id %s not found in %s" % (orig_id, src))
with open(out, "w") as o:
    o.write(">" + orig_id + "\n" + "\n".join(seq) + "\n")
PY
        mmseqs createdb {output[0]} {output[1]}
        '''

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
        qdb="work/mmseqs/{qid}/queryDB",
        tdb="work/mmseqs/{sp}.protDB"
    output:
        "work/mmseqs/{qid}/{sp}.alnDB"
    threads:
        config["threads"]["mmseqs"]
    conda:
        "workflow/envs/search.yaml"
    shell:
        "mkdir -p work/mmseqs/tmp_{wildcards.qid}_{wildcards.sp} && mmseqs search {input.qdb} {input.tdb} {output} work/mmseqs/tmp_{wildcards.qid}_{wildcards.sp} --threads {threads}"

rule mmseqs_tsv:
    input:
        qdb="work/mmseqs/{qid}/queryDB",
        tdb="work/mmseqs/{sp}.protDB",
        adb="work/mmseqs/{qid}/{sp}.alnDB"
    output:
        "work/mmseqs/{qid}/{sp}.hits.tsv"
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
        hits="work/mmseqs/{qid}/{sp}.hits.tsv"
    output:
        best="orthologs/{qid}/{sp}.best_id.txt",
        topk="orthologs/{qid}/{sp}.mmseqs_topk.tsv",
        qc="orthologs/{qid}/{sp}.qc.tsv"
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
        best="orthologs/{qid}/{sp}.best_id.txt"
    output:
        "orthologs/{qid}/{sp}.protein.faa"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/extract_fasta_record.py --fasta {input.fasta} --id-file {input.best} --out {output} --strict"

rule extract_cds:
    input:
        fasta=lambda wc: CDS_BY_SPECIES[wc.sp],
        best="orthologs/{qid}/{sp}.best_id.txt"
    output:
        "orthologs/{qid}/{sp}.cds.fna"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/extract_fasta_record.py --fasta {input.fasta} --id-file {input.best} --out {output} --strict"

rule combined_proteins_all:
    input:
        lambda wc: expand(f"orthologs/{wc.qid}" + "/{sp}.protein.faa", sp=SPECIES)
    output:
        "results/{qid}/combined.proteins.all.faa"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/combine_fasta.py --inputs {input} --out {output} --require-at-least-one"

rule combined_cds_all:
    input:
        lambda wc: expand(f"orthologs/{wc.qid}" + "/{sp}.cds.fna", sp=SPECIES)
    output:
        "results/{qid}/combined.cds.all.fna"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/combine_fasta.py --inputs {input} --out {output} --require-at-least-one"

rule mafft_protein_alignment_all:
    input:
        "results/{qid}/combined.proteins.all.faa"
    output:
        "results/{qid}/protein.all.aln.faa"
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
        prot_aln="results/{qid}/protein.all.aln.faa",
        cds="results/{qid}/combined.cds.all.fna"
    output:
        "results/{qid}/cds.codon_aware.all.aln.fasta"
    params:
        pal2nal_output=config["alignment"]["pal2nal_output"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "pal2nal.pl {input.prot_aln} {input.cds} -output {params.pal2nal_output} > {output}"

rule build_manifest_and_keep:
    input:
        topk=lambda wc: expand(f"orthologs/{wc.qid}" + "/{sp}.mmseqs_topk.tsv", sp=SPECIES),
        qc=lambda wc: expand(f"orthologs/{wc.qid}" + "/{sp}.qc.tsv", sp=SPECIES)
    output:
        manifest="results/{qid}/manifest.all.tsv",
        keep="results/{qid}/keep_species.txt"
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

rule combined_proteins_rejected:
    input:
        manifest="results/{qid}/manifest.all.tsv",
        proteins=lambda wc: expand(f"orthologs/{wc.qid}" + "/{sp}.protein.faa", sp=SPECIES)
    output:
        "results/{qid}/combined.proteins.rejected.faa"
    params:
        min_qcov=config["final_filter"]["min_qcov"],
        min_tcov=config["final_filter"]["min_tcov"],
        min_len_ratio=config["final_filter"]["min_len_ratio"],
        max_len_ratio=config["final_filter"]["max_len_ratio"]
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/collect_rejected_by_manifest.py --manifest {input.manifest} --orthologs-dir orthologs/{wildcards.qid} --suffix protein.faa --out {output} --min-qcov {params.min_qcov} --min-tcov {params.min_tcov} --min-len-ratio {params.min_len_ratio} --max-len-ratio {params.max_len_ratio}"

rule combined_proteins_filtered:
    input:
        keep="results/{qid}/keep_species.txt",
        manifest="results/{qid}/manifest.all.tsv",
        proteins=lambda wc: expand(f"orthologs/{wc.qid}" + "/{sp}.protein.faa", sp=SPECIES)
    output:
        "results/{qid}/combined.proteins.filtered.faa"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/filter_combined_by_keep.py --keep-species {input.keep} --orthologs-dir orthologs/{wildcards.qid} --suffix protein.faa --out {output} --require-at-least-one"

rule combined_cds_filtered:
    input:
        keep="results/{qid}/keep_species.txt",
        manifest="results/{qid}/manifest.all.tsv",
        cds=lambda wc: expand(f"orthologs/{wc.qid}" + "/{sp}.cds.fna", sp=SPECIES)
    output:
        "results/{qid}/combined.cds.filtered.fna"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/filter_combined_by_keep.py --keep-species {input.keep} --orthologs-dir orthologs/{wildcards.qid} --suffix cds.fna --out {output} --require-at-least-one"

rule mafft_protein_alignment_filtered:
    input:
        "results/{qid}/combined.proteins.filtered.faa"
    output:
        "results/{qid}/protein.filtered.aln.faa"
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
        prot_aln="results/{qid}/protein.filtered.aln.faa",
        cds="results/{qid}/combined.cds.filtered.fna"
    output:
        "results/{qid}/cds.codon_aware.filtered.aln.fasta"
    params:
        pal2nal_output=config["alignment"]["pal2nal_output"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "pal2nal.pl {input.prot_aln} {input.cds} -output {params.pal2nal_output} > {output}"
