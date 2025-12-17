import os, glob

configfile: "config.yaml"

PROTEIN_DIR = config["paths"]["protein_dir"]
CDS_DIR = config["paths"]["cds_dir"]
QUERY_PROTEIN = config["paths"]["query_protein"]
OVERRIDES = "overrides/overrides.tsv"

def species_from_filename(fn: str) -> str:
    return os.path.basename(fn).split(".", 1)[0]

def require_unique_species_file(dirpath, species, suffix):
    patterns = [
        os.path.join(dirpath, f"{species}.{suffix}"),
        os.path.join(dirpath, f"{species}.*.{suffix}")
    ]
    matches = []
    for pat in patterns:
        matches.extend(glob.glob(pat))
    matches = sorted(set(matches))
    if len(matches) == 0:
        raise ValueError(f"Missing {suffix} file for species {species}")
    if len(matches) > 1:
        raise ValueError(f"Non-unique {suffix} file for species {species}: {matches}")
    return matches[0]


#def require_unique(pattern: str, what: str) -> str:
#    matches = sorted(glob.glob(pattern))
#    if len(matches) == 0:
#        raise ValueError(f"Missing {what}: pattern={pattern}")
#    if len(matches) > 1:
#        raise ValueError(f"Non-unique {what}: pattern={pattern} matches={matches}")
#    return matches[0]

protein_files = sorted(glob.glob(os.path.join(PROTEIN_DIR, "*.protein.fa")))
if not protein_files:
    raise ValueError(f"No protein files found in {PROTEIN_DIR} matching *.protein.fa")

SPECIES = sorted({species_from_filename(p) for p in protein_files})

PROT_BY_SPECIES = {
    sp: require_unique_species_file(PROTEIN_DIR, sp, "protein.fa")
    for sp in SPECIES
}

CDS_BY_SPECIES = {
    sp: require_unique_species_file(CDS_DIR, sp, "cds.fa")
    for sp in SPECIES
}

rule all:
    input:
        "results/cds.codon_aware.aln.fasta"

rule codon_alignment:
    input:
        "results/cds.codon_aware.aln.fasta"

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
        "mmseqs convertalis {input.qdb} {input.tdb} {input.adb} {output} --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,qlen,tlen'"

rule pick_best_hit:
    input:
        hits="work/mmseqs/{sp}.hits.tsv"
    output:
        best="orthologs/{sp}.best_id.txt",
        topk="orthologs/{sp}.mmseqs_topk.tsv",
        qc="orthologs/{sp}.qc.tsv"
    params:
        top_k=config["search"]["top_k"],
        min_qcov=config["search"]["min_query_cov"],
        min_tcov=config["search"]["min_target_cov"],
        min_fident=config["search"]["min_fident"],
        min_bits=config["search"]["min_bits"],
        min_len_ratio=config["search"]["min_len_ratio"],
        max_len_ratio=config["search"]["max_len_ratio"],
        overrides=OVERRIDES
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/pick_best_hit.py --hits {input.hits} --species {wildcards.sp} --top-k {params.top_k} --min-qcov {params.min_qcov} --min-tcov {params.min_tcov} --min-fident {params.min_fident} --min-bits {params.min_bits} --min-len-ratio {params.min_len_ratio} --max-len-ratio {params.max_len_ratio} --overrides {params.overrides} --out-best-id {output.best} --out-topk {output.topk} --out-qc {output.qc}"

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

rule combined_proteins:
    input:
        expand("orthologs/{sp}.protein.faa", sp=SPECIES)
    output:
        "results/combined.proteins.faa"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/combine_fasta.py --inputs {input} --out {output} --require-at-least-one"

rule combined_cds:
    input:
        expand("orthologs/{sp}.cds.fna", sp=SPECIES)
    output:
        "results/combined.cds.fna"
    conda:
        "workflow/envs/search.yaml"
    shell:
        "python workflow/scripts/combine_fasta.py --inputs {input} --out {output} --require-at-least-one"

rule mafft_protein_alignment:
    input:
        "results/combined.proteins.faa"
    output:
        "results/protein.aln.faa"
    threads:
        config["threads"]["mafft"]
    params:
        mafft_args=config["alignment"]["mafft_args"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "mafft {params.mafft_args} --thread {threads} {input} > {output}"

rule pal2nal_codon_alignment:
    input:
        prot_aln="results/protein.aln.faa",
        cds="results/combined.cds.fna"
    output:
        "results/cds.codon_aware.aln.fasta"
    params:
        pal2nal_output=config["alignment"]["pal2nal_output"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        "pal2nal.pl {input.prot_aln} {input.cds} -output {params.pal2nal_output} > {output}"
