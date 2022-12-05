import itertools as it
from pathlib import Path
def add_ext(p, *args):
    base = str(p)
    for x in args:
        base += str(x)
    return Path(base)
    # return Path(str(p)+str(e))
# bcftools view -S na.samples -q 0.05 -Q 0.95 -m 2 -M 2 -v snps -Oz -o vcf/na.chr5.vcf.gz --threads 4 /emc/kristian/1kg_20220422/1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gza
P = config["python"]
VCF = config["vcf"]
HAPLONET = config["haplonet"]
MASK_ANCESTRY = config["mask_anc"]
PLOTTER = config["plotter"]
BCFTOOLS = config["bcftools"]
MAF_FILTER = config["filter"]
K = config["K"]
ALPHA = config["alpha"]
SEED = config["seed"]
SUBSPLIT_LST = config["subsplit"]
SUBSPLIT_MAX = max(SUBSPLIT_LST)
SUBSPLIT_ARG = "" if SUBSPLIT_MAX == 0 else f"--subsplit {SUBSPLIT_MAX}" 
SAMPLES = config["samples"]

PC_ITERATIVE = config["iterative"]
MAX_MISSINGNESS = config["max_missingness"]
HET_HOM = "--het_mask" if config["het_hom"] == "het" else ""
POP_LABELS = config["labels"]

res = Path(config["output"])

pcs = ["PC1", "PC2", "PC3"]
pcs_comb = list(it.combinations(pcs, 2))
t_wc = ["prob.npy", "path"]
k_seed = f"K{K}_s{SEED}"

with open(config["chromlist"], 'r') as fh:
    CHROMLIST = [x.rstrip() for x in fh]

MASKERTYPES = ["posterior", "decoding"]

wildcard_constraints:
    t = "|".join(t_wc),
    chrom = "|".join(CHROMLIST),
    masktype = "|".join(MASKERTYPES),
    maf_filter = "|".join(map(str, MAF_FILTER))

rule all:
    input:
        expand(res / "pca" / k_seed / "MAF{maf_filter}" / 
            "sub{subsplit}_{masktype}_{pc_comb[0]}_{pc_comb[1]}.png",
            maf_filter = MAF_FILTER,
            subsplit = SUBSPLIT_LST,
            masktype = MASKERTYPES,
            pc_comb = pcs_comb
            )

rule pca_decoding:
    input:
        expand(res / "pca" / k_seed / "MAF{maf_filter}" / 
            "sub{subsplit}_decoding_{pc_comb[0]}_{pc_comb[1]}.png",
            maf_filter = MAF_FILTER,
            subsplit = SUBSPLIT_LST,
            pc_comb = pcs_comb
            )

rule pca_posterior:
    input:
        expand(res / "pca" / k_seed / "MAF{maf_filter}" / 
            "sub{subsplit}_posterior_{pc_comb[0]}_{pc_comb[1]}.png",
            maf_filter = MAF_FILTER,
            subsplit = SUBSPLIT_LST,
            pc_comb = pcs_comb
            )

rule extract_samples:
    input:
        VCF
    output:
        res / "vcf" / "{chrom}.vcf.gz"
    threads: 4
    shell:
        "{BCFTOOLS} view --threads {threads} "
        "-S {SAMPLES} -q 0.05 -Q 0.95 -r {wildcards.chrom} "
        "-m 2 -M 2 -v snps -Oz -o {output} {input}"

rule train:
    input:
        res / "vcf" / "{chrom}.vcf.gz"
    output:
        res / "haplonet" / "{chrom}" / "temp.loglike.npy",
    log:
        res / "haplonet" / "{chrom}" / "temp.loglike.npy.log",  
    threads: 4
    params:
        out = lambda wc, output: output[0][:-12],
    shell:
        "{HAPLONET} train --vcf {input} "
        "--out {params.out} {SUBSPLIT_ARG} "
        "--threads {threads} > {log}"

rule make_subsplits:
    input:
        nosplit = res / "haplonet" / "{chrom}" / f"temp.loglike.npy",  
    output:
        res / "haplonet_split" / "{chrom}" / "sub{subsplit}.loglike.npy"
    run:
        import numpy as np
        outdir = os.path.dirname(output[0])
        shell("mkdir -p {outdir}")
        maxsplit = input.nosplit.replace(".loglike.npy", ".split.loglike.npy")
        if int(wildcards.subsplit) == 0:
            shell("cp {input.nosplit} {output[0]}")
        elif int(wildcards.subsplit) == SUBSPLIT_MAX:
            shell("cp {maxsplit} {output[0]}")
        else:
            out = np.load(maxsplit)
            temp_start = SUBSPLIT_MAX
            while temp_start != int(wildcards.subsplit):
                out = out[0::2, :, :] + out[1::2, :, :]
                temp_start -= temp_start // 2
            np.save(output[0], out)

# rule gen_filelist:
#     input:
#         expand(res / "haplonet_split" / "{chrom}" / "sub{subsplit}.loglike.npy"
#             chrom = CHROMLIST, allow_missing=True)
#     output:
#         res / "haplonet_split" / "sub{subsplit}.filelist"
#     run:
#        with open(output[0], 'w') as fh:
#             for x in input:
#                 print(x, file=fh)

rule concat:
    input:
        expand(res / "haplonet_split" / "{chrom}" / "sub{subsplit}.loglike.npy",
            chrom = CHROMLIST, allow_missing=True)
    output:
        protected(res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy") 
    log:
        res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy.log" 
    shell:
         "{P} ./scripts/concat_npyV2.py {output} {input} > {log}"
        
rule admix:
    input:
        res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy" 
    output:
        res / "admix" / k_seed / "sub{subsplit}.f.npy",
        res / "admix" / k_seed / "sub{subsplit}.q"
    params:
        out = lambda wc, output: output[0][:-6]
    threads: 20
    shell:
        "{HAPLONET} admix --like {input} "
        "--K {K} --seed {SEED} --out {params.out} "
        "--threads {threads}"

def set_alpha(subsplit):
    if str(ALPHA).lower() == "est":
        return "--alpha_save"
    elif int(subsplit) == 0:
        return f"--no_optim --alpha_save --alpha {ALPHA}"
    else:
        scaled_alpha = float(ALPHA) / float(subsplit)
        return f"--no_optim --alpha_save --alpha {scaled_alpha}"

rule fatass:
    input:
        l = res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
        f = res / "admix" / k_seed / "sub{subsplit}.f.npy",
        q = res / "admix" / k_seed / "sub{subsplit}.q"
    output:
        res / "fatass" / k_seed / "sub{subsplit}.prob.npy",
        res / "fatass" / k_seed / "sub{subsplit}.path",
        res / "fatass" / k_seed / "sub{subsplit}.alpha"
    params:
        out = lambda wc, output: output[0][:-9],
        alpha = lambda wc: set_alpha(wc.subsplit)
    threads: 10
    shell:
        "{HAPLONET} fatash {params.alpha} "
        "--like {input.l} "
        "--prop {input.q} --freq {input.f} "
        "--out {params.out} "
        "--threads {threads}"

rule mask_posterior:
    input:
        l = res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
        masker = res / "fatass" / k_seed / "sub{subsplit}.prob.npy",
    output:
        multiext(str(res / "masked" / k_seed / "sub{subsplit}_posterior"), 
                    ".npy", ".missingness", ".names", ".labels")
    params:
        outbase = lambda wc, output: output[0][:-4]
    shell:
        "{P} {MASK_ANCESTRY} --posterior {input.masker} "
        "--loglike {input.l} --labels {POP_LABELS}  --names {SAMPLES} "
        "--out {params.outbase} --max_missing {MAX_MISSINGNESS} {HET_HOM}"

rule mask_decoding:
    input:
        l = res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
        masker = res / "fatass" / k_seed / "sub{subsplit}.path",
    output:
        multiext(str(res / "masked" / k_seed / "sub{subsplit}_decoding"), 
                    ".npy", ".missingness", ".names", ".labels")
    params:
        outbase = lambda wc, output: output[0][:-4]
    shell:
        "{P} {MASK_ANCESTRY} --decoding {input.masker} "
        "--loglike {input.l} --labels {POP_LABELS}  --names {SAMPLES} "
        "--out {params.outbase} --max_missing {MAX_MISSINGNESS} {HET_HOM}"

rule pca:
    input:
        masker = res / "masked" / k_seed / "sub{subsplit}_{masktype}.npy",
    output:
        res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}.eigenvecs",
        res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}.eigenvals",
    params:
        outbase = lambda wc, output: output[0][:-10]
    threads: 60
    shell: """
    {HAPLONET} pca -l {input} -t {threads} --iterative {PC_ITERATIVE} \
        -o {params.outbase} --filter {wildcards.maf_filter}
    """

rule plot:
    input:
        res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}.eigenvecs",
        res / "masked" / k_seed / "sub{subsplit}_{masktype}.labels",
    output:
        res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}_{pc1}_{pc2}.png",
    shell:
        "Rscript {PLOTTER} {input} {wildcards.pc1} {wildcards.pc2} {output}"
