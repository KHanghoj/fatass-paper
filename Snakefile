import itertools as it
from pathlib import Path
def add_ext(p, *args):
    base = str(p)
    for x in args:
        base += str(x)
    return Path(base)
    # return Path(str(p)+str(e))
# bcftools view -S na.samples -q 0.05 -Q 0.95 -m 2 -M 2 -v snps -Oz -o vcf/na.chr5.vcf.gz --threads 4 /emc/kristian/1kg_20220422/1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gza
SCRIPT_DIR = os.path.join(workflow.basedir, "scripts")
P = config["python"]
VCF = config["vcf"]
HAPLONET = config["haplonet"]
# MASK_ANCESTRY = config["mask_anc"]
# PLOTTER = config["plotter"]
BCFTOOLS = config["bcftools"]
MAF_FILTER = config["filter"]
WINSIZE_SUB0 = config["winsize_sub0"]
LOGLIKE_DIFF = config["loglike_diff"]
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

pcs = ["PC1", "PC2", "PC3", "PC4", "PC5"]
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
            ),
        expand(res / "basepos" / "allchrom.sub{subsplit}.positions.txt",
            subsplit = SUBSPLIT_LST
        )

rule pca_decoding:
    input:
        expand(res / "pca" / k_seed / "MAF{maf_filter}" / 
            "sub{subsplit}_decoding_{pc_comb[0]}_{pc_comb[1]}.png",
            maf_filter = MAF_FILTER,
            subsplit = SUBSPLIT_LST,
            pc_comb = pcs_comb
            ),
        expand(res / "basepos" / "allchrom.sub{subsplit}.positions.txt",
            subsplit = SUBSPLIT_LST
        )

rule pca_posterior:
    input:
        expand(res / "pca" / k_seed / "MAF{maf_filter}" / 
            "posterior_{pc_comb[0]}_{pc_comb[1]}.png",
            maf_filter = MAF_FILTER,
            subsplit = SUBSPLIT_LST,
            pc_comb = pcs_comb
            ),
        res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt.counts.hist.png",

rule all_admixture:
    input:
        expand(res / "admix" / k_seed / "sub{subsplit}.q", subsplit=SUBSPLIT_LST),
        expand(res / "admix" / k_seed / "sub{subsplit}.q.png", subsplit=SUBSPLIT_LST),

# rule pca_posterior:
#     input:
#         expand(res / "pca" / k_seed / "MAF{maf_filter}" / 
#             "sub{subsplit}_posterior_{pc_comb[0]}_{pc_comb[1]}.png",
#             maf_filter = MAF_FILTER,
#             subsplit = SUBSPLIT_LST,
#             pc_comb = pcs_comb
#             ),
#         expand(res / "basepos" / "allchrom.sub{subsplit}.positions.txt",
#             subsplit = SUBSPLIT_LST
#         )
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

rule query_pos:
    input:
        res / "vcf" / "{chrom}.vcf.gz"
    output:
        res / "vcf" / "{chrom}.positions",
    shell:
        "{BCFTOOLS} query -f '%POS\\n' {input} > {output}"

rule get_bp:
    input:
        res / "vcf" / "{chrom}.positions",
    output:
        res / "basepos" / "sub{subsplit}" / "{chrom}.positions.txt",
    params:
        subsplit = lambda wc: wc.subsplit if int(wc.subsplit)>0 else 1, 
        windowsize_nosub = WINSIZE_SUB0
    run:
        import sys
        import numpy as np
        from math import ceil

        pos = input[0]
        chrom = wildcards.chrom
        subsplit = int(params.subsplit)
        windowsize_nosub = params.windowsize_nosub
        out = output[0]
        windowsize_sub = windowsize_nosub // subsplit
        S = np.loadtxt(pos, dtype=int)
        if (S.shape[0] % windowsize_nosub) < windowsize_nosub//2:
        	nSeg = S.shape[0]//windowsize_nosub
        else:
        	nSeg = ceil(S.shape[0]/windowsize_nosub)

        M = np.empty((nSeg*subsplit, 3), dtype=object)

        M[:,0] = chrom
        for i in range(nSeg):
        	positions = S[i*windowsize_nosub:i*windowsize_nosub+windowsize_nosub]
        	splits = np.array_split(positions, subsplit)
        	for subw, w in enumerate(splits):
        		M_idx = i*subsplit + subw
        		if i == 0 and subw==0:
        			M[M_idx, 1]=w[0]
        		else:
        			M[M_idx-1, 2] = w[0]
        			M[M_idx, 1] = w[0]
        M[-1, 2] = S[-1]
        np.savetxt(out, M, delimiter="\t", fmt="%s")

# rule get_bp:
#     input:
#         res / "vcf" / "{chrom}.vcf.gz"
#     output:
#         res / "basepos" / "sub0" / "{chrom}.positions.txt",
#     params:
#         bp = 1024,
#         out = lambda wc, output: output[0][:-14]
#     shell:
#         "{HAPLONET} convert -l {params.bp} -w -o {params.out} --vcf {input}"

# rule get_bp_subsplit:
#     input:
#         vcf = res / "vcf" / "{chrom}.vcf.gz",
#         nosub = res / "basepos" / "sub0" / "{chrom}.positions.txt",
#     output:
#         t = temp(res / "basepos" / "sub{subsplit}" / "{chrom}.temp.positions.txt"),
#         fixed = res / "basepos" / "sub{subsplit}" / "{chrom}.positions.txt",
#     params:
#         bp = lambda wc: 1024 // int(wc.subsplit),
#         out = lambda wc, output: output.t[:-14]
#     wildcard_constraints:
#         subsplit = "|".join([str(x) for x in SUBSPLIT_LST if int (x) != 0])
#     run:
#         shell("{HAPLONET} convert -l {params.bp} -w -o {params.out} --vcf {input.vcf}")
#         with open(input.nosub, 'r') as fh:
#             nlines = 0
#             for line in fh:
#                 nlines+=1
#             maxlines = nlines * int(wildcards.subsplit)
#             maxval = line.rstrip().split()[-1]    
#         with open(output.t, 'r') as fhin, open(output.fixed, 'w') as fhout:
#             for idx, line in enumerate(fhin, start=1):
#                 chrom, start, end = line.rstrip().split()
#                 if idx == maxlines:
#                     print(f"{chrom}\t{start}\t{maxval}", file=fhout)
#                     break 
#                 else:
#                     print(f"{chrom}\t{start}\t{end}", file=fhout)

rule concat_chrom_bp:
    input:
        expand(res / "basepos" / "sub{subsplit}" / "{chrom}.positions.txt",
            chrom=CHROMLIST, allow_missing=True)
    output:
        res / "basepos" / "allchrom.sub{subsplit}.positions.txt"
    shell:
        "cat {input} > {output}"
        
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
        windowsize_nosub = WINSIZE_SUB0
    shell:
        "{HAPLONET} train -x {params.windowsize_nosub} --vcf {input} "
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
         "{P} {SCRIPT_DIR}/concat_npyV2.py {output} {input} > {log}"
        
rule admix:
    input:
        res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy" 
    output:
        res / "admix" / k_seed / "sub{subsplit}.f.npy",
        res / "admix" / k_seed / "sub{subsplit}.q"
    log:
        res / "admix" / k_seed / "sub{subsplit}.log",
    params:
        out = lambda wc, output: output[0][:-6]
    threads: 20
    shell:
        "{HAPLONET} admix --like {input} "
        "--K {K} --seed {SEED} --out {params.out} "
        "--threads {threads} > {log}"

rule plot_admix:
    input:
        res / "admix" / k_seed / "sub{subsplit}.q"
    output:
        res / "admix" / k_seed / "sub{subsplit}.q.png"
    shell:
        "Rscript {SCRIPT_DIR}/plot_admixture.R {POP_LABELS} {input} {output}"

rule get_cool_windows:
    input:
        l = expand(res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
                    subsplit = SUBSPLIT_LST),
        w = res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.positions.txt",
    output:
        w = res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt",
        c = res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt.counts",
        b = res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt.bool",
    shell:
        "{P} {SCRIPT_DIR}/prep_windows_fatassV2.py -w {input.w} -L {input.l} -o {output.w} --loglike_diff {LOGLIKE_DIFF}" 

rule plot_win_dist:
    input:
        res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt.counts",
    output:
        res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt.counts.hist.png",
        res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt.counts.perpop.png",
    shell:
        "Rscript {SCRIPT_DIR}/plot_windowsubsplit.R {input} {POP_LABELS} {output}"

def set_alpha():
    if str(ALPHA).lower() == "est":
        return "--alpha_save --optim"
    else:
        return f"--alpha_save --alpha {ALPHA}"

rule fatass:
    input:
        l = res / "haplonet_split" / f"allchrom.sub{SUBSPLIT_MAX}.loglike.npy", 
        f = res / "admix" / k_seed / f"sub{SUBSPLIT_MAX}.f.npy",
        q = res / "admix" / k_seed / f"sub{SUBSPLIT_MAX}.q",
        w = res / "basepos" / f"allchrom.sub{SUBSPLIT_MAX}.cool.positions.txt",
    output:
        res / "fatass" / k_seed / f"cool.prob.npy",
        res / "fatass" / k_seed / f"cool.path",
        res / "fatass" / k_seed / f"cool.alpha",
        res / "fatass" / k_seed / f"cool.windows",
    params:
        out = lambda wc, output: output[0][:-9],
        alpha = set_alpha()
    threads: 10
    shell:
        "{HAPLONET} fatash {params.alpha} "
        "-w {input.w} "
        "--window_save "
        "--like {input.l} "
        "--prop {input.q} --freq {input.f} "
        "--out {params.out} "
        "--threads {threads}"

rule mask_posterior:
    input:
        l = res / "haplonet_split" / f"allchrom.sub{SUBSPLIT_MAX}.loglike.npy", 
        masker = res / "fatass" / k_seed / "cool.prob.npy",
    output:
        multiext(str(res / "masked" / k_seed / "posterior"), 
                    ".npy", ".missingness", ".names", ".labels")
    params:
        outbase = lambda wc, output: output[0][:-4]
    shell:
        "{P} {SCRIPT_DIR}/mask_ancV3.py --posterior {input.masker} "
        "--loglike {input.l} --labels {POP_LABELS}  --names {SAMPLES} "
        "--out {params.outbase} --max_missing {MAX_MISSINGNESS} {HET_HOM}"

rule mask_decoding:
    input:
        l = res / "haplonet_split" / f"allchrom.sub{SUBSPLIT_MAX}.loglike.npy", 
        masker = res / "fatass" / k_seed / "cool.path",
    output:
        multiext(str(res / "masked" / k_seed / "decoding"), 
                    ".npy", ".missingness", ".names", ".labels")
    params:
        outbase = lambda wc, output: output[0][:-4]
    shell:
        "{P} {SCRIPT_DIR}/mask_ancV3.py --decoding {input.masker} "
        "--loglike {input.l} --labels {POP_LABELS}  --names {SAMPLES} "
        "--out {params.outbase} --max_missing {MAX_MISSINGNESS} {HET_HOM}"

rule pca:
    input:
        masker = res / "masked" / k_seed / "{masktype}.npy",
    output:
        res / "pca" / k_seed / "MAF{maf_filter}" / "{masktype}.eigenvecs",
        res / "pca" / k_seed / "MAF{maf_filter}" / "{masktype}.eigenvals",
    log:
        res / "pca" / k_seed / "MAF{maf_filter}" / "{masktype}.log",
    params:
        outbase = lambda wc, output: output[0][:-10]
    threads: 28
    shell: """
    {HAPLONET} pca -l {input} -t {threads} --iterative {PC_ITERATIVE} \
        -o {params.outbase} --iterations 1000 --filter {wildcards.maf_filter} > {log}
    """

rule plot:
    input:
        res / "pca" / k_seed / "MAF{maf_filter}" / "{masktype}.eigenvecs",
        res / "masked" / k_seed / "{masktype}.labels",
        res / "masked" / k_seed / "{masktype}.names",
    output:
        res / "pca" / k_seed / "MAF{maf_filter}" / "{masktype}_{pc1}_{pc2}.png",
    shell:
        "Rscript {SCRIPT_DIR}/plot.R {input} {wildcards.pc1} {wildcards.pc2} {output}"

# rule fatass:
#     input:
#         l = res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
#         f = res / "admix" / k_seed / "sub{subsplit}.f.npy",
#         q = res / "admix" / k_seed / "sub{subsplit}.q",
#         w = res / "basepos" / "allchrom.sub{subsplit}.positions.txt",
#     output:
#         res / "fatass" / k_seed / "sub{subsplit}.prob.npy",
#         res / "fatass" / k_seed / "sub{subsplit}.path",
#         res / "fatass" / k_seed / "sub{subsplit}.alpha",
#         res / "fatass" / k_seed / "sub{subsplit}.windows",
#     params:
#         out = lambda wc, output: output[0][:-9],
#         alpha = set_alpha()
#     threads: 10
#     shell:
#         "{HAPLONET} fatash {params.alpha} "
#         "-w <( cut -f 2 {input.w} ) "
#         "--window_save "
#         "--like {input.l} "
#         "--prop {input.q} --freq {input.f} "
#         "--out {params.out} "
#         "--threads {threads}"

# rule mask_posterior:
#     input:
#         l = res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
#         masker = res / "fatass" / k_seed / "sub{subsplit}.prob.npy",
#     output:
#         multiext(str(res / "masked" / k_seed / "sub{subsplit}_posterior"), 
#                     ".npy", ".missingness", ".names", ".labels")
#     params:
#         outbase = lambda wc, output: output[0][:-4]
#     shell:
#         "{P} {SCRIPT_DIR}/mask_ancV3.py --posterior {input.masker} "
#         "--loglike {input.l} --labels {POP_LABELS}  --names {SAMPLES} "
#         "--out {params.outbase} --max_missing {MAX_MISSINGNESS} {HET_HOM}"

# rule mask_decoding:
#     input:
#         l = res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
#         masker = res / "fatass" / k_seed / "sub{subsplit}.path",
#     output:
#         multiext(str(res / "masked" / k_seed / "sub{subsplit}_decoding"), 
#                     ".npy", ".missingness", ".names", ".labels")
#     params:
#         outbase = lambda wc, output: output[0][:-4]
#     shell:
#         "{P} {SCRIPT_DIR}/mask_ancV3.py --decoding {input.masker} "
#         "--loglike {input.l} --labels {POP_LABELS}  --names {SAMPLES} "
#         "--out {params.outbase} --max_missing {MAX_MISSINGNESS} {HET_HOM}"

# rule pca:
#     input:
#         masker = res / "masked" / k_seed / "sub{subsplit}_{masktype}.npy",
#     output:
#         res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}.eigenvecs",
#         res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}.eigenvals",
#     params:
#         outbase = lambda wc, output: output[0][:-10]
#     threads: 40
#     shell: """
#     {HAPLONET} pca -l {input} -t {threads} --iterative {PC_ITERATIVE} \
#         -o {params.outbase} --filter {wildcards.maf_filter}
#     """

# rule plot:
#     input:
#         res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}.eigenvecs",
#         res / "masked" / k_seed / "sub{subsplit}_{masktype}.labels",
#     output:
#         res / "pca" / k_seed / "MAF{maf_filter}" / "sub{subsplit}_{masktype}_{pc1}_{pc2}.png",
#     shell:
#         "Rscript {SCRIPT_DIR}/plot.R {input} {wildcards.pc1} {wildcards.pc2} {output}"
