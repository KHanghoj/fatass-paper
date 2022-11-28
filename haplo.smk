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
BCFTOOLS = config["bcftools"]
MAF_FILTER = config["filter"]
K = config["K"]
SEED = config["seed"]
SUBSPLIT_LST = config["subsplit"]
SUBSPLIT_MAX = max(SUBSPLIT_LST)
SUBSPLIT_ARG = "" if SUBSPLIT_MAX == 0 else f"--subsplit {SUBSPLIT_MAX}" 
SAMPLES = config["samples"]
# PC_ITERATIVE = config["iterative"]
# MAX_MISSINGNESS = config["max_missingness"]
# HET_HOM = "--het_mask" if config["het_hom"] == "het" else ""
# POP_LABELS = config["labels"]
# POPULATION = config["population"]
# # LOGLIKE = "na.loglike{sp}allchrom.npy"
# LOGLIKE = "na.new.loglike{sp}allchrom.npy"
# LOGLIKE = "na.new.loglike{sp}allchrom_masked.npy"

res = Path(config["output"])

# # bname = f"{POPULATION}{{sp}}"
# bname = f"{POPULATION}.new{{sp}}"
# fs = bname + "K2.fatass"
# # fs = bname + "K2.s1.fatass"
# mk = bname + "allchrom.{t}"

with open(config["chromlist"], 'r') as fh:
    CHROMLIST = [x.rstrip() for x in fh]

rule all:
    input:
        expand(res / "fatass" / f"K{K}_s{SEED}" / "sub{subsplit}.prob.npy", subsplit=SUBSPLIT_LST)

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
    threads: 10
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
        res / "admix" / f"K{K}_s{SEED}" / "sub{subsplit}.f.npy",
        res / "admix" / f"K{K}_s{SEED}" / "sub{subsplit}.q"
    params:
        out = lambda wc, output: output[0][:-6]
    threads: 20
    shell:
        "{HAPLONET} admix --like {input} "
        "--K {K} --seed {SEED} --out {params.out} "
        "--threads {threads}"

rule fatass:
    input:
        l = res / "haplonet_split" / "allchrom.sub{subsplit}.loglike.npy", 
        f = res / "admix" / f"K{K}_s{SEED}" / "sub{subsplit}.f.npy",
        q = res / "admix" / f"K{K}_s{SEED}" / "sub{subsplit}.q"
    output:
        res / "fatass" / f"K{K}_s{SEED}" / "sub{subsplit}.prob.npy",
        res / "fatass" / f"K{K}_s{SEED}" / "sub{subsplit}.path",
    params:
        out = lambda wc, output: output[0][:-9]
    threads: 10
    shell:
        "{HAPLONET} fatash --like {input.l} "
        "--prop {input.q} --freq {input.f} "
        "--out {params.out} "
        "--threads {threads}"
