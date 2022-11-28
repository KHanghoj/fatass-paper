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
        expand(res / "fatass" / "sub{subsplit}" / "allchrom.prob.npy", subsplit=SUBSPLIT_LST)

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
        res / "haplonet" / "sub{subsplit}" / "{chrom}.loglike.npy"  
    threads: 20
    params:
        out = lambda wc, output: output[0][:-4],
        subsplit = lambda wc: "" if wc.subsplit == 0 else f"--subsplit {wc.subsplit}" 
    shell:
        "{HAPLONET} train --vcf {input} "
        "--out {params.out} {params.subsplit} "
        "--threads {threads}"

rule concat:
    input:
        expand(res / "haplonet" / "sub{subsplit}" / "{chrom}.loglike.npy" ,
            chrom = CHROMLIST, allow_missing=True)
    output:
        res / "haplonet" / "sub{subsplit}" / "allchrom.loglike.npy" 
    shell:
         "{P} ./scripts/concat_npyV2.py {output} {input}"
        
rule admix:
    input:
        res / "haplonet" / "sub{subsplit}" / "allchrom.loglike.npy" 
    output:
        res / "admix" / "sub{subsplit}" / f"K{K}.s{SEED}.f.npy",
        res / "admix" / "sub{subsplit}" / f"K{K}.s{SEED}.q"
    params:
        out = lambda wc, output: output[0][:-6]
    threads: 20
    shell:
        "{HAPLONET} admix --like {input} "
        "--K {K} --seed {SEED} --out {params.out} "
        "--threads {threads}"

rule fatass:
    input:
        l = res / "haplonet" / "sub{subsplit}" / "allchrom.loglike.npy",
        f = res / "admix" / "sub{subsplit}" / f"K{K}.s{SEED}.f.npy",
        q = res / "admix" / "sub{subsplit}" / f"K{K}.s{SEED}.q",
    output:
        res / "fatass" / "sub{subsplit}" / "allchrom.prob.npy",
        res / "fatass" / "sub{subsplit}" / "allchrom.path",
    params:
        out = lambda wc, output: output[0][:-9]
    threads: 10
    shell:
        "{HAPLONET} fatash --like {input.l} "
        "--prop {input.q} --freq {input.f} "
        "--out {params.out} "
        "--threads {threads}"
