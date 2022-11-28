import itertools as it
from pathlib import Path
def add_ext(p, *args):
    base = str(p)
    for x in args:
        base += str(x)
    return Path(base)
    # return Path(str(p)+str(e))

HAPLONET = config["haplonet"]
P = config["python"]
MASK_ANCESTRY = config["mask_anc"]
PLOTTER = config["plotter"]

MAF_FILTER = config["filter"]
PC_ITERATIVE = config["iterative"]
MAX_MISSINGNESS = config["max_missingness"]
HET_HOM = "--het_mask" if config["het_hom"] == "het" else ""
SAMPLES = config["samples"]
POP_LABELS = config["labels"]
POPULATION = config["population"]
# LOGLIKE = "na.loglike{sp}allchrom.npy"
LOGLIKE = "na.new.loglike{sp}allchrom.npy"
LOGLIKE = "na.new.loglike{sp}allchrom_masked.npy"

res = Path(config["output"])

# bname = f"{POPULATION}{{sp}}"
bname = f"{POPULATION}.new{{sp}}"
fs = bname + "K2.fatass"
# fs = bname + "K2.s1.fatass"
mk = bname + "allchrom.{t}"

with open(config["chromlist"], 'r') as fh:
    CHROMLIST = [x.rstrip() for x in fh]

pcs = ["PC1", "PC2", "PC3"]
pcs_comb = list(it.combinations(pcs, 2))
sp_wc  = [".subsplit.", "."]
t_wc = ["prob.npy", "path"]


wildcard_constraints:
    sp = "|".join(sp_wc),
    t = "|".join(t_wc),
    chrom = "|".join(CHROMLIST),
    maf_filter = "|".join(map(str, MAF_FILTER))

def get_concatter(t):
    if t == "prob.npy":
        return "scripts/concat_npy.py"    
    else:
        return "scripts/concat_paths.py" 

rule all:
    input:
        expand(res / "pca" / "MAF{maf_filter}" / add_ext(mk, ".{pc_comb[0]}.{pc_comb[1]}", ".png"), 
            maf_filter = MAF_FILTER,
            sp = sp_wc,
            t = t_wc,
            pc_comb = pcs_comb
            )

rule gen_filelist:
    input:
        expand("fatass" / add_ext(fs, ".{chrom}.{t}"), chrom=CHROMLIST, allow_missing=True)
    output:
        res / "filelist" / add_ext(fs, ".{t}.filelist")
    run:
        with open(output[0], 'w') as fh:
            for line in input:
                print(line, file=fh)
    
rule concat_post_decode:
    input:
        res / "filelist" / add_ext(fs, ".{t}.filelist")
    output:
        res / "concat" / add_ext(fs, ".allchrom.{t}")
    params:
        software = lambda wc: get_concatter(wc.t)
    shell:
        "{P} {params.software} {input} {output}"

rule mask:
    input:
        loglike = LOGLIKE,
        masker = res / "concat" / add_ext(fs, ".allchrom.{t}")
    output:
        multiext(str(res / "masked" / mk), 
                    ".npy", ".missingness", ".names", ".labels")
    params:
        dt = lambda wc: "--posterior" if wc.t == "prob.npy" else "--decoding",
        outbase = lambda wc, output: output[0][:-4]
    shell: """
    {P} {MASK_ANCESTRY} {params.dt} {input.masker} --loglike {input.loglike} \
        --labels {POP_LABELS}  --names {SAMPLES}  --out {params.outbase} \
        --max_missing {MAX_MISSINGNESS} {HET_HOM}
    """

rule pca:
    input:
        res / "masked" / add_ext(mk, ".npy")
    output:
        res / "pca" / "MAF{maf_filter}" / add_ext(mk, ".eigenvecs"),
        res / "pca" / "MAF{maf_filter}" / add_ext(mk, ".eigenvals"),
    params:
        outbase = lambda wc, output: output[0][:-10]
    threads: 60
    shell: """
    {HAPLONET} pca -l {input} -t {threads} --iterative {PC_ITERATIVE} \
        -o {params.outbase} --filter {wildcards.maf_filter}
    """

rule plot:
    input:
        res / "pca" / "MAF{maf_filter}" / add_ext(mk, ".eigenvecs"),
        res / "masked" / add_ext(mk, ".labels")
    output:
        res / "pca" / "MAF{maf_filter}" / add_ext(mk, ".{pc1}.{pc2}", ".png"),
    shell:
        "Rscript {PLOTTER} {input} {wildcards.pc1} {wildcards.pc2} {output}"
    