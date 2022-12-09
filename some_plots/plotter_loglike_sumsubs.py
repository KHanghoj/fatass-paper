import numpy as np
import matplotlib.pyplot as plt

d = [np.load(f"../res_inclafr/haplonet_split/allchrom.sub{x}.loglike.npy") for x in [0,2,4]]
pos = [np.loadtxt(f"../res_inclafr/basepos/allchrom.sub{x}.positions.txt", usecols=1) for x in [0,2,4]]
for idx, x in enumerate(d):
    print(idx, x.shape)
names = np.loadtxt("../configs/na_inclafr.samples", dtype=str)
names_dup = np.repeat(names, 2)
maxs = [x.max(axis=2) for x in d]
for name in ("HG00096", "HG02089", "HG02006"):

    K = d[0].shape[2]
    assert d[0].shape[1] == names_dup.shape[0]
    idx_d = np.array([name in x for x in names_dup])
    hap1 = np.where(idx_d)[0][0]
    ## sum logs
    ys = np.zeros((d[0].shape[0], len(d)), dtype=float)
    print(ys.shape)
    ys[:, 0] = maxs[0][:, hap1]
    for x in range(2):
        ys[:,1] += maxs[1][x::2, hap1]
    for x in range(4):
        ys[:,2] += maxs[2][x::4, hap1]

    # chrom1
    N_obss = 525
    megamin = ys[:N_obss].min()
    fig, ax = plt.subplots(1,1, figsize=(24,8), sharex=True)
    fig.suptitle(f"{name} hap1")
    ax.margins(x=0)
    for idx, sub in enumerate([0,2,4]):
        ax.plot(pos[0][:N_obss], ys[:N_obss, idx], label=f"lab{sub}")
    ax.set_ylim([megamin, 0])
    ax.legend()
    ax.grid()
    fig.tight_layout()
    fig.savefig(f"chr1_{name}_loglikes.png", bbox_inches='tight', dpi=500)
