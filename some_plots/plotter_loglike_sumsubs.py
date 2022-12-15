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
    pos_sub = pos[0][:N_obss]
    ys_sub = ys[:N_obss]
    megamin = ys_sub.min()
    fig, ax = plt.subplots(1,1, figsize=(24,8), sharex=True)
    fig.suptitle(f"{name} hap1")
    ax.margins(x=0)
    for idx, sub in enumerate([0,2,4]):
        ax.plot(pos_sub, ys_sub[:, idx], label=f"lab{sub}")
    ax.set_ylim([megamin, 0])
    ax.legend()
    ax.grid()
    fig.tight_layout()
    fig.savefig(f"chr1_{name}_loglikes.png", bbox_inches='tight', dpi=500)
    plt.close(fig)

    ll_ratio_02 = -2 * (ys_sub[:, 0] - ys_sub[:,1])
    ll_ratio_04 = -2 * (ys_sub[:, 0] - ys_sub[:,2])

    fig, ax = plt.subplots(1,1, figsize=(24,8), sharex=True)
    fig.suptitle(f"{name} hap1")
    ax.margins(x=0)
    ax.plot(pos_sub, ll_ratio_02, label="lab02")
    ax.plot(pos_sub, ll_ratio_04, label="lab04")    
    ax.legend()
    ax.grid()
    fig.tight_layout()
    fig.savefig(f"chr1_{name}_llr.png", bbox_inches='tight', dpi=500)
    plt.close(fig)
    

    ys_sub_norm = ys_sub[:, 1:] - ys_sub[:,0].reshape(N_obss,-1)
    fig, ax = plt.subplots(1,1, figsize=(24,8), sharex=True)
    fig.suptitle(f"{name} hap1")
    ax.margins(x=0)
    ax.plot(pos_sub, ys_sub_norm[:,0], label="lab02")
    ax.plot(pos_sub, ys_sub_norm[:,1], label="lab04")    
    ax.legend()
    ax.grid()
    fig.tight_layout()
    fig.savefig(f"chr1_{name}_lldiff.png", bbox_inches='tight', dpi=500)
    plt.close(fig)
