import numpy as np
import matplotlib.pyplot as plt

for name in ("HG00096", "HG02089", "HG02006"):
    d = [np.exp(np.load(f"../res_inclafr/fatass/K3_s1/sub{x}.prob.npy")) for x in [0,2,4]]
    pos = [np.loadtxt(f"../res_inclafr/basepos/allchrom.sub{x}.positions.txt", usecols=1) for x in [0,2,4]]
    for idx, x in enumerate(d):
        print(idx, x.shape)
    names = np.loadtxt("../configs/na_inclafr.samples", dtype=str)
    names_dup = np.repeat(names, 2)

    K = d[0].shape[2]
    assert d[0].shape[1] == names_dup.shape[0]
    idx_d = np.array([name in x for x in names_dup])
    d_test = [x[:, idx_d, :] for x in d]
    print(np.where(idx_d)[0][0])

    N_obss = [525*x for x in [1,2,4]]
    fig, axes = plt.subplots(3,1, figsize=(24,8), sharex=True)
    fig.suptitle(f"{name} hap1")
    for a in axes:
        a.margins(x=0)
        a.axhline(y=0.95, color='black', alpha=0.5)

    for idx, sub in enumerate([0,2,4]):
        ax = axes[idx]
        ax.set_title(f"sub{sub}")
        for lab in range(K):
            ax.scatter(pos[idx][:N_obss[idx]], d_test[idx][:N_obss[idx], 0, lab], s=3)
            ax.step(pos[idx][:N_obss[idx]], d_test[idx][:N_obss[idx], 0, lab], where="post", label=f"q{lab}")
            # ax.plot(range(N_obss[idx]), d_test[idx][:N_obss[idx], 0, lab], label=f"q{lab}")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(f"chr1_{name}_across_sub.png", bbox_inches='tight', dpi=500)

    N_obss = [100*x for x in [1,2,4]]
    fig, axes = plt.subplots(3,1, figsize=(24,8), sharex=True)
    fig.suptitle(f"{name} hap1")
    for a in axes:
        a.margins(x=0)
        a.axhline(y=0.95, color='black', alpha=0.5)

    for idx, sub in enumerate([0,2,4]):
        ax = axes[idx]
        ax.set_title(f"sub{sub}")
        for lab in range(K):
            ax.scatter(pos[idx][:N_obss[idx]], d_test[idx][:N_obss[idx], 0, lab], s=3)
            ax.step(pos[idx][:N_obss[idx]], d_test[idx][:N_obss[idx], 0, lab], where="post", label=f"q{lab}")
            # ax.plot(range(N_obss[idx]), d_test[idx][:N_obss[idx], 0, lab], label=f"q{lab}")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(f"chr1_100_{name}_across_sub.png", bbox_inches='tight', dpi=500)


