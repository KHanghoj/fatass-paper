import numpy as np
import matplotlib.pyplot as plt

for name in ("HG00096", "HG02089", "HG02006"):
    d = np.exp(np.load(f"../res_inclafr/fatass/K3_s1/cool.prob.npy"))
    pos = np.loadtxt(f"../res_inclafr/basepos/allchrom.sub4.positions.txt", usecols=[1,2]).mean(1)
    names = np.loadtxt("../configs/na_inclafr.samples", dtype=str)
    names_dup = np.repeat(names, 2)

    K = d.shape[2]
    assert d.shape[1] == names_dup.shape[0]
    idx_d = np.array([name in x for x in names_dup])
    d_test = d[:,idx_d,:]
    print(np.where(idx_d)[0][0])

    N_obss = 525*4
    fig, ax = plt.subplots(1,1, figsize=(24,3), sharex=True)
    fig.suptitle(f"{name} hap1")
    ax.margins(x=0)
    ax.axhline(y=0.95, color='black', alpha=0.5)
    for lab in range(K):
        ax.scatter(pos[:N_obss], d_test[:N_obss, 0, lab], s=3)
        ax.step(pos[:N_obss], d_test[:N_obss, 0, lab], where="post", label=f"q{lab}")
        # ax.plot(range(N_obss[idx]), d_test[idx][:N_obss[idx], 0, lab], label=f"q{lab}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"chr1_{name}_across_sub_temp.png", bbox_inches='tight', dpi=500)

