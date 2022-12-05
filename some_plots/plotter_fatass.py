import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

d = np.load("sub0.prob.npy")
d = np.exp(d)
names = np.loadtxt("na_inclafr.samples", dtype=str)
# n_hap1 = [x+"_hap1" for x in names]
# n_hap2 = [x+"_hap2" for x in names]
# names_dup = np.array(list(zip(n_hap1, n_hap2))).flatten()
names_dup = np.repeat(names, 2)
haps = np.array([["hap1", "hap2"] for _ in names]).flatten()

K = d.shape[2]
assert d.shape[1] == names_dup.shape[0]

idx_d = np.array(["HG02006" in x for x in names_dup])
haps_test = haps[idx_d]
d_test = d[:, idx_d, :]

N_obs = 525


fig, ax = plt.subplots(2,1, figsize=(24,8), sharex=True)
for a in ax:
    a.margins(x=0)
    a.axhline(y=0.95, color='black', alpha=0.5)
ax[0].set_title("hap1")
for lab in range(K):
    ax[0].plot(range(N_obs), d_test[:N_obs, 0, lab], label=f"K{lab}")

ax[1].set_title("hap2")
for lab in range(K):
    ax[1].plot(range(N_obs), d_test[:N_obs, 1, lab], label=f"K{lab}")
ax[0].legend()
fig.savefig("test3.png", bbox_inches='tight', dpi=500)

