import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

d = np.load("sub2_posterior.npy")
names = np.loadtxt("sub2_posterior.names", dtype=str)
# n_hap1 = [x+"_hap1" for x in names]
# n_hap2 = [x+"_hap2" for x in names]
# names_dup = np.array(list(zip(n_hap1, n_hap2))).flatten()
names_dup = np.repeat(names, 2)
haps = np.array([["hap1", "hap2"] for _ in names]).flatten()

assert d.shape[1] == names_dup.shape[0]

idx_d = np.array(["HG02006" in x for x in names_dup])
haps_test = haps[idx_d]
d_test = d[:, idx_d, :]

n_k = [(x[:-2], x[-2:]) for x in names_dup[idx_d]]
column_names = [(hap,k) for (_,k),hap in zip(n_k, haps_test)]
C = d.shape[2]
missing_val = np.log(1.0/C)
data = pd.DataFrame(~np.all(np.isclose(d_test, missing_val),axis=2)+0)
data.columns = pd.MultiIndex.from_tuples(column_names)

hap1 = data.loc[:2000, "hap1"]
hap2 = data.loc[:2000, "hap2"]


fig, ax = plt.subplots(2,1, figsize=(24,8))
ax[0].set_title("hap1")
for lab in hap1.columns:
    ax[0].plot(hap1.index, hap1[lab], label=lab)
ax[1].set_title("hap2")
for lab in hap2.columns:
    ax[1].plot(hap1.index, hap2[lab], label=lab)
ax[0].legend()
fig.savefig("test3.png", bbox_inches='tight', dpi=500)

