import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#data[["counts"]].plot(kind='hist', xlabel="N worst fit per individual", bins=20).get_figure().savefig("test5.png")
d = np.load("na.new.loglike.subsplit.allchrom.npy")
test = d.max(axis=2).argmin(axis=1)
cts = np.unique(test, return_counts=True)
labels = np.loadtxt("na.new.samples", dtype=str)
pops = np.loadtxt("na.new.pop", dtype=str)
data = pd.DataFrame(cts).T
data.columns = ["index", "counts"]
data["samples"] = labels[data[["index"]] // 2]
data["pops"] = pops[data[["index"]] // 2]

fig, ax = plt.subplots(1,1, figsize=(8,8))
g = sns.histplot(data[["counts"]], ax=ax)
ax.set_xlabel("Worst fit per haplotype across all windows")
ax.set_title(f"{data.counts.sum()} Windows")
for idx, (index, row) in enumerate(data.query("counts>20").iterrows()):
    ax.annotate(f"{row.samples}_{row.pops}_{row.counts}", (row.counts-2,40 + idx*5))
fig.savefig("worst_haplotype_per_window.png")

d_max = d.max(axis=2)
d_max_norm = (d_max - d_max.mean(axis=1, keepdims=True)) / d_max.std(axis=1, keepdims=True)
to_plot = d_max_norm.min(axis=1)

fig, ax = plt.subplots(1,3, figsize=(24,8))
g = sns.histplot(to_plot, ax=ax[0])
g.set_xlabel("SD")
g.set_ylabel("Count")
g.set_title("Worst SD per window")
g = sns.histplot(d_max_norm.flatten(), ax=ax[1])
g.set_ylabel("Count")
g.set_xlabel("SD")
g.set_title(f"W:{d_max_norm.shape[0]}, Haplotypes:{d_max_norm.shape[1]}, Nobs: {np.product(d_max_norm.shape)}")
g = sns.histplot(d_max_norm.flatten(), ax=ax[2])
g.set_ylabel("Count (log)")
g.set_xlabel("SD")
g.set_yscale("log")
g.set_title(f"W:{d_max_norm.shape[0]}, Haplotypes:{d_max_norm.shape[1]}, Nobs: {np.product(d_max_norm.shape)}")
fig.savefig("all_windows.png", bbox_inches='tight')

## careful.
d[d_max_norm< - d_max_norm.max(),:] = np.log(1.0/d.shape[2])
np.save("na.new.loglike.subsplit.allchrom_masked.npy", d)
# def norm(x):
#     return (x - x.mean()) / x.std()

# def pick_low_sd(x):
#     x_norm = norm(x)
#     return x_norm[x_norm==x_norm.min()]