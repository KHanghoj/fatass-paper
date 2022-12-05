import numpy as np
import matplotlib.pyplot as plt

d = [np.load(f"sub{x}.loglike.npy") for x in [0,2,4]]
for idx, x in enumerate(d):
    print(idx, x.shape)

names = np.loadtxt("na_inclafr.samples", dtype=str)

names_dup = np.repeat(names, 2)
idx_d = np.array(["HG02006" in x for x in names_dup])
print(np.where(idx_d))
hap1_idx = np.where(idx_d)[0][0]
print(hap1_idx)
N_obss = [100*x for x in [1,2,4]]

def sum_subs(x1, x2):
    return x1[0::2,:,:] + x2[1::2,:,:]

d_sums_0 = sum_subs(d[1], d[1])
d_sums_2 = sum_subs(d[2], d[2])
d_sums_20 = sum_subs(d_sums_2, d_sums_2)


print(d_sums_0.shape)
print(d_sums_2.shape)
print(d_sums_20.shape)
print(d[0][0, hap1_idx])
print(d_sums_0[0, hap1_idx])
print(d_sums_20[0, hap1_idx])


fig, axes = plt.subplots(3,1, figsize=(18,6), sharex=True, sharey=True)
fig.suptitle(f"HG02006 hap1 max(loglike) per window")
for a in axes:
    a.margins(x=0)
axes[0].plot(range(100), d[0][:100, hap1_idx].max(1), c='red', label=f"sub0")
axes[0].set_title("sub0")
axes[1].plot(range(100), d_sums_0[:100, hap1_idx].max(1), c='green', label=f"sub2_0")
axes[1].set_title("sub2_0")
axes[2].plot(range(100), d_sums_20[:100, hap1_idx].max(1), c='blue', label=f"sub4_20")
axes[2].set_title("sub4_20")
fig.tight_layout()
fig.savefig("chr1_100_sumloglike_max.png", bbox_inches='tight', dpi=300)

fig, axes = plt.subplots(3,1, figsize=(18,6), sharex=True, sharey=True)
fig.suptitle(f"HG02006 hap1 sum(loglike) per window")
for a in axes:
    a.margins(x=0)
axes[0].plot(range(100), d[0][:100, hap1_idx].sum(1), c='red', label=f"sub0")
axes[0].set_title("sub0")
axes[1].plot(range(100), d_sums_0[:100, hap1_idx].sum(1), c='green', label=f"sub2_0")
axes[1].set_title("sub2_0")
axes[2].plot(range(100), d_sums_20[:100, hap1_idx].sum(1), c='blue', label=f"sub4_20")
axes[2].set_title("sub4_20")
fig.tight_layout()
fig.savefig("chr1_100_sumloglike_sum.png", bbox_inches='tight', dpi=300)

# fig, axes = plt.subplots(3,1, figsize=(24,8))
# fig.suptitle(f"HG02006 hap1")
# for a in axes:
#     a.margins(x=0)
#     a.axhline(y=0.95, color='black', alpha=0.5)

# for idx, sub in enumerate([0,2,4]):
#     ax = axes[idx]
#     ax.set_title(f"sub{sub}")
#         # ax.scatter(range(N_obss[idx]), d_test[idx][:N_obss[idx], 0, lab], s=3)
#     ax.plot(range(N_obss[idx]), d_test[idx][:N_obss[idx], 0, lab], label=f"K{lab}")
# axes[0].legend()
# fig.savefig("chr1_100_sumloglike.png", bbox_inches='tight', dpi=600)


# n_k = [(x[:-2], x[-2:]) for x in names_dup[idx_d]]
# column_names = [(hap,k) for (_,k),hap in zip(n_k, haps_test)]
# C = d.shape[2]
# missing_val = np.log(1.0/C)
# data = pd.DataFrame(~np.all(np.isclose(d_test, missing_val),axis=2)+0)
# data.columns = pd.MultiIndex.from_tuples(column_names)

# hap1 = data.loc[:2000, "hap1"]
# hap2 = data.loc[:2000, "hap2"]


# fig, ax = plt.subplots(2,1, figsize=(24,8))
# ax[0].set_title("hap1")
# for lab in hap1.columns:
#     ax[0].plot(hap1.index, hap1[lab], label=lab)
# ax[1].set_title("hap2")
# for lab in hap2.columns:
#     ax[1].plot(hap1.index, hap2[lab], label=lab)
# ax[0].legend()
# fig.savefig("test3.png", bbox_inches='tight', dpi=500)

