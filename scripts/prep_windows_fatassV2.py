
import sys
import numpy as np
import argparse

# Argparse
parser = argparse.ArgumentParser(prog="Prepare windows for Fatass")
parser.add_argument("-w", "--windows", help="Genomic position per window")
parser.add_argument("-L", help="LogLikes in incremental order (0,2,4,8,16 ... )",
                    default=None, nargs="+")
parser.add_argument("-o", "--out", help="path to output")
parser.add_argument("--normalizer", help="Scaling basepair distance between windows",
                    default=1e6,type=float)
parser.add_argument("--loglike_diff", help="Minimum loglike units improvement for one subsplit to the next",
                    default=100, type=float)
args = parser.parse_args()

windows = np.loadtxt(args.windows, usecols=[1,2]).astype(np.float32).mean(1)

if args.L is None:
    print("Same windows for all haplotypes")
    windows = np.diff(windows, prepend=0,n=1).astype(np.float32)
    chrom_skip = windows<0
    ## put chromsoome jumps to large value
    windows[windows<0] = 1e9
    ## normalizing the windows on the same chromosomes
    windows[~chrom_skip] = windows[~chrom_skip]/args.normalizer
    np.savetxt(args.out, windows)
else:
    print("individual windows per haplotype")
    nsplits = [x for x in range(0, 1000) if x & (x-1)==0 and x!=1]
    split_keep = nsplits[:len(args.L)]
    max_split = max(split_keep)
    chrom_skip = np.diff(windows, prepend=0,n=1).astype(np.float32)<0

    ls_temp = [np.load(x).max(axis=2) for x in args.L]
    (Wmax, N) = ls_temp[-1].shape
    (Wmin, _) = ls_temp[0].shape
    nll = len(ls_temp)
    assert Wmax == windows.shape[0], "number of windows must match highest subsplit loglike"

    sumlogs = [np.zeros((Wmin, N)) for _ in range(len(split_keep))]
    for idx, (sub, data) in enumerate(zip(split_keep, ls_temp)):
        print(sub)
        print(data.shape)
        if idx == 0:
            sumlogs[idx] = data        
        else:
            for x in range(sub):
                sumlogs[idx][:] += data[x::sub]

    windows_per_hap_bool = np.zeros((Wmax,N), dtype=int)
    windows_per_hap = np.zeros((Wmax,N), dtype=float)
    # windows_per_hap_bool[0] = 1 # must start with a distance
    windows_per_hap_bool[(max_split-1)::max_split] = 1 # must allow transition after each window
    windows_per_hap_bool[chrom_skip] = 1 # must allow transition at chromosome breaks
    log = np.zeros((N,nll), dtype=int)
    res = np.zeros(Wmin, dtype=float)
    for i in range(N):
        for (sl1, sl2) in zip(sumlogs, sumlogs[1:]):
            res[:] += (sl2[:,i] > (sl1[:,i]+args.loglike_diff)).astype(int)
        print(f"{i+1}/{N}", end='\r')
        log[i] = np.unique(res, return_counts=1)[-1] # count the differences
        for idx, x in enumerate(res):
            start_p = idx * 4
            if x == 1: # break between sub2 windows
                windows_per_hap_bool[start_p+1, i] = 1
            elif x== 2: # break between sub4 windows
                for x_idx in range(start_p, start_p+4):
                    windows_per_hap_bool[x_idx, i] = 1
        keep = windows_per_hap_bool[:,i]==1
        windows_per_hap[keep,i] = np.diff(windows[keep], prepend=0,n=1)/args.normalizer
        res[:] = 0
    windows_per_hap[chrom_skip] = 1e9
    print(f"dumping data into {args.out}")
    np.savetxt(args.out+".counts", log, fmt="%d")
    np.savetxt(args.out, windows_per_hap, fmt="%.5f")
    np.savetxt(args.out+".bool", windows_per_hap_bool, fmt="%d")
