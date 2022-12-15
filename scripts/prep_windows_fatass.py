import sys
import numpy as np
import argparse

# Argparse
parser = argparse.ArgumentParser(prog="Prepare windows for Fatass")
parser.add_argument("-w", "--windows", help="Genomic position per window")
parser.add_argument("-L0", help="LogLike no subsplit", default=None)
parser.add_argument("-L2", help="LogLike 2 subsplit", default=None)
parser.add_argument("-L4", help="LogLike 4 subsplit", default=None)
parser.add_argument("-o", "--out", help="path to output")
parser.add_argument("--normalizer", help="Scaling basepair distance between windows", default=1e6,type=float)
args = parser.parse_args()
if any([args.L0, args.L2, args.L4])==True and all([args.L0, args.L2, args.L4])==False:
    print("Paths to all Loglike arguments (L{0,2,4}) to enable fancy masking or none of them") 
    sys.exit(1)


if args.L0 is None:
    print("Same windows for all haplotypes")
    windows = np.loadtxt(args.windows, usecols=1).astype(np.float32)
    windows = np.diff(windows, prepend=0,n=1).astype(np.float32)
    chrom_skip = windows<0
    ## put chromsoome jumps to large value
    windows[windows<0] = 1e9
    ## normalizing the windows on the same chromosomes
    windows[~chrom_skip] = windows[~chrom_skip]/args.normalizer
    np.savetxt(args.out, windows)
else:
    print("individual windows per haplotype")
    windows = np.loadtxt(args.windows, usecols=1).astype(np.float32)
    chrom_skip = np.diff(windows, prepend=0,n=1).astype(np.float32)<0
    l0 = np.load(args.L0).max(axis=2)
    
    l2_temp = np.load(args.L2).max(axis=2)
    l2 = np.zeros(l0.shape)
    for x in range(2):
        l2[:] += l2_temp[x::2]
    del l2_temp
    l4_temp = np.load(args.L4).max(axis=2)
    (W4, N) = l4_temp.shape
    l4 = np.zeros(l0.shape)
    for x in range(4):
        l4[:] += l4_temp[x::4]
    del l4_temp

    assert W4 == windows.shape[0], "number of windows must match subsplit4 loglike"
    windows_per_hap_bool = np.zeros((W4,N), dtype=int)
    windows_per_hap = np.zeros((W4,N), dtype=float)
    # windows_per_hap_bool[0] = 1 # must start with a distance
    windows_per_hap_bool[3::4] = 1 # must allow transition after each window
    windows_per_hap_bool[chrom_skip] = 1 # must allow transition at chromosome breaks
    log = np.zeros((N,3), dtype=int)
    for i in range(N):
        res = (l2[:,i]>l0[:,i]+200).astype(int) + (l4[:,i]>l2[:,i]+200).astype(int)
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
    windows_per_hap[chrom_skip] = 1e9
    print(f"dumping data into {args.out}")
    np.savetxt(args.out+".log", log, fmt="%d")
    np.savetxt(args.out, windows_per_hap, fmt="%.5f")
    np.savetxt(args.out+".bool", windows_per_hap_bool, fmt="%d")
