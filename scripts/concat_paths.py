import os
import sys
import numpy as np

def load_files(f):
    res = []
    with open(f, 'r') as fh:
        for line in fh:
            res.append(line.rstrip())
    return res


path_filelist = load_files(sys.argv[1])
out= sys.argv[2]

L_list = []
for f in path_filelist:
    L_list.append(np.loadtxt(f, dtype=int).T)
    print("path", L_list[-1].shape)
L = np.concatenate(L_list, axis=0)
print(L.shape)
np.savetxt(fname=out, X=L, fmt="%d")

