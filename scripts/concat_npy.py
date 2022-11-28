import os
import sys
import numpy as np

def load_files(f):
    res = []
    with open(f, 'r') as fh:
        for line in fh:
            res.append(line.rstrip())
    return res


npy_filelist = load_files(sys.argv[1])
out= sys.argv[2]


# jonas_dir = "/home/jonas/results/haplonet/greenland"
L_list = []
for f in npy_filelist:
    L_list.append(np.load(f))
    print("npy", L_list[-1].shape)
L = np.concatenate(L_list, axis=0)
print(L.shape)
np.save(file=out, arr=L)

