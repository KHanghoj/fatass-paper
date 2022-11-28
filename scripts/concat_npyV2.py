import os
import sys
import numpy as np

out = sys.argv[1]
npy_filelist = sys.argv[2:]

L_list = []
for f in npy_filelist:
    L_list.append(np.load(f))
    print("npy", L_list[-1].shape)
L = np.concatenate(L_list, axis=0)
print(L.shape)
np.save(file=out, arr=L)

