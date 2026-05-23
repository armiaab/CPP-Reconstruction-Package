"""
Extract data from .mat files and write as binary double arrays
for the C++ test program to load.
"""
import scipy.io
import numpy as np
import struct, os

root = os.path.dirname(os.path.abspath(__file__))

def save_doubles(filename, arr):
    arr = np.asarray(arr, dtype=np.float64).ravel()
    with open(filename, 'wb') as f:
        f.write(struct.pack('<Q', len(arr)))
        f.write(arr.tobytes())

# OUdata1D — used by Examples 1,2,3,7,8 in tutorial
mat = scipy.io.loadmat(os.path.join(root, 'OUdata1D.mat'))
ou_data = mat['data'].ravel()
print(f"OUdata1D: {len(ou_data)} points, min={ou_data.min():.4f}, max={ou_data.max():.4f}")
save_doubles(os.path.join(root, 'OUdata1D.bin'), ou_data)

# MayData1D — used by Examples 4,5 in tutorial
mat = scipy.io.loadmat(os.path.join(root, 'MayData1D.mat'))
may_data = mat['data'].ravel()
print(f"MayData1D: {len(may_data)} points, min={may_data.min():.4f}, max={may_data.max():.4f}")
save_doubles(os.path.join(root, 'MayData1D.bin'), may_data)

print("Done.")
