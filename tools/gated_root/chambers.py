# Merge CADS heart labels (task 553) into the s1230 label volume, same rules as /home/claude/seg/chambers3.py.
import numpy as np, nibabel as nib, sys
from scipy.ndimage import distance_transform_edt, label as cc
D = '/home/claude/s1230'; k0 = int(sys.argv[1])
tov = lambda x: np.ascontiguousarray(x.transpose(2, 1, 0)[::-1, ::-1, ::-1])
W = np.load(D + '/W.npy').astype(np.float32); lab = np.load(D + '/lab_ts.npy')
c = np.asanyarray(nib.load(D + '/cads.nii.gz').dataobj).astype(np.uint8)
full = np.zeros(W.shape[::-1], np.uint8); full[:, :, k0:k0 + c.shape[2]] = c; C = tov(full)
MAP = {5: 4, 4: 96, 6: 94, 7: 95, 3: 97}
new = lab.copy(); new[lab == 22] = 0
heartish = np.isin(lab, [0, 22])
for k, i in MAP.items():
    m = (C == k) & heartish
    l, n = cc(m)
    if n > 1: sz = np.bincount(l.ravel()); sz[0] = 0; m = l == sz.argmax()
    new[m] = i
pa = (C == 8) & np.isin(lab, [0, 22]); new[pa] = 66
left = (lab == 22) & (new == 0); new[left & (W > -20)] = 97
zz, yy, xx = np.where(np.isin(new, [4, 94, 95, 96, 97]))
sl = (slice(max(zz.min() - 6, 0), zz.max() + 7), slice(max(yy.min() - 6, 0), yy.max() + 7), slice(max(xx.min() - 6, 0), xx.max() + 7)); o = new[sl]
mm = lambda m: distance_transform_edt(~m) * 1.5
dpa = mm(o == 66)
o[(o == 95) & (dpa < 14)] = 98
dLV, dRV, dRA, dLA, dRVOT = mm(o == 4), mm(o == 95), mm(o == 94), mm(o == 96), mm(o == 98)
o[(o == 97) & (dLV <= 9) & ((dRV <= 9) | (dRVOT <= 9))] = 99
o[((o == 94) & (dRV <= 2.5)) | ((o == 95) & (dRA <= 2.5))] = 100
o[(o == 98) & (dpa <= 2.5)] = 101
o[((o == 96) & (dLV <= 2.5)) | ((o == 4) & (dLA <= 2.5))] = 102
new[sl] = o
for i, n in ((4, 'LV'), (96, 'LA'), (94, 'RA'), (95, 'RV'), (98, 'RVOT'), (97, 'myo'), (99, 'IVS'), (66, 'PA'), (1, 'aorta')):
    m = new == i; print(f'{n:5s} {m.sum() * 3.375 / 1000:7.1f} mL', f'HU median {int(np.median(W[m]))}' if m.any() else '')
np.save(D + '/lab.npy', new)
