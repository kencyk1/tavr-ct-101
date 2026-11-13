# Stage 1 for the real-root MPR case (case s1230, TotalSegmentator dataset, CC BY 4.0).
# Puts the CT and automatic labels into viewer orientation (z head->feet, y anterior->posterior, x patient right->left,
# world = (index + 0.5) * 1.5 mm) and merges TotalSegmentator (fast) + CADS heart labels into page label ids.
import numpy as np, nibabel as nib, json, sys
from totalsegmentator.map_to_binary import class_map
import os; D = os.environ.get('CASE_DIR', '/home/claude/s1230')
tov = lambda x: np.ascontiguousarray(x.transpose(2, 1, 0)[::-1, ::-1, ::-1])
im = nib.load(D + '/ct.nii.gz'); assert nib.aff2axcodes(im.affine) == ('R', 'A', 'S')
W = tov(np.asanyarray(im.dataobj).astype(np.int16)); np.save(D + '/W.npy', W)
ts = tov(np.asanyarray(nib.load(D + '/ts_fast.nii.gz').dataobj).astype(np.uint8))
# group table from the main pipeline (same ids and names as the page)
src = open('/home/claude/seg/merge.py').read(); src = src[src.index('G=[]'):src.index('inv={v:k')]
ns = {}; exec(src, ns); G = ns['G']
inv = {v: k for k, v in class_map['total'].items()}
lab = np.zeros(W.shape, np.uint8)
for gid, n, s, c, names, show in G:
    for nm in names:
        if nm in inv: lab[ts == inv[nm]] = gid
lab[ts == inv['aorta']] = 1
for nm, gid in (('iliac_artery_right', 7), ('iliac_artery_left', 6)):
    if nm in inv: lab[ts == inv[nm]] = gid
np.save(D + '/lab_ts.npy', lab)
print(W.shape, {int(k): int(v) for k, v in zip(*np.unique(lab, return_counts=True)) if k in (1, 22, 66, 6, 7)})
z = np.where((lab == 22).any((1, 2)))[0]; print('heart z', z.min(), z.max(), ' RAS k range', W.shape[0] - 1 - z.max(), W.shape[0] - 1 - z.min())
