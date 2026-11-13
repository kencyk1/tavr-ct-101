# Semi-manual left main: seed placed on the vessel seen just outside the left sinus (h 4 mm, ~134 deg, 25 mm from the axis),
# then grown through bright non-lumen voxels; the ostium is where it meets the sinus lumen.
import numpy as np, json
from scipy.ndimage import distance_transform_edt, label as cc
from skimage.graph import MCP_Geometric
V = np.load('root_vol.npy').astype(np.float32); L = np.load('root_lab.npy'); G = np.load('root_geom.npz')
H = G['H'].astype(np.float32); dT = G['dT'].astype(np.float32); TH = G['TH'].astype(np.float32)
T = json.load(open('root_real_truth.json')); F = json.load(open('frame_ann.json')); O = np.array(T['cube_origin_world'])
c = np.array(F['c']) - O; A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2']); thr = F['thr']
seed_w = c + 4 * A + (-17.8) * u1 + 18.2 * u2; si = np.round(seed_w / 0.5 - .5).astype(int)[::-1]
print('seed HU', V[tuple(si)], 'label', L[tuple(si)])
lum = dT == 0; metal = V > 1100; nearM = distance_transform_edt(~metal) * 0.5 < 1.5
br = (V > thr) & ~lum & ~nearM
# move seed to the brightest voxel within 2 mm
r = 4; sub = V[si[0]-r:si[0]+r+1, si[1]-r:si[1]+r+1, si[2]-r:si[2]+r+1]; k = np.unravel_index(np.argmax(sub), sub.shape); si = si - r + np.array(k); print('seed moved to HU', V[tuple(si)])
cost = np.where(br, 1.0, np.inf); lo = np.maximum(si - 70, 0); hi = np.minimum(si + 70, 240); sl = tuple(slice(lo[i], hi[i]) for i in range(3))
dist, _ = MCP_Geometric(cost[sl]).find_costs([tuple(si - lo)]); keep = np.zeros(V.shape, bool); keep[sl] = np.isfinite(dist) & (dist * 0.5 < 22)
lumS = lum & (H > 0) & (H < T['stj_h'] + 6); dS = distance_transform_edt(~lumS) * 0.5
ct = keep & (dS <= 0.75)
print('vessel voxels', int(keep.sum()), 'contact voxels', int(ct.sum()))
if ct.any():
    hb = np.arange(0, 30, 0.5); cnt = [int((ct & (H >= h) & (H < h + 0.5)).sum()) for h in hb]
    print('contact by h:', ' '.join(f'{h:.1f}:{n}' for h, n in zip(hb, cnt) if n), ' angle', round(float(np.degrees(np.angle(np.exp(1j * np.radians(TH[ct]))).mean())) % 360))
np.save('lm_keep.npy', keep)
