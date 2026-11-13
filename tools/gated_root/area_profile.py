# Lumen area/perimeter (SCCT-style contour) on planes parallel to the frame plane at several offsets.
import numpy as np, json, os
from scipy.ndimage import map_coordinates, gaussian_filter, median_filter
D = os.environ['CASE_DIR']; F = json.load(open(D + '/frame_ann.json')); hs = [float(x) for x in os.environ.get('HS', '-3,-2.5,-2,-1.5,-1,-0.5,0,1').split(',')]
W = np.load(D + '/W.npy').astype(np.float32); src = gaussian_filter(W, 0.5); thr = F['thr']; c0 = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2'])
TT = np.radians(np.arange(0, 360, 1.0)); RR = np.arange(0, 25, 0.02); dd = np.cos(TT)[:, None] * u1 + np.sin(TT)[:, None] * u2
for h in hs:
    cen = c0 + A * h
    for it in range(4):
        P = cen[None, None] + dd[:, None] * RR[None, :, None]; v = map_coordinates(src, [P[..., k] / 1.5 - .5 for k in (2, 1, 0)], order=1, cval=-1000)
        b = (v < thr); b[:, :50] = False; r = np.where(b.any(1), RR[b.argmax(1)], RR[-1]); r = np.minimum(r, np.median(r) * 1.15); r = median_filter(r, 15, mode='wrap')
        x = r * np.cos(TT); y = r * np.sin(TT); xs, ys = np.roll(x, -1), np.roll(y, -1); cr = x * ys - xs * y; Ar = cr.sum() / 2
        cx = np.sum((x + xs) * cr) / (6 * Ar); cy = np.sum((y + ys) * cr) / (6 * Ar); cen = cen + cx * u1 + cy * u2
    per = np.sum(np.hypot(np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0]))))
    print(f'h {h:+5.1f}: area {Ar:6.0f} mm2  perimeter {per:5.1f} mm  min/max radius {r.min():4.1f}/{r.max():4.1f}')
