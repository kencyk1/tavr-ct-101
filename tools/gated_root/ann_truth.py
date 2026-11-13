# Annulus answer key on the real 1.5 mm scan: lumen contour at the half-maximum edge, with calcium kept outside the contour
import numpy as np, json, os, sys
from scipy.ndimage import map_coordinates, gaussian_filter, median_filter
D = os.environ['CASE_DIR']; CA = float(os.environ.get('CA_HU', 650))
W = np.load(D + '/W.npy').astype(np.float32); src = gaussian_filter(W, 0.5)
F = json.load(open(D + '/frame_ann.json')); thr = F['thr']; c_w = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2'])
T = json.load(open(D + '/root_real_truth.json')); O = np.array(T['cube_origin_world'])
TT = np.radians(np.arange(0, 360, 0.5)); RR = np.arange(0, 25, 0.02); dd_ = np.cos(TT)[:, None] * u1 + np.sin(TT)[:, None] * u2
def contour(excl):
    cen = c_w.copy()
    for it in range(6):
        P = cen[None, None] + dd_[:, None] * RR[None, :, None]
        val = map_coordinates(src, [P[..., k] / 1.5 - .5 for k in (2, 1, 0)], order=3, cval=-1000)
        b = (val < thr); b[:, :50] = False; r = np.where(b.any(1), RR[b.argmax(1)], RR[-1]); r = np.minimum(r, np.median(r) * 1.15)
        if excl:   # SCCT: draw the contour 'as if no calcium is present' -> rays whose edge meets calcium are bridged from their neighbours
            ie = np.clip((r / 0.02).astype(int), 0, len(RR) - 1); hit = np.array([(val[i, max(ie[i] - 100, 0):min(ie[i] + 100, len(RR))] > CA).any() for i in range(len(TT))])
            ok = ~hit
            if ok.sum() > 20:
                ang = np.arange(len(TT)); r = np.interp(ang, np.concatenate([ang[ok] - len(TT), ang[ok], ang[ok] + len(TT)]), np.tile(r[ok], 3))
            calc_frac = hit.mean()
        r = median_filter(r, 31, mode='wrap'); Fr = np.fft.rfft(r); Fr[9:] = 0; r = np.fft.irfft(Fr, len(r))
        x = r * np.cos(TT); y = r * np.sin(TT); xs, ys = np.roll(x, -1), np.roll(y, -1); cr = x * ys - xs * y; Ar = cr.sum() / 2
        cx = np.sum((x + xs) * cr) / (6 * Ar); cy = np.sum((y + ys) * cr) / (6 * Ar); cen = cen + cx * u1 + cy * u2
    X_, Y_ = x - cx, y - cy; Xs, Ys = np.roll(X_, -1), np.roll(Y_, -1); a_ = X_ * Ys - Xs * Y_
    Iyy = np.sum(a_ * (Y_**2 + Y_ * Ys + Ys**2)) / 12 / Ar; Ixx = np.sum(a_ * (X_**2 + X_ * Xs + Xs**2)) / 12 / Ar; Ixy = np.sum(a_ * (X_ * Ys + 2 * X_ * Y_ + 2 * Xs * Ys + Xs * Y_)) / 24 / Ar
    w_, _ = np.linalg.eigh(np.array([[Ixx, Ixy], [Ixy, Iyy]])); sa, sb = 2 * np.sqrt(w_[1]), 2 * np.sqrt(w_[0]); k = np.sqrt(Ar / (np.pi * sa * sb)); sa *= k; sb *= k
    per = float(np.sum(np.hypot(np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0])))))
    return dict(c=[round(float(q), 2) for q in cen - O], max=round(2 * sa, 1), min=round(2 * sb, 1), area=round(float(Ar)), per=round(per, 1)), (x - cx, y - cy)
inc, _ = contour(False); exc, xy = contour(True)
print('calcium inside contour :', inc); print('contour as if no calcium (SCCT):', exc)
if len(sys.argv) > 1 and sys.argv[1] == 'write':
    T.update(exc); T['annulus_rule'] = f'half-maximum lumen edge ({round(thr)} HU); where the edge meets calcium (>{CA:.0f} HU, within 2 mm) the contour is bridged smoothly as if no calcium were present (SCCT 2019)'
    json.dump(T, open(D + '/root_real_truth.json', 'w'), indent=1); np.save(D + '/ann_contour.npy', np.array(xy)); print('written')
