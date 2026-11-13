# Valve-in-valve answer key for s0972 (stented surgical bioprosthesis): frame inner contour at the sewing ring,
# post height, coronary ostium heights above the inflow (lowest metal) and their radial distance from the valve axis.
import numpy as np, json, os
from scipy.ndimage import map_coordinates, gaussian_filter, median_filter, distance_transform_edt
from PIL import Image
D = '/home/claude/s0972'; OUT = '/home/claude/s0972/data'
W = np.load(D + '/W.npy').astype(np.float32); src = gaussian_filter(W, 0.5)
F = json.load(open(D + '/frame_ann.json')); T = json.load(open(D + '/root_real_truth.json'))
c_w = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2']); lum = F['lum']; O = np.array(T['cube_origin_world'])
P = W > 1200; pts = (np.argwhere(P)[:, ::-1] + .5) * 1.5; q = pts - c_w; hp = q @ A; rp = np.hypot(q @ u1, q @ u2); fr = (hp > -2) & (hp < 20) & (rp < 18)
met_peak = float(np.percentile(W[P][fr[: len(W[P])]] if False else src[tuple(np.argwhere(P)[fr].T)], 75)); edge = (lum + met_peak) / 2
print('metal peak (75th pct, smoothed)', round(met_peak), ' inner-edge threshold', round(edge))
TT = np.radians(np.arange(0, 360, 1.0)); RR = np.arange(0, 18, 0.02); dd = np.cos(TT)[:, None] * u1 + np.sin(TT)[:, None] * u2
rin = []
for h in np.arange(0.5, 4.01, 0.5):
    Pp = (c_w + A * h)[None, None] + dd[:, None] * RR[None, :, None]
    v = map_coordinates(src, [Pp[..., k] / 1.5 - .5 for k in (2, 1, 0)], order=3, cval=-1000); b = v > edge; b[:, :100] = False
    rin.append(np.where(b.any(1), RR[b.argmax(1)], np.nan))
rin = np.array(rin); r = np.nanmedian(rin, 0); ok = np.isfinite(r); ang = np.arange(len(TT))
r = np.interp(ang, np.concatenate([ang[ok] - 360, ang[ok], ang[ok] + 360]), np.tile(r[ok], 3)); r = median_filter(r, 15, mode='wrap')
Fr = np.fft.rfft(r); Fr[9:] = 0; r = np.fft.irfft(Fr, len(r))
x = r * np.cos(TT); y = r * np.sin(TT); xs, ys = np.roll(x, -1), np.roll(y, -1); cr = x * ys - xs * y; Ar = cr.sum() / 2
cx = np.sum((x + xs) * cr) / (6 * Ar); cy = np.sum((y + ys) * cr) / (6 * Ar)
X_, Y_ = x - cx, y - cy; Xs, Ys = np.roll(X_, -1), np.roll(Y_, -1); a_ = X_ * Ys - Xs * Y_
Iyy = np.sum(a_ * (Y_**2 + Y_ * Ys + Ys**2)) / 12 / Ar; Ixx = np.sum(a_ * (X_**2 + X_ * Xs + Xs**2)) / 12 / Ar; Ixy = np.sum(a_ * (X_ * Ys + 2 * X_ * Y_ + 2 * Xs * Ys + Xs * Y_)) / 24 / Ar
w_, _ = np.linalg.eigh(np.array([[Ixx, Ixy], [Ixy, Iyy]])); sa, sb = 2 * np.sqrt(w_[1]), 2 * np.sqrt(w_[0]); k = np.sqrt(Ar / (np.pi * sa * sb)); sa *= k; sb *= k
per = float(np.sum(np.hypot(np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0])))))
cen = c_w + cx * u1 + cy * u2
print(f'frame inner contour: {2*sa:.1f} x {2*sb:.1f} mm, area {Ar:.0f} mm2, perimeter {per:.1f} mm, area-derived diameter {2*np.sqrt(Ar/np.pi):.1f} mm')
# coronary ostia: lower edge height and radial distance from the valve axis at that level
V = np.load(D + '/root_vol.npy').astype(np.float32); L = np.load(D + '/root_lab.npy'); G = np.load(D + '/root_geom.npz'); H = G['H'].astype(np.float32); dT = G['dT'].astype(np.float32)
keep = np.load(D + '/lm_keep.npy'); lumS = (dT == 0) & (H > 0) & (H < T['stj_h'] + 6); dS = distance_transform_edt(~lumS) * 0.5
N = V.shape[0]; ax = (np.arange(N) + .5) * .5
def ostium(mask, h_lo):
    ct = mask & (dS <= 0.75); idx = np.argwhere(ct); Pw = np.stack([ax[idx[:, 2]], ax[idx[:, 1]], ax[idx[:, 0]]], 1) + O
    qq = Pw - cen; hh = qq @ A; rr = np.hypot(qq @ u1, qq @ u2); sel = (hh >= h_lo) & (hh < h_lo + 1.5)
    return float(np.median(rr[sel])), float(np.median(np.degrees(np.arctan2(qq[sel] @ u2, qq[sel] @ u1)) % 360))
LM_h = 6.5; RCA_h = float(T['RCA'])
rL, aL = ostium(keep, LM_h); rR, aR = ostium(L == 9, RCA_h)
Rv = float(F['frame_r'])
print(f'LM: lower edge {LM_h} mm, ostium {rL:.1f} mm from the valve axis at {aL:.0f} deg -> VTC (to a {2*Rv:.0f} mm virtual valve) {rL-Rv:.1f} mm')
print(f'RCA: lower edge {RCA_h} mm, ostium {rR:.1f} mm from the axis at {aR:.0f} deg -> VTC {rR-Rv:.1f} mm')
T.update(c=[round(float(q), 2) for q in cen - O], max=round(2 * sa, 1), min=round(2 * sb, 1), area=round(float(Ar)), per=round(per, 1), LM=LM_h, RCA=RCA_h,
         frame_centreline_d=round(2 * Rv, 1), post_h=round(float(F['post_h']), 1), LM_r=round(rL, 1), RCA_r=round(rR, 1), VTC_LM=round(rL - Rv, 1), VTC_RCA=round(rR - Rv, 1),
         annulus_rule=f'inner edge of the metal frame (threshold {round(edge)} HU, half-way between lumen and metal) in the sewing-ring band 0.5-4 mm above the lowest metal',
         plane_rule='inflow plane = lowest metal of the frame (3rd percentile along the valve axis); axis from the frame geometry')
json.dump(T, open(D + '/root_real_truth.json', 'w'), indent=1)
# labels: replace the automatic left coronary with the seeded trace; frame metal and calcium = 10
lab = L.copy(); lab[lab == 8] = 0; lab[keep & (lab != 10)] = 8; lab[V > 1100] = 10
np.save(D + '/root_lab.npy', lab)
u = (np.load(D + '/root_vol.npy').astype(np.int32) + 1024).astype(np.uint16); files = []
for i, z0 in enumerate(range(0, N, 64)):
    z1 = min(z0 + 64, N); cu = u[z0:z1].reshape(-1, N); ll = lab[z0:z1].reshape(-1, N)
    rgb = np.stack([(cu >> 8).astype(np.uint8), (cu & 255).astype(np.uint8), ll], -1); fn = f'{OUT}/root3_{i:02d}.png'
    Image.fromarray(rgb, 'RGB').save(fn, optimize=True, compress_level=9); files.append([z0, z1, os.path.getsize(fn)])
M = json.load(open(D + '/root_real_files.json')); M['root_files'] = files; json.dump(M, open(D + '/root_real_files.json', 'w')); print(files)
print({k: T[k] for k in ['c', 'n_toLV', 'max', 'min', 'area', 'per', 'LM', 'RCA', 'post_h', 'VTC_LM', 'VTC_RCA', 'sinus_d', 'stj_h', 'stj_d']})
