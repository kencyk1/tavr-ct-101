"""Real gated root case for the TAVR CT 101 MPR lab (case s1230, TotalSegmentator dataset, CC BY 4.0).
A 120 mm cube around the aortic root, resampled from the real 1.5 mm ECG-gated TAVR CT to 0.5 mm (cubic, mild sharpening,
light grain). Nothing is modelled: cusps, calcium and coronaries are the patient's own. Labels are automatic.
Writes data/root_*.png (0.5 mm cube), data/rootfb_*.png (1.5 mm surroundings) and root_real_truth.json."""
import numpy as np, json, os, sys
from scipy.ndimage import map_coordinates, gaussian_filter, gaussian_filter1d, median_filter, label as cc, binary_dilation
from PIL import Image
import os; D = os.environ.get('CASE_DIR', '/home/claude/s1230'); OUT = sys.argv[1] if len(sys.argv) > 1 else '/home/claude/ctp/data'
W = np.load(D + '/W.npy').astype(np.float32); LAB = np.load(D + '/lab.npy')
F = json.load(open(D + '/frame_ann.json')); thr = F['thr']
c_w = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2'])
cus = json.load(open(D + '/cusps.json'))['cus']
_ad = lambda a, b: abs(((a - b + 180) % 360) - 180)
cusR = min(cus, key=lambda q: _ad(q, 0)); cusL = min([q for q in cus if q != cusR], key=lambda q: _ad(q, 120))
print('cusp centres', [round(q) for q in cus], 'RCC', round(cusR), 'LCC', round(cusL))
N, VS = 240, 0.5; Lc = N * VS
O = c_w + A * 14 - Lc / 2; c = c_w - O
ax = ((np.arange(N) + .5) * VS).astype(np.float32); Z, Y, X = np.meshgrid(ax, ax, ax, indexing='ij')
cw = [(Z + O[2]) / 1.5 - .5, (Y + O[1]) / 1.5 - .5, (X + O[0]) / 1.5 - .5]
up = map_coordinates(W, cw, order=3, mode='nearest').astype(np.float32)
labU = map_coordinates(LAB, cw, order=0, mode='nearest').astype(np.uint8); del cw
rng = np.random.default_rng(7)
nz = gaussian_filter(rng.standard_normal(up.shape).astype(np.float32), 0.7); nz /= nz.std()
vol = up + 0.45 * (up - gaussian_filter(up, 2.0)) + 7 * nz

# ---- lumen envelope along the curved centreline (labels and answer key only; the image is not altered) ----
src = gaussian_filter(W, 0.5)
def sampW(P): P = np.asarray(P, float); return map_coordinates(src, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=1, cval=-1000)
DT = 1.0; TH = np.radians(np.arange(0, 360, DT)); HS = np.arange(-20, 82.01, 0.25); RS = np.arange(0, 32, 0.2)
dirs = np.cos(TH)[:, None] * u1 + np.sin(TH)[:, None] * u2
off = np.zeros((len(HS), 2)); Rraw = np.zeros((len(HS), len(TH))); o_prev = np.zeros(2)
for i, h in enumerate(HS):
    o = o_prev.copy() if h >= 3 else np.zeros(2)
    for it in range(5):
        cen = c_w + A * h + o[0] * u1 + o[1] * u2
        val = sampW(cen[None, None, :] + dirs[:, None, :] * RS[None, :, None]); b = val < thr; b[:, 0] = False
        r = np.where(b.any(1), RS[b.argmax(1)], RS[-1]); r = np.minimum(r, np.median(r) * (1.15 if h < 30 else 1.3))
        if h < 3: break
        o = o + np.stack([r * np.cos(TH), r * np.sin(TH)], 1).mean(0) * 0.8
    off[i] = o; Rraw[i] = r; o_prev = o
off = gaussian_filter1d(off, 16, axis=0, mode='nearest') * np.clip((HS - 2) / 8, 0, 1)[:, None]
Rs = median_filter(Rraw, (5, 15), mode='wrap'); Rs = gaussian_filter1d(gaussian_filter1d(Rs, 4, axis=1, mode='wrap'), 3, axis=0, mode='nearest')
iH = lambda h: int(np.argmin(abs(HS - h)))
mR = Rs.mean(1); hs_sinus = HS[iH(5) + int(np.argmax(mR[iH(5):iH(22)]))]; STJ_H = None
for h in np.arange(hs_sinus + 4, 45, 0.25):
    if mR[iH(h)] <= mR[iH(h + 0.25)] + 0.005: STJ_H = float(h); break
print('sinus max at h', hs_sinus, 'diam', round(2 * mR[iH(hs_sinus)], 1), ' STJ h', STJ_H, 'diam', round(2 * mR[iH(STJ_H)], 1))

dX, dY, dZ = X - c[0], Y - c[1], Z - c[2]
H = (dX * A[0] + dY * A[1] + dZ * A[2]).astype(np.float32)
hi = np.clip((H - HS[0]) / 0.25, 0, len(HS) - 1)
o1 = np.interp(hi.ravel(), np.arange(len(HS)), off[:, 0]).reshape(H.shape).astype(np.float32)
o2 = np.interp(hi.ravel(), np.arange(len(HS)), off[:, 1]).reshape(H.shape).astype(np.float32)
Q1 = (dX * u1[0] + dY * u1[1] + dZ * u1[2] - o1).astype(np.float32); Q2 = (dX * u2[0] + dY * u2[1] + dZ * u2[2] - o2).astype(np.float32)
del dX, dY, dZ, o1, o2
Rr = np.hypot(Q1, Q2).astype(np.float32); THd = (np.degrees(np.arctan2(Q2, Q1)) % 360).astype(np.float32)
d = (Rr - map_coordinates(Rs, [hi, THd / DT], order=1, mode='grid-wrap')).astype(np.float32)
wh = np.clip(np.minimum((H + 18) / 4, (54 - H) / 7), 0, 1).astype(np.float32)

# ---- labels ----
# true lumen: bright voxels, morphologically opened (2.5 mm) to cut coronaries and partial-volume bridges, connected to the root
from scipy.ndimage import distance_transform_edt, maximum_filter
br_l = gaussian_filter(up, 0.6) > thr
core = distance_transform_edt(br_l) * VS > 2.5
opened = (distance_transform_edt(~core) * VS <= 2.5) & br_l
lc0, _ = cc(opened); cidx = tuple(np.round((c + A * 8) / VS - .5).astype(int)[::-1]); lumT = lc0 == lc0[cidx]
dT = (distance_transform_edt(~lumT) * VS).astype(np.float32)
lab = labU.copy(); lab[np.isin(lab, [1, 2, 3, 5, 7])] = 0
inl = lumT & (wh > 0.5) & (Rr < 32) & (H > -15) & (d < 3.5) & ~np.isin(labU, [94, 95, 96, 98, 66, 49, 48, 40])
lab[(labU == 1) & ~inl & (H > 20)] = 1
lab[inl & (H > STJ_H)] = 1; lab[inl & (H >= 0) & (H <= STJ_H)] = 2; lab[inl & (H < 0)] = 7
leaf = ~lumT & (H > 0.4) & (H < STJ_H - 2) & (d < -0.8) & (up < thr)          # real cusp tissue inside the root
lab[leaf] = 3
lab[inl & (np.abs(H) < 0.6)] = 5
CA_T = max(600.0, float(F.get('lum', 350)) + 350); ca = (up > CA_T) & (d < 3) & (H > -8) & (H < STJ_H) & (Rr < 30); lab[ca] = 10

# ---- coronaries: a bright tube leaving its sinus (reaches >= 4 mm outside the true lumen), then followed <= 40 mm ----
from skimage.graph import MCP_Geometric
chamb = np.isin(labU, [4, 94, 95, 96, 98, 66, 48, 49, 40, 41, 99, 100, 101, 102])
br = up > 200
thick = maximum_filter(distance_transform_edt(br) * VS, size=9)
thin = br & (thick <= 2.8) & ~chamb & ~lumT
if os.environ.get('EXCL_METAL'):   # prosthesis cases: ignore bright blooming around the metal frame
    nearM = distance_transform_edt(~(up > 1100)) * VS < 3.5; thin &= ~nearM; chamb = chamb | nearM
angd = lambda a, b: np.abs(((a - b + 180) % 360) - 180)
_hr = lambda k, d: tuple(float(q) for q in os.environ.get(k, d).split(','))
SIDES = {'L': dict(th=cusL, tw=40, h=_hr('LM_H', '0,24'), lv=8), 'R': dict(th=cusR - 10, tw=45, h=_hr('RCA_H', '4,40'), lv=9)}
cor = np.zeros(up.shape, np.uint8); heights = {}
for side, P_ in SIDES.items():
    sect = (angd(THd, P_['th']) < P_['tw']) & (H > P_['h'][0]) & (H < P_['h'][1])
    chan = (up > thr) & ~lumT & (dT < 8) & sect & ~chamb
    lc, nl = cc(chan, structure=np.ones((3, 3, 3)))
    best = None
    for k in range(1, nl + 1):
        m = lc == k; dm = dT[m]
        if dm.min() <= 1.0 and dm.max() >= 4 and (best is None or m.sum() > best[1]): best = (k, int(m.sum()))
    if best is None: print('no ostium found for', side); continue
    m = lc == best[0]; ost = m & (dT <= 1.0)
    hv = H[ost]; print(side, 'ostium voxels', int(ost.sum()), 'H range', hv.min().round(1), hv.max().round(1), 'angle', round(float(np.degrees(np.angle(np.exp(1j*np.radians(THd[ost]))).mean())) % 360))
    hmin = float(hv.min()); heights[side] = hmin
    P_['ost_th'] = float(np.degrees(np.angle(np.exp(1j * np.radians(THd[ost]))).mean())) % 360
    seeds = np.argwhere(m & (dT > 2))
    cost = np.where(thin | m, 1.0, np.inf)
    lo = np.maximum(seeds.min(0) - 90, 0); hi_ = np.minimum(seeds.max(0) + 90, N)
    sl = tuple(slice(lo[i], hi_[i]) for i in range(3))
    dist, _ = MCP_Geometric(cost[sl]).find_costs(seeds - lo)
    keep = np.zeros(up.shape, bool); keep[sl] = np.isfinite(dist) & (dist * VS < 40)
    keep |= m
    cor[keep & (cor == 0)] = P_['lv']
    print(f'coronary {side}: ostium lower edge h {hmin:.2f}  traced voxels {int(keep.sum())}')
lab[(cor > 0) & (lab != 10)] = cor[(cor > 0) & (lab != 10)]
# ostium lower edge: lowest 0.5 mm level of the continuous run where the traced vessel opens into the sinus/STJ lumen
lumS = lumT & (H > 0) & (H < STJ_H + 8); dS = distance_transform_edt(~lumS) * VS
for side, P_ in SIDES.items():
    if side not in heights: continue
    m = (cor == P_['lv']) & (dS <= 0.75) & (angd(THd, P_.get('ost_th', P_['th'])) < 25) & (H > 2); hb = np.arange(-2, 40, 0.5); cnt = np.array([int((m & (H >= h) & (H < h + 0.5)).sum()) for h in hb])
    run = cnt >= 8; best = None; i = 0
    while i < len(hb):
        if run[i]:
            j = i
            while j + 1 < len(hb) and cnt[j + 1] >= 3: j += 1
            if best is None or cnt[i:j + 1].sum() > best[1]: best = (hb[i], cnt[i:j + 1].sum())
            i = j + 1
        else: i += 1
    if best is None: print('run rule found nothing for', side); continue
    heights[side] = float(best[0]); print(f'coronary {side}: ostium lower edge (run rule) h {best[0]:.1f}')
vol = np.clip(vol, -1024, 3071).astype(np.int16)

# ---- answer key: annulus contour at h = 0 on the real image ----
TT = np.radians(np.arange(0, 360, 0.5)); RR = np.arange(0, 25, 0.02); dd_ = np.cos(TT)[:, None] * u1 + np.sin(TT)[:, None] * u2
cen = c_w.copy()
for it in range(6):
    val = map_coordinates(src, [((cen[None, None] + dd_[:, None] * RR[None, :, None])[..., k] / 1.5 - .5) for k in (2, 1, 0)], order=3, cval=-1000)
    b = val < thr; b[:, :50] = False; r = np.where(b.any(1), RR[b.argmax(1)], RR[-1])
    r = np.minimum(r, np.median(r) * 1.15); r = median_filter(r, 31, mode='wrap')
    Fr = np.fft.rfft(r); Fr[9:] = 0; r = np.fft.irfft(Fr, len(r))
    x = r * np.cos(TT); y = r * np.sin(TT); xs, ys = np.roll(x, -1), np.roll(y, -1); cr = x * ys - xs * y; Ar = cr.sum() / 2
    cx = np.sum((x + xs) * cr) / (6 * Ar); cy = np.sum((y + ys) * cr) / (6 * Ar); cen = cen + cx * u1 + cy * u2
X_, Y_ = x - cx, y - cy; Xs, Ys = np.roll(X_, -1), np.roll(Y_, -1); a_ = X_ * Ys - Xs * Y_
Iyy = np.sum(a_ * (Y_**2 + Y_ * Ys + Ys**2)) / 12 / Ar; Ixx = np.sum(a_ * (X_**2 + X_ * Xs + Xs**2)) / 12 / Ar; Ixy = np.sum(a_ * (X_ * Ys + 2 * X_ * Y_ + 2 * Xs * Ys + Xs * Y_)) / 24 / Ar
w_, Vv = np.linalg.eigh(np.array([[Ixx, Ixy], [Ixy, Iyy]])); sa, sb = 2 * np.sqrt(w_[1]), 2 * np.sqrt(w_[0]); kk = np.sqrt(Ar / (np.pi * sa * sb)); sa *= kk; sb *= kk
perim = float(np.sum(np.hypot(np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0])))))
# Feret (caliper) diameters of the real contour, for reference
angs = np.radians(np.arange(0, 180, 0.5)); proj = np.cos(angs)[:, None] * x[None] + np.sin(angs)[:, None] * y[None]; fer = proj.max(1) - proj.min(1)
cc_cube = cen - O
T = dict(case=os.path.basename(D), cube_origin_world=[float(q) for q in O], shape=[N, N, N], voxel_mm=VS, c=[round(float(q), 2) for q in cc_cube],
         n_toLV=[round(float(-q), 4) for q in A], max=round(2 * sa, 1), min=round(2 * sb, 1), area=round(float(Ar)), per=round(perim, 1),
         feret_max=round(float(fer.max()), 1), feret_min=round(float(fer.min()), 1),
         LM=round(heights.get('L', np.nan), 1), RCA=round(heights.get('R', np.nan), 1),
         sinus_d=round(float(2 * mR[iH(hs_sinus)]), 1), stj_h=STJ_H, stj_d=round(float(2 * mR[iH(STJ_H)]), 1))
json.dump(T, open(D + '/root_real_truth.json', 'w'), indent=1); print(json.dumps(T))
np.save(D + '/root_vol.npy', vol); np.save(D + '/root_lab.npy', lab); np.savez_compressed(D + '/root_geom.npz', H=H.astype(np.float16), d=d.astype(np.float16), TH=THd.astype(np.float16), dT=dT.astype(np.float16))

def strips(v16, ll, prefix, nzs=64):
    u = (v16.astype(np.int32) + 1024).astype(np.uint16); files = []; Nz, Ny, Nx = v16.shape
    for i, z0 in enumerate(range(0, Nz, nzs)):
        z1 = min(z0 + nzs, Nz); cu = u[z0:z1].reshape(-1, Nx); l_ = ll[z0:z1].reshape(-1, Nx)
        rgb = np.stack([(cu >> 8).astype(np.uint8), (cu & 255).astype(np.uint8), l_], -1); fn = f'{OUT}/{prefix}_{i:02d}.png'
        Image.fromarray(rgb, 'RGB').save(fn, optimize=True, compress_level=9); files.append([z0, z1, os.path.getsize(fn)])
    return files
PFX = os.environ.get('PREFIX', 'root'); f_root = strips(vol, lab, PFX)
# ---- 1.5 mm surroundings for panning beyond the cube (same patient) ----
ctr = (O + Lc / 2); half = 130
lo = np.maximum(np.floor((ctr - half) / 1.5).astype(int), 0); hi_ = np.minimum(np.ceil((ctr + half) / 1.5).astype(int), np.array(W.shape[::-1]))
fbv = np.clip(W[lo[2]:hi_[2], lo[1]:hi_[1], lo[0]:hi_[0]], -1024, 3071).astype(np.int16); fbl = LAB[lo[2]:hi_[2], lo[1]:hi_[1], lo[0]:hi_[0]].copy()
fbl[np.isin(fbl, [5, 7, 3])] = 0
f_fb = strips(fbv, fbl, PFX + 'fb')
Fo = lo * 1.5
meta = dict(root_files=f_root, fb_files=f_fb, fb_shape=list(fbv.shape), fb_off=[round(float(q), 3) for q in (O - Fo)])
json.dump(meta, open(D + '/root_real_files.json', 'w')); print(json.dumps(meta)); print('MB', sum(f[2] for f in f_root) / 1e6, sum(f[2] for f in f_fb) / 1e6)
