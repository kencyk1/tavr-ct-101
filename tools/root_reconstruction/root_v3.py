"""Root reconstruction v3 for the TAVR CT 101 MPR lab.
Builds a 0.5 mm, 120 mm cube around the aortic root of case ct1349 from the real 1.5 mm CTA:
fine-grain texture, curved centreline, crisp lumen with surrounding real tissue pulled in to meet it,
real LVOT, basal-ring ellipse at the annulus, modelled thickened calcified cusps, routed proximal coronaries.
Writes root_v3_vol.npy, root_v3_lab.npy, root_v3_truth.json and data/root_*.png."""
import numpy as np, json, os, sys
from scipy.ndimage import map_coordinates, gaussian_filter, gaussian_filter1d, median_filter
from PIL import Image

S = '/tmp/claude-0/-home-claude/41fb56d0-d570-55ae-901b-9576833dc018/scratchpad'
OUT = sys.argv[1] if len(sys.argv) > 1 else '/home/claude/ctp/data'
v = np.load(S + '/v.npy').astype(np.float32)                       # CTA [z,y,x], 1.5 mm, world = (i+.5)*1.5
labC = np.load('/home/claude/seg/labels_chambers.npy')
O_old = np.load(S + '/hyb_meta.npy')[:3]
u1, u2, A, c_old = np.load(S + '/hyb_frame.npy')                   # A: LV -> aorta, u1 anterior, u2 patient left
c_w = c_old + O_old                                                 # annulus centre, world mm

N, VS = 240, 0.5; Lc = N * VS
O = c_w + A * 14 - Lc / 2                                           # cube origin, world mm
c = c_w - O                                                         # annulus centre, cube mm
ax = ((np.arange(N) + .5) * VS).astype(np.float32)
Z, Y, X = np.meshgrid(ax, ax, ax, indexing='ij')
src = gaussian_filter(v, 0.55)
cw = [(Z + O[2]) / 1.5 - .5, (Y + O[1]) / 1.5 - .5, (X + O[0]) / 1.5 - .5]
up = map_coordinates(src, cw, order=3, mode='nearest').astype(np.float32)
labU = map_coordinates(labC, cw, order=0, mode='nearest').astype(np.uint8)
del cw

# ---- base texture: sharpened edges + fine isotropic grain ----
rng = np.random.default_rng(11)
nz = gaussian_filter(rng.standard_normal(up.shape).astype(np.float32), 0.7); nz /= nz.std()
lp = gaussian_filter(up, 2.5)
base = (up + 0.8 * (up - lp) + 15 * nz).astype(np.float32)
del lp

# ---- centreline + lumen radius Rs(h, theta), measured on the full CTA in world coords ----
def sampW(P, vol=src): P = np.asarray(P, float); return map_coordinates(vol, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=1, cval=-1000)
DT = 1.0; TH = np.radians(np.arange(0, 360, DT)); HS = np.arange(-20, 82.01, 0.25); RS = np.arange(0, 30, 0.2)
dirs = np.cos(TH)[:, None] * u1 + np.sin(TH)[:, None] * u2
off = np.zeros((len(HS), 2)); Rraw = np.zeros((len(HS), len(TH))); o_prev = np.zeros(2)
for i, h in enumerate(HS):
    thr = 250 if h >= -3 else (205 if h <= -7 else 250 - (250 - 205) * (-3 - h) / 4)
    o = o_prev.copy() if h >= 3 else np.zeros(2)
    for it in range(5):
        cen = c_w + A * h + o[0] * u1 + o[1] * u2
        P = cen[None, None, :] + dirs[:, None, :] * RS[None, :, None]
        val = sampW(P); b = val < thr; b[:, 0] = False
        r = np.where(b.any(1), RS[b.argmax(1)], RS[-1]); r = np.minimum(r, np.median(r) * (1.12 if h < 30 else 1.3))
        if h < 3: break
        pts = np.stack([r * np.cos(TH), r * np.sin(TH)], 1); o = o + pts.mean(0) * 0.8
    off[i] = o; Rraw[i] = r; o_prev = o
# offsets: zero at the annulus, measured above 10 mm
wo = np.clip((HS - 2) / 8, 0, 1)[:, None]
off = gaussian_filter1d(off, 16, axis=0, mode='nearest') * wo
Rs = median_filter(Rraw, (5, 15), mode='wrap')
Rs = gaussian_filter1d(gaussian_filter1d(Rs, 5, axis=1, mode='wrap'), 6, axis=0, mode='nearest')
iH = lambda h: int(np.argmin(abs(HS - h)))

# ---- basal-ring ellipse fitted to the real contour at h in [-2, 1] ----
r0 = Rs[[iH(h) for h in np.arange(-2, 1.01, 0.25)]].mean(0); F = np.fft.rfft(r0); F[7:] = 0; r0 = np.fft.irfft(F, len(r0))
x = r0 * np.cos(TH); y = r0 * np.sin(TH); xs, ys = np.roll(x, -1), np.roll(y, -1); cr = x * ys - xs * y; Ar = cr.sum() / 2
cx = np.sum((x + xs) * cr) / (6 * Ar); cy = np.sum((y + ys) * cr) / (6 * Ar)
X_, Y_ = x - cx, y - cy; Xs, Ys = np.roll(X_, -1), np.roll(Y_, -1); a = X_ * Ys - Xs * Y_
Iyy = np.sum(a * (Y_**2 + Y_ * Ys + Ys**2)) / 12 / Ar; Ixx = np.sum(a * (X_**2 + X_ * Xs + Xs**2)) / 12 / Ar
Ixy = np.sum(a * (X_ * Ys + 2 * X_ * Y_ + 2 * Xs * Ys + Xs * Y_)) / 24 / Ar
w, Vv = np.linalg.eigh(np.array([[Ixx, Ixy], [Ixy, Iyy]])); sa, sb = 2 * np.sqrt(w[1]), 2 * np.sqrt(w[0]); ea = Vv[:, 1]
k = np.sqrt(Ar / (np.pi * sa * sb)); sa *= k; sb *= k; phi = np.arctan2(ea[1], ea[0])
def ell_r(th):
    d = np.stack([np.cos(th), np.sin(th)]); cp, sp = np.cos(phi), np.sin(phi)
    tr = lambda q: np.stack([(q[0] * cp + q[1] * sp) / sa, (-q[0] * sp + q[1] * cp) / sb])
    D = tr(d); Oo = tr(np.array([cx, cy])[:, None] * np.ones_like(th))
    aa = (D**2).sum(0); bb = -2 * (D * Oo).sum(0); cq = (Oo**2).sum(0) - 1
    return (-bb + np.sqrt(bb * bb - 4 * aa * cq)) / (2 * aa)
re_ = ell_r(TH)
for i, h in enumerate(HS):
    we = np.clip(1 - abs(h) / 3.0, 0, 1); Rs[i] = we * re_ + (1 - we) * Rs[i]
np.save(S + '/v3_Rs.npy', Rs); np.save(S + '/v3_off.npy', off)

# ---- per-voxel root coordinates (curved centreline) ----
dX, dY, dZ = X - c[0], Y - c[1], Z - c[2]
H = (dX * A[0] + dY * A[1] + dZ * A[2]).astype(np.float32)
hi = np.clip((H - HS[0]) / 0.25, 0, len(HS) - 1)
o1 = np.interp(hi.ravel(), np.arange(len(HS)), off[:, 0]).reshape(H.shape).astype(np.float32)
o2 = np.interp(hi.ravel(), np.arange(len(HS)), off[:, 1]).reshape(H.shape).astype(np.float32)
Q1 = (dX * u1[0] + dY * u1[1] + dZ * u1[2] - o1).astype(np.float32); Q2 = (dX * u2[0] + dY * u2[1] + dZ * u2[2] - o2).astype(np.float32)
del dX, dY, dZ
Rr = np.hypot(Q1, Q2).astype(np.float32); THd = (np.degrees(np.arctan2(Q2, Q1)) % 360).astype(np.float32)
RsV = map_coordinates(Rs, [hi, THd / DT], order=1, mode='grid-wrap').astype(np.float32)
d = (Rr - RsV).astype(np.float32)
wh = np.clip(np.minimum((H + 18) / 4, (54 - H) / 7), 0, 1).astype(np.float32)
STJ_H = None
mR = Rs.mean(1); hs_sinus = HS[iH(6) + int(np.argmax(mR[iH(6):iH(18)]))]
for h in np.arange(hs_sinus + 6, 40, 0.25):
    if mR[iH(h)] <= mR[iH(h + 0.25)]: STJ_H = float(h); break

# ---- outside tissue pulled in to meet the new edge ----
band = (d >= 0) & (d <= 6) & (wh > 0)
idx = np.nonzero(band)
dd = d[idx]; r_src = RsV[idx] + dd + 3.5 * (1 - dd / 6)
radial_x = (Q1[idx] * u1[0] + Q2[idx] * u2[0]) / np.maximum(Rr[idx], 1e-3)
radial_y = (Q1[idx] * u1[1] + Q2[idx] * u2[1]) / np.maximum(Rr[idx], 1e-3)
radial_z = (Q1[idx] * u1[2] + Q2[idx] * u2[2]) / np.maximum(Rr[idx], 1e-3)
shift = r_src - Rr[idx]
Px = X[idx] + radial_x * shift; Py = Y[idx] + radial_y * shift; Pz = Z[idx] + radial_z * shift
outv = base.copy()
outv[idx] = map_coordinates(base, [Pz / VS - .5, Py / VS - .5, Px / VS - .5], order=1, mode='nearest')
del Px, Py, Pz, radial_x, radial_y, radial_z, shift, r_src, dd

# ---- compose lumen / wall / outside ----
lum = 358.0
lowv = np.clip(gaussian_filter(up, 4) - lum, -30, 30) * 0.5
fl = np.clip((H + 14) / 8, 0, 1)
lumen_val = fl * (lum + lowv + 18 * nz) + (1 - fl) * (gaussian_filter(up, 1.2) + 14 * nz)
tw = np.clip((H + 1) / 2, 0, 1) * 1.2
wall_val = 0.5 * np.clip(outv, -60, 130) + 0.5 * 95 + 6 * nz
sig = lambda t: 1 / (1 + np.exp(-np.clip(t, -60, 60)))
a_l = sig(-d / 0.3); a_w = np.clip(sig(-(d - tw) / 0.3) - a_l, 0, 1); a_w = np.where(tw > 0.05, a_w, 0); a_o = 1 - a_l - a_w
vol = a_l * lumen_val + a_w * wall_val + a_o * outv
vol = wh * vol + (1 - wh) * base
realCa = (up > 600) & (d > -1.5) & (d < 4) & (wh > 0)
vol = np.where(realCa, np.maximum(vol, base), vol).astype(np.float32)
del a_l, a_w, a_o, wall_val, lowv, outv

# ---- labels ----
lab = labU.copy(); lab[np.isin(lab, [1, 2, 3, 5, 7])] = 0
inl = (d < 0) & (wh > 0.5)
lab[inl & (H > STJ_H)] = 1; lab[inl & (H >= 0) & (H <= STJ_H)] = 2; lab[inl & (H < 0) & (H > -15)] = 7
lab[inl & (np.abs(H) < 0.6)] = 5

# ---- coronaries (routes from the earlier build, moved into this cube) ----
rt = np.load(S + '/routes.npz'); sh_ = O_old - O
def tubes(paths):
    dc = np.full(X.shape, 99, np.float32); lc = np.zeros(X.shape, np.uint8)
    for P, r0_, r1_, lv in paths:
        P = np.asarray(P, float); seg = np.linalg.norm(np.diff(P, axis=0), axis=1); s = np.concatenate([[0], np.cumsum(seg)]); L = s[-1]
        rad = lambda t: r0_ + (r1_ - r0_) * np.clip(t / L, 0, 1)
        for k in range(len(P) - 1):
            a_, b_ = P[k], P[k + 1]; rm = max(rad(s[k]), rad(s[k + 1])) + 3
            lo = np.maximum((np.minimum(a_, b_) - rm) / VS, 0).astype(int); hi_ = np.minimum((np.maximum(a_, b_) + rm) / VS + 1, N).astype(int)
            sl = (slice(lo[2], hi_[2]), slice(lo[1], hi_[1]), slice(lo[0], hi_[0]))
            vv = b_ - a_; wx, wy, wz = X[sl] - a_[0], Y[sl] - a_[1], Z[sl] - a_[2]
            t = np.clip((wx * vv[0] + wy * vv[1] + wz * vv[2]) / max(vv @ vv, 1e-9), 0, 1)
            dist = np.sqrt((wx - t * vv[0])**2 + (wy - t * vv[1])**2 + (wz - t * vv[2])**2) - rad(s[k] + t * seg[k])
            bt = dist < dc[sl]; dc[sl] = np.where(bt, dist, dc[sl]); lc[sl] = np.where(bt, lv, lc[sl])
    return dc, lc
LAD = np.vstack([rt['lm_in'], rt['LAD']]) + sh_; LCX = np.vstack([rt['lm_in'], rt['LCX']]) + sh_; RCA = np.vstack([rt['rca_in'], rt['RCA']]) + sh_
dco, lco = tubes([(LAD, 2.1, 1.5, 8), (LCX, 2.1, 1.4, 8), (RCA, 1.9, 1.5, 9)])
outc = (dco < 2) & (d > 0)
ac_l = sig(-dco / 0.3); ac_w = np.clip(sig(-(dco - 0.7) / 0.3) - ac_l, 0, 1)
vol = np.where(outc, vol * (1 - ac_l - ac_w) + ac_l * (lum + 18 * nz) + ac_w * (55 + 6 * nz), vol)
lab[(dco < 0) & (d > 0)] = lco[(dco < 0) & (d > 0)]

# ---- cusps: thickened, restricted, fitted to the real sinus walls ----
def Rlook(h, thd): return map_coordinates(Rs, [np.clip((h - HS[0]) / 0.25, 0, len(HS) - 1), (thd % 360) / DT], order=1, mode='grid-wrap')
Hc = min(21.0, 0.8 * STJ_H); Hf = 15.5
cen = {'RCC': -23.0, 'LCC': 118.0, 'NCC': -105.0}
comm = {('RCC', 'LCC'): 47.5, ('LCC', 'NCC'): 186.5, ('NCC', 'RCC'): 296.0}
bounds = {'RCC': (comm[('NCC', 'RCC')] - 360, comm[('RCC', 'LCC')]), 'LCC': (comm[('RCC', 'LCC')], comm[('LCC', 'NCC')]), 'NCC': (comm[('LCC', 'NCC')] - 360, comm[('NCC', 'RCC')] - 360)}
roi = (H > -1.5) & (H < Hc + 1) & (Rr < 27); ix = np.nonzero(roi); h_ = H[ix]; th_ = THd[ix]; r_ = Rr[ix]
sd = np.full(h_.shape, 99.0, np.float32); par = {}; ssb = np.full(h_.shape, 9.0, np.float32); cub = np.zeros(h_.shape, np.int8)
for nm in ('RCC', 'LCC', 'NCC'):
    (lo_, hi_), ctr = bounds[nm], cen[nm]; ddg = ((th_ - ctr + 180) % 360) - 180; inside = (ddg >= lo_ - ctr) & (ddg <= hi_ - ctr)
    half = np.where(ddg < 0, ctr - lo_, hi_ - ctr); tt = np.clip(np.abs(ddg) / half, 0, 1); hatt = Hc * (1 - np.cos(np.pi / 2 * tt)); ok = inside & (hatt < Hf - 1.0)
    ss = (h_ - hatt) / np.maximum(Hf - hatt, 1.0); ok &= (ss >= 0) & (ss <= 1); ratt = Rlook(hatt, th_); rfree = 4.2 + 1.0 * np.cos(np.pi * tt)
    sc = np.clip(ss, 0, 1) ** 0.75; rl = ratt * (1 - sc) + rfree * sc
    TL = 1.9 + 0.9 * (1 - np.clip(ss, 0, 1))**3 + 0.8 * np.clip(ss, 0, 1)**4
    cand = np.abs(r_ - rl) - TL / 2; better = ok & (cand < sd)
    sd = np.where(better, cand, sd); ssb = np.where(better, ss, ssb); cub = np.where(better, {'RCC': 1, 'LCC': 2, 'NCC': 3}[nm], cub); par[nm] = (lo_, hi_, ctr)
aL = np.zeros(vol.shape, np.float32); aL[ix] = np.clip(0.5 - sd / 0.35, 0, 1); aL = np.clip(gaussian_filter(aL, 0.7) * 1.1, 0, 1)
vol = vol * (1 - aL) + (112 + 8 * nz) * aL; lab[aL > 0.4] = 3

# ---- calcium: irregular plaques that follow the leaflet surface, hinge points, LVOT extension ----
rng = np.random.default_rng(5)
gn = gaussian_filter(rng.standard_normal(vol.shape).astype(np.float32), 2.2); gn /= gn.std()
caL = np.zeros(vol.shape, np.float32); hu = np.zeros(vol.shape, np.float32)
thr_c = {1: 0.75, 2: 0.95, 3: 0.45}                                            # NCC heaviest, then RCC, LCC
g_ix = gn[ix]
for cid, t_ in thr_c.items():
    m = (cub == cid) & (sd < 0.75) & (ssb > 0.04) & (ssb < 0.75) & (g_ix > t_)
    caL[tuple(np.array(ix)[:, m])] = 1; hu[tuple(np.array(ix)[:, m])] = 750 + 550 * np.clip((g_ix[m] - t_) / 0.9, 0, 1)
ncc_c, rcc_c = 255.0, 337.0
angd = lambda a_, b_: np.abs(((a_ - b_ + 180) % 360) - 180)
hinge = (d > -1.6) & (d < 0.4) & (np.abs(H) < 1.6) & (((angd(THd, ncc_c) < 16) & (gn > 0.2)) | ((angd(THd, rcc_c) < 12) & (gn > 0.5)))
lvot = (d > -1.4) & (d < 0.4) & (H > -6) & (H < 0.5) & (THd > 182) & (THd < 214) & (gn > 0.1 + 0.12 * np.clip(-H, 0, 6))
for m, base_hu in ((hinge, 820), (lvot, 760)):
    caL[m] = 1; hu[m] = np.maximum(hu[m], base_hu + 450 * np.clip(gn[m], 0, 1))
aC = np.clip(gaussian_filter(caL, 0.55) * 1.25, 0, 1)
hC = gaussian_filter(hu, 0.55) / np.maximum(gaussian_filter(caL, 0.55), 1e-3)
vol = np.maximum(vol, vol * (1 - aC) + hC * aC); lab[aC > 0.45] = 10
del gn, caL, hu

vol = np.clip(vol, -1024, 3071).astype(np.int16)

# ---- truth ----
def ostium(P, r0_):
    rel = P - c; h = rel @ A; hi_ = np.clip((h - HS[0]) / 0.25, 0, len(HS) - 1)
    oo = np.stack([np.interp(hi_, np.arange(len(HS)), off[:, j]) for j in (0, 1)], 1)
    q1 = rel @ u1 - oo[:, 0]; q2 = rel @ u2 - oo[:, 1]; th = np.degrees(np.arctan2(q2, q1)) % 360; rr = np.hypot(q1, q2)
    dd_ = rr - map_coordinates(Rs, [hi_, th / DT], order=1, mode='grid-wrap'); k = int(np.argmax(dd_ > 0)); f = -dd_[k - 1] / (dd_[k] - dd_[k - 1])
    return round(float(h[k - 1] + f * (h[k] - h[k - 1]) - r0_), 1)
per = float(np.pi * (3 * (sa + sb) - np.sqrt((3 * sa + sb) * (sa + 3 * sb))))
cc = c + cx * u1 + cy * u2
T = dict(cube_origin_world=[float(q) for q in O], shape=[N, N, N], voxel_mm=VS, c=[round(float(q), 2) for q in cc], n_toLV=[round(float(-q), 4) for q in A],
         max=round(2 * sa, 1), min=round(2 * sb, 1), area=round(float(np.pi * sa * sb)), per=round(per, 1),
         LM=ostium(LAD, 2.1), RCA=ostium(RCA, 1.9), sinus_d=round(float(2 * mR[iH(hs_sinus)]), 1), stj_h=STJ_H, stj_d=round(float(2 * mR[iH(STJ_H)]), 1))
json.dump(T, open(S + '/root_v3_truth.json', 'w'), indent=1); print(json.dumps(T))
np.save(S + '/root_v3_vol.npy', vol); np.save(S + '/root_v3_lab.npy', lab)
u = (vol.astype(np.int32) + 1024).astype(np.uint16); files = []
for i, z0 in enumerate(range(0, N, 64)):
    z1 = min(z0 + 64, N); cu = u[z0:z1].reshape(-1, N); ll = lab[z0:z1].reshape(-1, N)
    rgb = np.stack([(cu >> 8).astype(np.uint8), (cu & 255).astype(np.uint8), ll], -1); fn = f'{OUT}/root_{i:02d}.png'
    Image.fromarray(rgb, 'RGB').save(fn, optimize=True, compress_level=9); files.append([z0, z1, os.path.getsize(fn)])
print(json.dumps(files), sum(f[2] for f in files) / 1e6, 'MB')
