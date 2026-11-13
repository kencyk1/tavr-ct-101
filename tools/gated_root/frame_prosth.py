# Frame for a case with a surgical (stented) bioprosthesis: the valve axis and sewing-ring plane come from the metal frame.
# Writes frame.json / frame_ann.json / cusps.json in CASE_DIR, plus diagnostic images.
import numpy as np, json, os
from scipy.ndimage import map_coordinates, gaussian_filter, binary_dilation, label as cc
from PIL import Image, ImageDraw
D = os.environ['CASE_DIR']; MET = float(os.environ.get('MET_HU', 1200))
W = np.load(D + '/W.npy').astype(np.float32); lab = np.load(D + '/lab.npy')
toW = lambda idx: (idx[:, ::-1] + .5) * 1.5
itf = (lab == 1) & binary_dilation(lab == 4, iterations=3)
c0 = toW(np.argwhere(itf)).mean(0) if itf.any() else toW(np.argwhere(lab == 1)).mean(0)
m = W > MET; l, n = cc(m, structure=np.ones((3, 3, 3)))
P_all = toW(np.argwhere(m)); near = np.linalg.norm(P_all - c0, axis=1) < 40
ids = l[m][near]; best = np.bincount(ids).argmax() if len(ids) else 0
P = toW(np.argwhere(l == best)); print('metal component voxels', len(P), 'distance of its centroid from aorta/LV interface', round(float(np.linalg.norm(P.mean(0) - c0)), 1))
ctr = P.mean(0); ev, V = np.linalg.eigh(np.cov((P - ctr).T)); A = V[:, 0]                 # smallest spread = valve axis
ao = toW(np.argwhere(lab == 1)); ao = ao[np.linalg.norm(ao - ctr, axis=1) < 45].mean(0)
lv = toW(np.argwhere(lab == 4)); lv = lv[np.linalg.norm(lv - ctr, axis=1) < 50].mean(0)
if A @ (ao - lv) < 0: A = -A
print('eigenvalues', ev.round(1), 'axis (x,y,z)', A.round(3))
# refine axis/centre: least-squares cylinder (iterate centre in the plane, keep axis from PCA, then re-estimate axis from ring/post geometry)
def frame(A):
    up_ = np.array([0., -1., 0.]); u1 = up_ - A * (up_ @ A); u1 /= np.linalg.norm(u1); u2 = np.cross(A, u1); u2 = u2 if u2[0] > 0 else -u2; return u1, u2
for it in range(5):
    u1, u2 = frame(A); q = P - ctr; x, y, h = q @ u1, q @ u2, q @ A
    # circle fit (algebraic) of the frame cross-section
    M_ = np.stack([x, y, np.ones_like(x)], 1); bb = x * x + y * y; sol = np.linalg.lstsq(M_, bb, rcond=None)[0]; cx, cy = sol[0] / 2, sol[1] / 2
    ctr = ctr + cx * u1 + cy * u2
    # axis refinement: points of the lowest 6 mm (sewing ring) define the ring plane
    q = P - ctr; h = q @ A; low = h < np.percentile(h, 5) + 6
    Q = q[low]; Q0 = Q.mean(0); _, _, Vt = np.linalg.svd(Q - Q0); nA = Vt[2] if Vt[2] @ A > 0 else -Vt[2]
    ang = np.degrees(np.arccos(min(1, nA @ A))); A = A + 0.5 * (nA - A); A /= np.linalg.norm(A); print(f' it{it}: centre shift {np.hypot(cx, cy):.2f} mm, ring-plane vs axis {ang:.1f} deg')
u1, u2 = frame(A); q = P - ctr; h = q @ A; r = np.hypot(q @ u1, q @ u2)
h0 = float(np.percentile(h, 3)); c = ctr + A * h0
hp = h - h0; print('frame: radius median', round(float(np.median(r)), 1), 'mm; post tips', round(float(np.percentile(hp, 99)), 1), 'mm above the inflow; inflow ring at h0 =', round(h0, 1))
src = gaussian_filter(W, 0.5)
lum = float(np.median(map_coordinates(src, [((c + A * 10 + np.random.default_rng(0).normal(0, 2, (300, 3)))[:, k] / 1.5 - .5) for k in (2, 1, 0)], order=1)))
thr = (lum + 60) / 2; print('lumen HU', round(lum), 'threshold', round(thr))
F = dict(c=c.tolist(), A=A.tolist(), u1=u1.tolist(), u2=u2.tolist(), thr=thr, lum=lum, frame_r=float(np.median(r)), post_h=float(np.percentile(hp, 99)), met=MET)
json.dump(F, open(D + '/frame.json', 'w')); json.dump(F, open(D + '/frame_ann.json', 'w')); json.dump(dict(com=[], cus=[120.0, 240.0, 0.0]), open(D + '/cusps.json', 'w'))
s = lambda Pp: map_coordinates(W, [Pp[..., 2] / 1.5 - .5, Pp[..., 1] / 1.5 - .5, Pp[..., 0] / 1.5 - .5], order=3, cval=-1000)
def img(cen, e1, e2, half, px, W_=1400, L_=500):
    t = np.arange(-half, half, px); a_, b_ = np.meshgrid(t, t, indexing='xy'); v = s(cen[None, None] + a_[..., None] * e1 + (-b_)[..., None] * e2)
    return Image.fromarray(np.clip((v - (L_ - W_ / 2)) / W_ * 255, 0, 255).astype(np.uint8)).convert('RGB')
ims = []
for hh in (-4, -2, 0, 2, 5, 8, 12, 16):
    im = img(c + A * hh, u2, u1, 22, 0.15); ImageDraw.Draw(im).text((3, 3), f'h {hh}', fill=(255, 255, 0)); ims.append(im)
for ang in (0, 60, 120):
    e1 = np.cos(np.radians(ang)) * u1 + np.sin(np.radians(ang)) * u2; im = img(c + A * 8, e1, A, 22, 0.15); d = ImageDraw.Draw(im)
    y = (22 + 8) / 0.15; d.line([(0, y), (293, y)], fill=(0, 255, 255)); d.text((3, 3), f'long {ang}', fill=(255, 255, 0)); ims.append(im)
g = Image.new('RGB', (295 * 4, 295 * 3)); [g.paste(im, ((i % 4) * 295, (i // 4) * 295)) for i, im in enumerate(ims)]; g.save(D + '/diag_prosth.png')
