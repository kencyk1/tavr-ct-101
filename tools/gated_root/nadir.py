import numpy as np, json, os
from scipy.ndimage import map_coordinates, gaussian_filter, gaussian_filter1d
W = np.load('W.npy').astype(np.float32); src = gaussian_filter(W, 0.5)
F = json.load(open('frame.json')); thr = F['thr']
s = lambda P: map_coordinates(src, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=1, cval=-1000)
def frame(A):
    up_ = np.array([0., -1., 0.]); u1 = up_ - A * (up_ @ A); u1 /= np.linalg.norm(u1); u2 = np.cross(A, u1); u2 = u2 if u2[0] > 0 else -u2; return u1, u2
TH = np.radians(np.arange(0, 360, 2.0)); RS = np.arange(0, 30, 0.1); HS = np.arange(-8, 12.01, 0.25)
def radii(c, A):
    u1, u2 = frame(A); dirs = np.cos(TH)[:, None] * u1 + np.sin(TH)[:, None] * u2; R = np.zeros((len(HS), len(TH)))
    for i, h in enumerate(HS):
        val = s((c + A * h)[None, None] + dirs[:, None] * RS[None, :, None]); b = val < thr; b[:, :5] = False
        r = np.where(b.any(1), RS[b.argmax(1)], RS[-1]); R[i] = np.minimum(r, np.median(r) * 1.25)
    return u1, u2, R
c = np.array(F['c']); A = np.array(F['A']); cus = json.load(open('cusps.json'))['cus']
for it in range(int(os.environ.get("NIT", 4))):
    u1, u2, R = radii(c, A)
    Rm = np.median(R, axis=1)
    pts = []
    for cd in cus:
        k = int(round(cd / 2)) % len(TH); sec = [(k + j) % len(TH) for j in range(-8, 9)]       # +-16 deg around cusp centre
        rc = gaussian_filter1d(np.median(R[:, sec], axis=1), 2)
        base = np.median(rc[(HS >= -7) & (HS <= -3)])
        hn = HS[np.argmax((rc > base + 1.5) & (HS > -5))]                                         # sinus starts to bulge
        dirc = np.cos(np.radians(cd)) * u1 + np.sin(np.radians(cd)) * u2
        r_at = float(np.interp(hn, HS, rc)); pts.append(c + A * hn + dirc * r_at)
        print(f'  it{it} cusp {cd:5.0f}: LVOT r {base:5.2f}  nadir h {hn:5.2f}  r {r_at:5.2f}')
    pts = np.array(pts); m = pts.mean(0); nrm = np.cross(pts[1] - pts[0], pts[2] - pts[0]); nrm /= np.linalg.norm(nrm); nrm = nrm if nrm @ A > 0 else -nrm
    ang = np.degrees(np.arccos(min(1, nrm @ A))); print(f' it{it}: plane tilt vs axis {ang:.1f} deg; centre shift along axis {(m - c) @ A:.2f} mm')
    damp = float(os.environ.get('DAMP', 1.0))
    A_new = A + damp * (nrm - A); A_new /= np.linalg.norm(A_new); c = c + damp * (m - c); A = A_new
u1, u2 = frame(A)
json.dump(dict(F, c=c.tolist(), A=A.tolist(), u1=u1.tolist(), u2=u2.tolist(), nadirs=pts.tolist()), open('frame_ann.json', 'w'))
print('annulus centre', c.round(2), 'normal (LV->aorta)', A.round(4))
