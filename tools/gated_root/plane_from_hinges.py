# Annulus plane through three hinge points read on cross-sections of the base frame.
# HINGES="h_cusp1,h_cusp2,h_cusp3" in the order of cusps.json 'cus'; each point sits on the lumen edge at that cusp's centre angle.
import numpy as np, json, os
from scipy.ndimage import map_coordinates, gaussian_filter
D = os.environ['CASE_DIR']; base = os.environ.get('BASE', 'frame_ann_nadir.json'); hk = [float(x) for x in os.environ['HINGES'].split(',')]
W = np.load(D + '/W.npy').astype(np.float32); src = gaussian_filter(W, 0.5)
F = json.load(open(D + '/' + base)); thr = F['thr']; c = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2'])
cus = json.load(open(D + '/cusps.json'))['cus']
s = lambda P: map_coordinates(src, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=1, cval=-1000)
RS = np.arange(0, 25, 0.05)
def edge(h, deg):   # lumen edge radius at height h and angle deg (median over +-10 deg), ignoring the first 4 mm
    rr = []
    for d in np.arange(deg - 10, deg + 10.1, 2.5):
        e = np.cos(np.radians(d)) * u1 + np.sin(np.radians(d)) * u2; v = s((c + A * h)[None] + RS[:, None] * e)
        b = (v < thr) & (RS > 4); rr.append(RS[np.argmax(b)] if b.any() else np.nan)
    return float(np.nanmedian(rr))
pts = []
for cd, h in zip(cus, hk):
    r = edge(min(h - 1.5, -0.5), cd); e = np.cos(np.radians(cd)) * u1 + np.sin(np.radians(cd)) * u2; pts.append(c + A * h + e * r); print(f'cusp {cd:5.1f}: hinge h {h:+.2f} mm, radius {r:.1f} mm')
pts = np.array(pts); m = pts.mean(0); n = np.cross(pts[1] - pts[0], pts[2] - pts[0]); n /= np.linalg.norm(n); n = n if n @ A > 0 else -n
up_ = np.array([0., -1., 0.]); v1 = up_ - n * (up_ @ n); v1 /= np.linalg.norm(v1); v2 = np.cross(n, v1); v2 = v2 if v2[0] > 0 else -v2
ang = np.degrees(np.arccos(min(1, n @ A)))
# direction of the correction: where the old plane sat lowest relative to the new one
th = np.radians(np.arange(0, 360, 5)); low = []
for t in th:
    e = np.cos(t) * u1 + np.sin(t) * u2; p = c + 12 * e; low.append((p - m) @ n)
k = int(np.argmin(low)); print(f'new plane differs by {ang:.1f} deg; the old plane sat lowest (toward the LV) at {np.degrees(th[k]):.0f} deg (0 = anterior, 90 = patient left), by {-min(low):.1f} mm at 12 mm radius')
F2 = dict(F, c=m.tolist(), A=n.tolist(), u1=v1.tolist(), u2=v2.tolist(), hinges=pts.tolist(), hinge_h_base=hk, plane_rule='three leaflet hinge points read on cross-sections (plane_from_hinges.py)')
json.dump(F2, open(D + '/frame_ann_hinge.json', 'w')); print('saved frame_ann_hinge.json; centre', m.round(2), 'normal', n.round(4))
