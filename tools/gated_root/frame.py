# Stage 3: root frame from the real scan. Initial plane = aorta/LV interface of the automatic labels,
# then the axis is refined from lumen centroids (LVOT -> sinus). Writes diagnostic long-axis and cross-section images.
import numpy as np, json, sys
from scipy.ndimage import map_coordinates, gaussian_filter, binary_dilation
from PIL import Image, ImageDraw
import os; D = os.environ.get('CASE_DIR', '/home/claude/s1230')
W = np.load(D + '/W.npy').astype(np.float32); lab = np.load(D + '/lab.npy')
src = gaussian_filter(W, 0.5)
def sampW(P, order=1):
    P = np.asarray(P, float); return map_coordinates(src, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=order, cval=-1000)
toW = lambda idx: (idx[:, ::-1] + .5) * 1.5            # (z,y,x) index -> (x,y,z) mm
itf = (lab == 1) & binary_dilation(lab == 4, iterations=2)
P = toW(np.argwhere(itf)); c0 = P.mean(0)
_, _, Vt = np.linalg.svd(P - c0); n = Vt[2]
ao = toW(np.argwhere(lab == 1)); ao = ao[np.linalg.norm(ao - c0, axis=1) < 30].mean(0)
lv = toW(np.argwhere(lab == 4)); lv = lv[np.linalg.norm(lv - c0, axis=1) < 40].mean(0)
if n @ (ao - lv) < 0: n = -n
print('interface pts', len(P), 'c0', c0.round(1), 'n', n.round(3))
lum = float(np.median(sampW(c0 + n * 12 + np.random.default_rng(0).normal(0, 2, (200, 3)))))
thr = (lum + 60) / 2; print('lumen HU', round(lum), 'threshold', round(thr))
TH = np.radians(np.arange(0, 360, 2.0)); RS = np.arange(0, 30, 0.2)
def frame(A):
    up_ = np.array([0., -1., 0.]); u1 = up_ - A * (up_ @ A); u1 /= np.linalg.norm(u1)   # anterior
    u2 = np.cross(A, u1); u2 = u2 if u2[0] > 0 else -u2                                  # patient left
    return u1, u2
def section(c, A, h, u1, u2, it=4):
    o = np.zeros(2)
    for _ in range(it):
        cen = c + A * h + o[0] * u1 + o[1] * u2
        dirs = np.cos(TH)[:, None] * u1 + np.sin(TH)[:, None] * u2
        val = sampW(cen[None, None] + dirs[:, None] * RS[None, :, None]); b = val < thr; b[:, 0] = False
        r = np.where(b.any(1), RS[b.argmax(1)], RS[-1]); r = np.minimum(r, np.median(r) * 1.3)
        pts = np.stack([r * np.cos(TH), r * np.sin(TH)], 1); o = o + pts.mean(0)
    return cen, r
c, A = c0.copy(), n.copy()
for rep in range(int(os.environ.get("REPS", 3))):
    u1, u2 = frame(A)
    cens = np.array([section(c, A, h, u1, u2)[0] for h in np.arange(float(os.environ.get('FIT_LO', -10)), float(os.environ.get('FIT_HI', 14)) + 0.1, 1)])
    m = cens.mean(0); _, _, Vt = np.linalg.svd(cens - m); A2 = Vt[0] if Vt[0] @ A > 0 else -Vt[0]
    c = m + A2 * ((c - m) @ A2); print('rep', rep, 'axis change deg', round(float(np.degrees(np.arccos(min(1, A2 @ A)))), 2)); A = A2
u1, u2 = frame(A)
HS = np.arange(-15, 40.1, 0.5)
prof = []
for h in HS:
    _, r = section(c, A, h, u1, u2); F = np.fft.rfft(r) / len(r)
    prof.append((h, float(r.mean()), float(2 * abs(F[3])), float(2 * abs(F[2]))))
prof = np.array(prof)
for row in prof[::2]: print(f'h {row[0]:6.1f}  meanR {row[1]:5.2f}  H3 {row[2]:4.2f}  H2 {row[3]:4.2f}')
json.dump(dict(c=c.tolist(), A=A.tolist(), u1=u1.tolist(), u2=u2.tolist(), thr=thr, lum=lum), open(D + '/frame.json', 'w'))
# diagnostics: long-axis views every 30 degrees and cross-sections
def plane_img(cen, e1, e2, half=36, px=0.25, W_=900, L_=250):
    t = np.arange(-half, half, px); a_, b_ = np.meshgrid(t, t, indexing='xy')
    Pp = cen[None, None] + a_[..., None] * e1 + (-b_)[..., None] * e2
    v = map_coordinates(W, [Pp[..., 2] / 1.5 - .5, Pp[..., 1] / 1.5 - .5, Pp[..., 0] / 1.5 - .5], order=3, cval=-1000)
    return Image.fromarray(np.clip((v - (L_ - W_ / 2)) / W_ * 255, 0, 255).astype(np.uint8))
tiles = []
for ang in range(0, 180, 30):
    e1 = np.cos(np.radians(ang)) * u1 + np.sin(np.radians(ang)) * u2
    im = plane_img(c + A * 8, e1, A).convert('RGB'); d = ImageDraw.Draw(im)
    y0 = 36 / 0.25 + 8 / 0.25                                     # h = 0 line (A points up in the image)
    for hh, col in ((0, (0, 255, 255)),):
        y = 36 / 0.25 - (hh - 8) / 0.25; d.line([(0, y), (288, y)], fill=col)
    for hh in range(-10, 31, 5):
        y = 36 / 0.25 - (hh - 8) / 0.25; d.text((2, y - 6), str(hh), fill=(255, 255, 0))
    d.text((200, 4), f'{ang} deg', fill=(255, 255, 0)); tiles.append(im)
row = Image.new('RGB', (288 * 3 + 6, 288 * 2 + 3))
for i, im in enumerate(tiles): row.paste(im, ((i % 3) * 291, (i // 3) * 291))
row.save(D + '/diag_long.png')
cs = []
for hh in range(-6, 22, 2):
    im = plane_img(c + A * hh, u2, u1, half=20, px=0.2).convert('RGB'); ImageDraw.Draw(im).text((3, 3), f'h {hh}', fill=(255, 255, 0)); cs.append(im)
g = Image.new('RGB', (200 * 7 + 12, 200 * 2 + 2))
for i, im in enumerate(cs): g.paste(im, ((i % 7) * 202, (i // 7) * 202))
g.save(D + '/diag_cross.png')
