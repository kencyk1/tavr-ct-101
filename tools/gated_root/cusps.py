import numpy as np, json
from scipy.ndimage import map_coordinates, gaussian_filter1d
from PIL import Image, ImageDraw
W = np.load('W.npy').astype(np.float32); F = json.load(open('frame.json')); c = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2'])
s = lambda P: map_coordinates(W, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=3, cval=-1000)
th = np.radians(np.arange(0, 360, 1.0))
prof = 0
for h in (12, 13, 14, 15):
    for r in (7, 9, 11):
        P = c + A * h + r * (np.cos(th)[:, None] * u1 + np.sin(th)[:, None] * u2); prof = prof + s(P)
prof = gaussian_filter1d(prof, 4, mode='wrap')
# three deepest minima at least 70 deg apart
order = np.argsort(prof); com = []
for i in order:
    if all(min(abs(i - j), 360 - abs(i - j)) > 70 for j in com): com.append(int(i))
    if len(com) == 3: break
com = sorted(com); print('commissures (deg from anterior toward patient left):', com)
cus = [((com[k] + com[(k + 1) % 3] + (360 if k == 2 else 0)) / 2) % 360 for k in range(3)]; print('cusp centres', cus)
json.dump(dict(com=com, cus=cus), open('cusps.json', 'w'))
def plane_img(cen, e1, e2, half=28, px=0.2, W_=900, L_=250):
    t = np.arange(-half, half, px); a_, b_ = np.meshgrid(t, t, indexing='xy')
    return s(cen[None, None] + a_[..., None] * e1 + (-b_)[..., None] * e2)
tiles = []
for cdeg in cus:
    e1 = np.cos(np.radians(cdeg)) * u1 + np.sin(np.radians(cdeg)) * u2
    v = plane_img(c + A * 6, e1, A); im = Image.fromarray(np.clip((v + 200) / 900 * 255, 0, 255).astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(im)
    for hh in range(-6, 21, 2):
        y = (28 - (hh - 6)) / 0.2; d.line([(270, y), (280, y)], fill=(255, 255, 0)); d.text((250, y - 6), str(hh), fill=(255, 255, 0))
    y0 = (28 + 6) / 0.2; d.line([(0, y0), (60, y0)], fill=(0, 255, 255))
    d.text((4, 4), f'cusp at {cdeg:.0f} deg (right side of image)', fill=(255, 255, 0)); tiles.append(im)
g = Image.new('RGB', (280 * 3 + 6, 280)); [g.paste(im, (i * 283, 0)) for i, im in enumerate(tiles)]; g.save('diag_cusps.png')
