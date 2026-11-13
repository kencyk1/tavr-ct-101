# Long-axis views through each cusp centre (radial direction to the right) to read the leaflet hinge (nadir) by eye.
import numpy as np, json, os, sys
from scipy.ndimage import map_coordinates
from PIL import Image, ImageDraw
D = os.environ['CASE_DIR']; offs = [float(x) for x in os.environ.get('OFFS', '-15,0,15').split(',')]
W = np.load(D + '/W.npy').astype(np.float32); F = json.load(open(D + '/' + os.environ.get('FRAME', 'frame_ann.json')))
c = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2']); cus = json.load(open(D + '/cusps.json'))['cus']
L_, W_ = float(os.environ.get('L', 250)), float(os.environ.get('WW', 900))
s = lambda P: map_coordinates(W, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=3, cval=-1000)
px = 0.1; xs = np.arange(-6, 22, px); hs = np.arange(-10, 18, px)[::-1]; X_, H_ = np.meshgrid(xs, hs, indexing='xy')
tiles = []
for cd in cus:
    for o in offs:
        th = np.radians(cd + o); e = np.cos(th) * u1 + np.sin(th) * u2
        v = s(c[None, None] + X_[..., None] * e + H_[..., None] * A)
        im = Image.fromarray(np.clip((v - (L_ - W_ / 2)) / W_ * 255, 0, 255).astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(im)
        for hh in range(-10, 18, 2):
            y = (17.9 - hh) / px; d.line([(0, y), (8 if hh % 10 else 16, y)], fill=(0, 255, 255)); d.text((18, y - 6), str(hh), fill=(0, 255, 255))
        for rr in range(0, 22, 2):
            x = (rr + 6) / px; d.line([(x, im.height - 1), (x, im.height - (8 if rr % 10 else 16))], fill=(255, 160, 0))
            if rr % 4 == 0: d.text((x - 4, im.height - 28), str(rr), fill=(255, 160, 0))
        d.text((60, 3), f'cusp {cd:.0f}{o:+.0f} deg', fill=(255, 255, 0)); tiles.append(im)
n = len(offs); g = Image.new('RGB', (280 * n + 4 * (n - 1), 280 * len(cus) + 4 * (len(cus) - 1)), (40, 40, 40))
for i, im in enumerate(tiles): g.paste(im, ((i % n) * 284, (i // n) * 284))
g.save(D + '/' + os.environ.get('OUT', 'hinge_views.png')); print(g.size)
