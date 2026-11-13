# Cross-sections perpendicular to the frame axis at a list of heights, cusp labels on the rim (for reading hinge levels).
import numpy as np, json, os
from scipy.ndimage import map_coordinates
from PIL import Image, ImageDraw
D = os.environ['CASE_DIR']; F = json.load(open(D + '/' + os.environ.get('FRAME', 'frame_ann.json'))); hs = [float(x) for x in os.environ['HS'].split(',')]
half = float(os.environ.get('HALF', 18)); L_, W_ = float(os.environ.get('L', 250)), float(os.environ.get('WW', 700))
W = np.load(D + '/W.npy').astype(np.float32); c = np.array(F['c']); A = np.array(F['A']); u1 = np.array(F['u1']); u2 = np.array(F['u2'])
cus = json.load(open(D + '/cusps.json'))['cus']; names = os.environ.get('NAMES', 'LCC,NCC,RCC').split(',')
s = lambda P: map_coordinates(W, [P[..., 2] / 1.5 - .5, P[..., 1] / 1.5 - .5, P[..., 0] / 1.5 - .5], order=3, cval=-1000)
px = half / 150; t = np.arange(-half, half, px); a_, b_ = np.meshgrid(t, t, indexing='xy'); ims = []
for hh in hs:
    v = s((c + A * hh)[None, None] + a_[..., None] * u2 + (-b_)[..., None] * u1)
    im = Image.fromarray(np.clip((v - (L_ - W_ / 2)) / W_ * 255, 0, 255).astype(np.uint8)).convert('RGB'); d = ImageDraw.Draw(im); d.text((3, 3), f'h {hh:g}', fill=(255, 255, 0))
    for cd, nm in zip(cus, names):
        x = 150 + 138 * np.sin(np.radians(cd)); y = 150 - 138 * np.cos(np.radians(cd)); d.text((x - 10, y - 6), nm, fill=(255, 120, 0))
    ims.append(im)
n = int(os.environ.get('NCOL', 6)); g = Image.new('RGB', (302 * n, 302 * ((len(ims) + n - 1) // n)))
for i, im in enumerate(ims): g.paste(im, ((i % n) * 302, (i // n) * 302))
g.save(D + '/' + os.environ.get('OUT', 'xsec.png')); print(g.size)
