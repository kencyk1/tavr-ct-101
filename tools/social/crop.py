# Crop the app captures (from cap.js, 2x) into the three panels used by og.html.
from PIL import Image
c = Image.open('cap-coronal.png').convert('RGB'); c.crop((655, 0, 1214, 1206)).save('img-body.png')   # body coronal, centred on the trunk
m = Image.open('cap-mpr1.png').convert('RGB'); m.crop((10, 10, 668, 668)).save('img-ann.png')         # gated root 1: axial pane on the annulus, virtual BEV
v = Image.open('cap-mpr3.png').convert('RGB'); v.crop((678, 10, 1336, 580)).save('img-viv.png')       # root 3 (valve-in-valve): coronal pane at the left main
