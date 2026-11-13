# Access reconstruction (Access check tab)

Centerlines and reference lumen diameters for the iliofemoral reconstruction on the body CTA (case ct1349, 1.5 mm).

1. `centerline.py` traces each side slice by slice from the distal aorta (label 1) along the iliac labels (6 = left, 7 = right), then follows the contrast lumen below the end of the labels to just past the femoral bifurcation; smoothed and resampled every 1 mm (world mm).
2. `diam.py` cuts a 30 × 30 mm plane square to the path at every step (0.2 mm pixels, cubic interpolation), segments the lumen at half maximum with calcium excluded, and records the area-derived, smallest and largest width.
3. Landmarks (aortic and iliac bifurcations, top of the femoral head, femoral bifurcation), the calcium flags and the narrowest point between the bifurcations are packed into the `RECON` object in `index.html`. Reference minimal lumen: right 8.7 mm (external iliac), left 9.2 mm (common femoral). Expect about ±1 mm at 1.5 mm voxels.
