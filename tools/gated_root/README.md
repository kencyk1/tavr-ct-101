# Gated root case (MPR lab)

Source: case s1230 of the TotalSegmentator dataset v2.0.1 (Zenodo record 10047292, CC BY 4.0), an ECG-gated TAVR planning CT.

Order of the scripts (paths in them point at the original build machine; edit them to rerun):

1. `TotalSegmentator -i ct.nii.gz -o ts_fast.nii.gz --fast --ml` (organs, aorta)
2. `prep.py` puts the CT and labels into viewer orientation (z head to feet, y anterior to posterior, x right to left, 1.5 mm).
3. `cads_any.py ct.nii.gz 250 442 cads.nii.gz` runs the CADS heart model on the chest; `chambers.py 250` merges it.
4. `frame.py` finds the root axis from lumen centroids; `cusps.py` finds the commissures and cusp centres.
5. `nadir.py` fits the annulus plane through the three cusp nadirs (where each sinus starts to bulge above the LVOT).
6. `build_root.py` resamples a 12 cm box to 0.5 mm, labels the lumen, cusps, calcium and coronaries, measures the answer key
   (`root_real_truth.json`) and writes `data/root_*.png` plus the 1.5 mm surroundings `data/rootfb_*.png`.

Each answer was checked visually: the annulus plane on cross-sections and long-axis views, and each coronary ostium on zoomed cross-sections.
