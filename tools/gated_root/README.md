# Gated root case (MPR lab)

Sources: cases s1230 (gated root 1) and s1105 (gated root 2) of the TotalSegmentator dataset v2.0.1 (Zenodo record 10047292, CC BY 4.0), ECG-gated TAVR planning CTs.
Set `CASE_DIR` to the case folder and `PREFIX` to `root` or `root2`.

Order of the scripts (paths in them point at the original build machine; edit them to rerun):

1. `TotalSegmentator -i ct.nii.gz -o ts_fast.nii.gz --fast --ml` (organs, aorta)
2. `prep.py` puts the CT and labels into viewer orientation (z head to feet, y anterior to posterior, x right to left, 1.5 mm).
3. `cads_any.py ct.nii.gz 250 442 cads.nii.gz` runs the CADS heart model on the chest; `chambers.py 250` merges it.
4. `frame.py` finds the root axis from lumen centroids (s1105 used `FIT_LO=-2 FIT_HI=22 REPS=6`, because its LVOT is short); `cusps.py` finds the commissures and cusp centres.
5. `nadir.py` fits the annulus plane through the three cusp nadirs (where each sinus starts to bulge above the LVOT); s1105 used a damped fit, `NIT=10 DAMP=0.4`.
6. `build_root.py` resamples a 12 cm box to 0.5 mm, labels the lumen, cusps, calcium and coronaries, measures the answer key
   and writes the 0.5 mm strips plus the 1.5 mm surroundings.
7. `ann_truth.py write` re-measures the annulus with the SCCT rule (contour bridged smoothly across calcium) and stores the
   answer key in `root_real_truth_<case>.json`.

Each answer was checked visually: the annulus plane on cross-sections and long-axis views, and each coronary ostium on zoomed cross-sections.

Note: gated root 1 (s1230) was built with an earlier copy of `build_root.py` that searched for the RCA ostium between 10 and 40 mm and took the left and right cusps from the list order; the current script picks the cusps by angle and searches from 4 mm. Its stored answer key is unchanged.
