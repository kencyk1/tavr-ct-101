# TAVR CT 101

Hands-on CT basics for structural heart advanced practice providers. Learners work on a real, publicly released CT angiogram, from the neck to the femoral bifurcation, in the browser. Nothing to install.

Companion to [LubDubAPP Guide](https://lubdubappguide.com/).

## Sections

| Section | What learners do |
|---|---|
| **Scroll the scan** | Scroll a real CTA in axial, coronal and sagittal planes. Change the window, probe Hounsfield units, measure vessels, and follow 16 TAVR-focused teaching stops from the carotids to the femoral bifurcation. About 80 organs, vessels, bones and muscles are named on the image, including all four heart chambers, the RVOT, the interventricular septum and the tricuspid, pulmonary and mitral valve planes. A self-test mode hides the names. |
| **MPR lab** | Three linked, rotatable planes, the way structural imagers work. Exercises: find the annular plane, measure the annulus, read the C-arm perpendicularity curve, and measure coronary heights. Each exercise is graded against known answers. |
| **How CT works** | Hounsfield units, windowing, calcium blooming, image orientation, and the TAVR CT protocol. |
| **Cardiac CT** | ECG gating, double-oblique alignment, annulus area/perimeter calculations, and coronary obstruction risk. |
| **Access check** | Sheath-to-femoral artery ratio (SFAR) with calcium and tortuosity modifiers. |
| **Quiz** | Nine questions drawn from the other sections. |

## Run it locally

The page loads its CT data with `fetch`, so it needs a web server rather than a double-click:

```bash
python3 -m http.server 8000
# open http://localhost:8000
```

## Repository layout

```
index.html                 the whole app (HTML, CSS and JavaScript in one file)
data/ct_00..08.png         body CTA, 517 slices at 1.5 mm, plus organ labels
data/root_00..03.png       gated root 1 (240³ at 0.5 mm), plus labels
data/rootfb_00..02.png     the same patient's 1.5 mm scan around that box
data/root2_*.png, root2fb_* gated root 2, same layout
data/root3_*.png, root3fb_* root 3 (valve-in-valve), same layout
tools/organ_labels/        scripts that produced the organ and heart chamber labels
tools/root_reconstruction/ scripts for the earlier reconstructed root case
tools/gated_root/          scripts that produced the gated root case
```

The PNG files are lossless containers, not pictures. In each pixel, red and green hold the high and low bytes of (HU + 1024), and blue holds a label ID.

## About the images

- **Body CTA:** case ct1349 from the TAVRP-PL dataset, which is derived from the TotalSegmentator dataset. It was cropped, re-oriented and compressed for the web. It is probably not ECG-gated.
- **Organ labels:** automatic segmentation with TotalSegmentator. The pulmonary trunk was relabeled from "heart" to "pulmonary artery".
- **Heart chambers:** RA, RV, LA, LV and LV myocardium come from the CADS open model (task 553), run on this scan. The RVOT (RV within 14 mm of the pulmonary artery), the interventricular septum and the tricuspid, mitral and pulmonary valve planes are derived from where those chambers meet, so they are approximate. Leftover heart wall of muscle density (mostly the thin RV free wall) is labeled myocardium. Labels may be imperfect at organ edges.
- **Gated roots (MPR lab):** cases s1230 and s1105 of the TotalSegmentator dataset v2, ECG-gated TAVR planning CTs of two other patients (released at 1.5 mm). For each, a 12 cm box around the aortic root is resampled to 0.5 mm with cubic interpolation, light sharpening and fine grain; nothing is modelled. The cusps, their calcium and both coronary ostia are the patients' own, and each patient's 1.5 mm scan fills a 26 cm box around it. Labels come from TotalSegmentator (fast), CADS and lumen geometry; the coronaries are traced from the image. Each set of reference values was measured on its scan. The annulus plane runs through the three leaflet hinge points, read on cross-sections, and sits where the cusps just disappear (the lumen area stops rising below it). The contour follows the half-maximum lumen edge, bridged smoothly across calcium as the SCCT consensus describes. Case 1 (s1230): annulus 31.6 × 23.7 mm, area 589 mm², perimeter 88.7 mm; left main 8.0 mm (low), RCA 20.0 mm (high in the sinus). Case 2 (s1105): annulus 26.1 × 20.2 mm, area 415 mm², perimeter 74.5 mm; bulky calcium at the left-coronary hinge reaching into the LVOT; sinotubular junction about 19 mm up; left main 8.0 mm, RCA 12.0 mm. An earlier version fitted the plane where each sinus starts to bulge; a structural-heart reviewer found it tilted about 10–15° toward the septum, and the hinge-based planes correct that (11.6° and 12.8°). Treat these values as careful automated reads, not core-lab measurements.
- **Root 3 · valve-in-valve (MPR lab):** case s0972 of the same dataset, a cardiac CT (heart and thoracic aorta protocol) of a patient with a stented surgical bioprosthesis in the aortic position, prepared the same way. The valve axis and inflow plane come from the metal frame (inflow = lowest metal). The dataset does not name the valve model, so sizes are CT measurements: frame inner contour 21.6 × 21.3 mm (area 361 mm², perimeter 67.5 mm; metal blooming makes this read slightly small), frame centreline about 23 mm, stent posts 14.5 mm high. Left main lower edge 6.5 mm and RCA 22.5 mm above the inflow; valve-to-coronary distance to a 23 mm virtual valve about 10.2 mm (left main) and 8.6 mm (RCA). The left main was traced from a seed placed by eye on the vessel beside the left sinus, because the frame's blooming confused the automatic tracing. A fifth exercise overlays a schematic virtual valve so learners can measure the valve-to-coronary distance themselves (left main reference: 10.1 mm with BEV 23, 11.1 mm with SEV 23).
- **Virtual valve overlay (MPR lab):** a toggle (BEV / SEV) with a size menu draws a schematic transcatheter valve on the case's reference plane in all three panes. Each case opens at the size its reference annulus fits: BEV by CT annulus area (SAPIEN 3 IFU ranges 20: 273–345, 23: 338–430, 26: 430–546, 29: 540–683 mm²; heights 15.5/18/20/22.5 mm, about 20% ventricular in native roots, 10% in valve-in-valve) and SEV by annulus perimeter (Evolut PRO+ IFU ranges 23: 56.5–62.8, 26: 62.8–72.3, 29: 72.3–81.7, 34: 81.7–94.2 mm; target depth 3 mm). Gated root 1 (area 589 mm², perimeter 88.7 mm) opens at BEV 29 and SEV 34; gated root 2 (415 mm², 74.5 mm) at BEV 23 and SEV 29; root 3 (valve-in-valve, valve model unknown) at 23 mm for both. The BEV outline is a plain cylinder; the SEV outline (inflow = size, waist size − 4 mm, outflow size + 6 mm, 45 mm tall) is a schematic, not a manufacturer drawing. Teaching only: follow the current IFU.

**For education only. Not for clinical decision-making.** Size devices from the manufacturer's instructions for use.

## Credits

Created by [Ken Chan](https://www.linkedin.com/in/ken-chan-6774411a5/), ACNPC-AG, CCK, RCIS, AACC. Built with the help of AI (Claude, Anthropic).

Data and tools:

- TAVRP-PL dataset, Zenodo [10.5281/zenodo.16274177](https://doi.org/10.5281/zenodo.16274177), CC BY 4.0.
- Wasserthal J. et al. TotalSegmentator dataset, Zenodo [10.5281/zenodo.6802613](https://doi.org/10.5281/zenodo.6802613), CC BY 4.0. MPR roots: cases s1230, s1105 and s0972 of dataset v2.0.1, Zenodo [record 10047292](https://zenodo.org/records/10047292), CC BY 4.0.
- Wasserthal J. et al. *TotalSegmentator: robust segmentation of 104 anatomic structures in CT images.* Radiology: Artificial Intelligence 2023. [Software](https://github.com/wasserth/TotalSegmentator), Apache-2.0.
- Xu M. et al. *CADS: a comprehensive anatomical dataset and segmentation for whole-body anatomy in computed tomography.* Nature Communications 2026, doi:[10.1038/s41467-026-78732-1](https://doi.org/10.1038/s41467-026-78732-1). [Software](https://github.com/murong-xu/CADS), Apache-2.0; open model weights CC BY-SA 4.0.

## License

See [LICENSE.md](LICENSE.md). Original content is free for non-commercial educational use with attribution. Third-party data keeps its own licence.
