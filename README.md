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
data/root_00..03.png       aortic root reconstruction (240³ at 0.5 mm), plus labels
tools/organ_labels/        scripts that produced the organ and heart chamber labels
tools/root_reconstruction/ scripts that produced the 0.5 mm root case
```

The PNG files are lossless containers, not pictures. In each pixel, red and green hold the high and low bytes of (HU + 1024), and blue holds a label ID.

## About the images

- **Body CTA:** case ct1349 from the TAVRP-PL dataset, which is derived from the TotalSegmentator dataset. It was cropped, re-oriented and compressed for the web. It is probably not ECG-gated.
- **Organ labels:** automatic segmentation with TotalSegmentator. The pulmonary trunk was relabeled from "heart" to "pulmonary artery".
- **Heart chambers:** RA, RV, LA, LV and LV myocardium come from the CADS open model (task 553), run on this scan. The RVOT (RV within 14 mm of the pulmonary artery), the interventricular septum and the tricuspid, mitral and pulmonary valve planes are derived from where those chambers meet, so they are approximate. Leftover heart wall of muscle density (mostly the thin RV free wall) is labeled myocardium. Labels may be imperfect at organ edges.
- **Root reconstruction (MPR lab):** this patient's real LVOT, aortic root and ascending aorta rebuilt in a 12 cm box at 0.5 mm. The rebuild follows the vessel's own curved centerline, with fine thin-slice grain, and the surrounding septum, atria, RVOT and fat are brought in to meet a crisp lumen edge. The basal ring is an ellipse fitted to the real annulus. The real scan cannot show the cusps, so thickened, calcified severe-AS cusps are modelled on the patient's own sinus walls, with hinge-point and LVOT calcium. The proximal coronaries are routed through the patient's epicardial fat planes. Beyond the box, the original 1.5 mm scan fills the view. This makes every answer key exact.

**For education only. Not for clinical decision-making.** Size devices from the manufacturer's instructions for use.

## Credits

Created by [Ken Chan](https://www.linkedin.com/in/ken-chan-6774411a5/), ACNPC-AG, CCK, RCIS, AACC. Built with the help of AI (Claude, Anthropic).

Data and tools:

- TAVRP-PL dataset, Zenodo [10.5281/zenodo.16274177](https://doi.org/10.5281/zenodo.16274177), CC BY 4.0.
- Wasserthal J. et al. TotalSegmentator dataset, Zenodo [10.5281/zenodo.6802613](https://doi.org/10.5281/zenodo.6802613), CC BY 4.0.
- Wasserthal J. et al. *TotalSegmentator: robust segmentation of 104 anatomic structures in CT images.* Radiology: Artificial Intelligence 2023. [Software](https://github.com/wasserth/TotalSegmentator), Apache-2.0.
- Xu M. et al. *CADS: a comprehensive anatomical dataset and segmentation for whole-body anatomy in computed tomography.* Nature Communications 2026, doi:[10.1038/s41467-026-78732-1](https://doi.org/10.1038/s41467-026-78732-1). [Software](https://github.com/murong-xu/CADS), Apache-2.0; open model weights CC BY-SA 4.0.

## License

See [LICENSE.md](LICENSE.md). Original content is free for non-commercial educational use with attribution. Third-party data keeps its own licence.
