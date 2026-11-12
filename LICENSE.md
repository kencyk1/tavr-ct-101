# License

## Original content

The page, its text, figures, exercises and code (`index.html`, `tools/`) are © 2026 Ken Chan. They are licensed under the **Creative Commons Attribution-NonCommercial 4.0 International** license (CC BY-NC 4.0): https://creativecommons.org/licenses/by-nc/4.0/

You may share and adapt this material for non-commercial purposes. You must give appropriate credit, link to the license, and indicate whether you made changes.

## Third-party material

- **CT images** (`data/ct_*.png`, and the image data underlying `data/root_*.png`) are derived from case ct1349 of the TAVRP-PL dataset (https://doi.org/10.5281/zenodo.16274177). That dataset is itself derived from the TotalSegmentator dataset (Wasserthal et al., https://doi.org/10.5281/zenodo.6802613). Both are licensed under **CC BY 4.0**: https://creativecommons.org/licenses/by/4.0/ The images here were cropped, re-oriented, resampled and compressed. The root case also contains modelled cusps, coronary arteries and wall borders.
- **Organ labels** were generated with the TotalSegmentator software (https://github.com/wasserth/TotalSegmentator), which is licensed under **Apache-2.0**.
- **Heart chamber labels** (right and left atria and ventricles, myocardium, and the RVOT, septum and valve-plane labels derived from them) were generated with the CADS open model weights (https://github.com/murong-xu/CADS), which are licensed under **CC BY-SA 4.0**: https://creativecommons.org/licenses/by-sa/4.0/ To respect that ShareAlike term, this label layer (the blue channel of the PNG files, for those label IDs) is shared under CC BY-SA 4.0. Citation: Xu M. et al., *Nature Communications* 2026, doi:10.1038/s41467-026-78732-1.

## Disclaimer

For education only. Not for clinical decision-making. Size devices from the manufacturer's instructions for use.
