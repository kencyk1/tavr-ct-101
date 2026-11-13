# Link-preview image and icons

`og.jpg` (1200 × 630, used for social link previews) and the icons in `icons/` are made from the app itself.

1. `npm i @fontsource-variable/bricolage-grotesque@5 @fontsource/ibm-plex-sans@5 @fontsource/ibm-plex-mono@5 playwright` in this folder (the fonts are served locally so the captures match the site).
2. Serve the repository root on port 8765 (`python3 -m http.server 8765`) and run `node cap.js`: it captures the body coronal with organ names, gated root 1 on its reference plane with a virtual BEV, and root 3 at the left main.
3. `python3 crop.py`, then `node render.js` (renders `og.html` at 2×; downscale to 1200 × 630 and save as JPEG, quality about 88, to keep it well under 300 KB for messaging apps).
4. `node icons.js` renders `icon.svg` (the LubDubAPP heart with a transcatheter valve frame on the annulus line) to the favicon, and `icon-square.svg` (white square, for iOS) to the apple-touch icon.
