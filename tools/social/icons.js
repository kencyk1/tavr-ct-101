const { chromium } = require('playwright'); const fs = require('fs');
(async () => {
  const b = await chromium.launch(); const q = await (await b.newContext({ viewport: { width: 512, height: 512 }, deviceScaleFactor: 1 })).newPage();
  for (const [src, sizes, pre] of [['icon.svg', [512, 32], 'icon-'], ['icon-square.svg', [180, 512], 'sq-']]) {
    const svg = fs.readFileSync(__dirname + '/' + src, 'utf8');
    for (const s of sizes) { await q.setViewportSize({ width: s, height: s }); await q.setContent(`<html><body style="margin:0;background:transparent">${svg.replace('<svg ', `<svg width="${s}" height="${s}" style="display:block" `)}</body></html>`); await q.screenshot({ path: `${pre}${s}.png`, omitBackground: true }); }
  }
  await b.close();
})();
