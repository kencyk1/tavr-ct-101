const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const p = await (await b.newContext({ viewport: { width: 1200, height: 630 }, deviceScaleFactor: 2 })).newPage();
  await p.goto('file://' + __dirname + '/og.html'); await p.evaluate(() => document.fonts.ready); await p.waitForTimeout(400);
  await p.screenshot({ path: 'og2x.png' });
  const q = await (await b.newContext({ viewport: { width: 512, height: 512 }, deviceScaleFactor: 1 })).newPage();
  for (const s of [512, 192, 180, 32]) { await q.setViewportSize({ width: s, height: s }); await q.setContent(`<html><body style="margin:0;background:transparent"><img src="file://${__dirname}/icon.svg" style="width:${s}px;height:${s}px;display:block"></body></html>`); await q.waitForTimeout(150); await q.screenshot({ path: `icon-${s}.png`, omitBackground: true }); }
  await b.close();
})();
