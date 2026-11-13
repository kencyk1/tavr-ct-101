const { chromium } = require('playwright'); const useFonts = require('./fonts.js');
(async () => {
  const b = await chromium.launch(); const errs = [];
  const ctx = await b.newContext({ viewport: { width: 1440, height: 900 }, deviceScaleFactor: 2 }); await useFonts(ctx);
  const p = await ctx.newPage(); p.on('pageerror', e => errs.push(e.message));
  await p.goto('http://localhost:8765/', { waitUntil: 'networkidle' }); await p.waitForTimeout(4000);
  await (await p.$('#screen')).screenshot({ path: 'cap-axial.png' });
  await p.click('.viewSeg.np [data-v="co"]'); await p.waitForTimeout(800); await (await p.$('#screen')).screenshot({ path: 'cap-coronal.png' });
  await p.goto('http://localhost:8765/#mpr'); await p.waitForTimeout(5000);
  await p.click('#exCard [data-act="auto"]'); await p.waitForTimeout(800);
  await p.click('#mTHV [data-thv="bev"]'); await p.waitForTimeout(800);
  await (await p.$('#panes')).screenshot({ path: 'cap-mpr1.png' });
  // valve-in-valve: go to the left main with the virtual valve
  await p.selectOption('#mCaseSel', 'root3'); await p.waitForTimeout(5000);
  const n = await p.$$eval('#exList button', bs => bs.length); await p.click(`#exList li:nth-child(${n}) button`); await p.waitForTimeout(400);
  await p.click('#exCard [data-act="goLM"]'); await p.waitForTimeout(1200);
  await (await p.$('#panes')).screenshot({ path: 'cap-mpr3.png' });
  console.log('errors:', errs.join(' | ') || 'none'); await b.close();
})();
