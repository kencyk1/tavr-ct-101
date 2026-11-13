const fs = require('fs'), path = require('path');
const F = n => path.join(__dirname, 'node_modules', n);
const files = {
  'brico.woff2': F('@fontsource-variable/bricolage-grotesque/files/bricolage-grotesque-latin-wght-normal.woff2'),
  'plex400.woff2': F('@fontsource/ibm-plex-sans/files/ibm-plex-sans-latin-400-normal.woff2'),
  'plex500.woff2': F('@fontsource/ibm-plex-sans/files/ibm-plex-sans-latin-500-normal.woff2'),
  'plex600.woff2': F('@fontsource/ibm-plex-sans/files/ibm-plex-sans-latin-600-normal.woff2'),
  'plex700.woff2': F('@fontsource/ibm-plex-sans/files/ibm-plex-sans-latin-700-normal.woff2'),
  'mono400.woff2': F('@fontsource/ibm-plex-mono/files/ibm-plex-mono-latin-400-normal.woff2'),
  'mono500.woff2': F('@fontsource/ibm-plex-mono/files/ibm-plex-mono-latin-500-normal.woff2'),
};
const css = `
@font-face{font-family:"Bricolage Grotesque";src:url(https://fonts.gstatic.com/local/brico.woff2) format("woff2");font-weight:200 800}
@font-face{font-family:"IBM Plex Sans";src:url(https://fonts.gstatic.com/local/plex400.woff2) format("woff2");font-weight:400}
@font-face{font-family:"IBM Plex Sans";src:url(https://fonts.gstatic.com/local/plex500.woff2) format("woff2");font-weight:500}
@font-face{font-family:"IBM Plex Sans";src:url(https://fonts.gstatic.com/local/plex600.woff2) format("woff2");font-weight:600}
@font-face{font-family:"IBM Plex Sans";src:url(https://fonts.gstatic.com/local/plex700.woff2) format("woff2");font-weight:700}
@font-face{font-family:"IBM Plex Mono";src:url(https://fonts.gstatic.com/local/mono400.woff2) format("woff2");font-weight:400}
@font-face{font-family:"IBM Plex Mono";src:url(https://fonts.gstatic.com/local/mono500.woff2) format("woff2");font-weight:500}`;
module.exports = async ctx => {
  await ctx.route('https://fonts.googleapis.com/**', r => r.fulfill({ status: 200, contentType: 'text/css', body: css }));
  await ctx.route('https://fonts.gstatic.com/local/*', r => { const f = files[r.request().url().split('/').pop()]; r.fulfill({ status: 200, contentType: 'font/woff2', body: fs.readFileSync(f), headers: { 'access-control-allow-origin': '*' } }); });
};
module.exports.css = css; module.exports.files = files;
