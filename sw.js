/* TAVR CT 101: service worker for the installable app.
   Network first: while online every request goes to the network, so a new version shows up on the next load exactly
   as it would without this file. Each successful same-origin response is also saved on the device, so once the page
   and a case's images have loaded they open without internet. Each deploy changes VERSION, which refreshes the saved
   copy of the core files and removes old caches. */
var VERSION = 'd12258955b62', CACHE = 'tavrct-site', FONTS = 'tavrct-fonts';
var ASSETS = ["./", "manifest.webmanifest", "icons/icon-192.png", "icons/icon-512.png", "icons/apple-touch-icon.png", "icons/favicon.svg"];
var PAGE = new URL('./', self.registration.scope).href;

self.addEventListener('install', function (e) {
  e.waitUntil(caches.open(CACHE).then(function (c) {
    return Promise.all(ASSETS.map(function (u) {
      var url = new URL(u, self.registration.scope).href;
      return fetch(new Request(url, { cache: 'reload' })).then(function (r) { if (r.ok && !r.redirected) return c.put(url, r); }).catch(function () {});
    }));
  }).then(function () { return self.skipWaiting(); }));
});

self.addEventListener('activate', function (e) {
  e.waitUntil(caches.keys().then(function (ks) {
    return Promise.all(ks.filter(function (k) { return k.indexOf('tavrct-') === 0 && k !== CACHE && k !== FONTS; }).map(function (k) { return caches.delete(k); }));
  }).then(function () { return self.clients.claim(); }));
});

function keep(c, key, r) {
  if (!r || !r.ok || r.type !== 'basic' || r.redirected) return Promise.resolve();
  return c.match(key).then(function (old) {
    var a = old && old.headers.get('etag'), b = r.headers.get('etag');
    if (a && b && a === b) return;
    return c.put(key, r);
  });
}

self.addEventListener('fetch', function (e) {
  var req = e.request; if (req.method !== 'GET') return;
  var url = new URL(req.url);
  if (url.origin === self.location.origin) {
    if (/\/sw\.js$/.test(url.pathname)) return;
    var page = req.mode === 'navigate', key = page && /\/(index\.html)?$/.test(url.pathname) ? PAGE : url.origin + url.pathname;
    e.respondWith(caches.open(CACHE).then(function (c) {
      return fetch(req).then(function (r) { e.waitUntil(keep(c, key, r.clone())); return r; }).catch(function () {
        return c.match(key).then(function (hit) {
          if (hit) return hit;
          return page ? c.match(PAGE).then(function (p) { return p || Response.error(); }) : Response.error();
        });
      });
    }));
    return;
  }
  if (url.hostname === 'fonts.googleapis.com' || url.hostname === 'fonts.gstatic.com') {
    e.respondWith(caches.open(FONTS).then(function (c) {
      return c.match(req).then(function (hit) {
        var net = fetch(req).then(function (r) { if (r.ok || r.type === 'opaque') c.put(req, r.clone()); return r; }).catch(function () { return hit || Response.error(); });
        return hit || net;
      });
    }));
  }
});
