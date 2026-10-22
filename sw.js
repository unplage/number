// @ts-nocheck -- service worker global scope typing differs from the DOM lib.
// Service worker: precache the full app (including the ~1 MB numbering engine)
// so the tool works completely offline after the first visit.

// Derive the project id from this worker's own path (e.g. /number/sw.js ->
// "number") so several PWAs can share one *.github.io origin without ever
// touching each other's caches.
const PROJECT = (() => {
  const dir = self.location.pathname.replace(/\/sw\.js$/, '');
  const segment = dir.split('/').filter(Boolean).pop();
  return segment || 'root';
})();
const CACHE_PREFIX = `pwa-cache-${PROJECT}--`;
const CACHE = `${CACHE_PREFIX}v8`;
// Pre-convention cache name used by earlier releases of this same app.
const LEGACY_PREFIX = 'antibody-number-';

const ASSETS = [
  './',
  './index.html',
  './manifest.webmanifest',
  './css/style.css',
  './js/app.js',
  './js/engine.js',
  './js/schemes.js',
  './js/data/blosum62.js',
  './js/validate.js',
  './js/analysis.js',
  './js/fasta.js',
  './js/kabat-light.js',
  './js/vhh.js',
  './js/export.js',
  './lib/anarcism/index.js',
  './lib/anarcism/worker-pool.js',
  './lib/anarcism/worker.js',
  './lib/anarcism/anarcism.wasm',
  './icons/icon.svg',
  './icons/icon-192.png',
  './icons/icon-512.png',
  './icons/apple-touch-icon.png',
];

self.addEventListener('install', (event) => {
  event.waitUntil(
    caches.open(CACHE).then((cache) => cache.addAll(ASSETS)).then(() => self.skipWaiting()),
  );
});

self.addEventListener('activate', (event) => {
  event.waitUntil(
    caches.keys().then((keys) =>
      Promise.all(
        keys
          .filter((k) => k !== CACHE)
          .filter((k) => k.startsWith(CACHE_PREFIX) || k.startsWith(LEGACY_PREFIX))
          .map((k) => caches.delete(k)),
      ),
    ).then(() => self.clients.claim()),
  );
});

self.addEventListener('fetch', (event) => {
  const req = event.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url);
  if (url.origin !== self.location.origin) return;
  // Stay inside this app's own scope (/number/), even if the worker is ever
  // registered more broadly.
  if (!url.pathname.startsWith(new URL(self.registration.scope).pathname)) return;

  event.respondWith(
    caches.open(CACHE).then((cache) =>
      cache.match(req).then((cached) => {
        if (cached) return cached;
        return fetch(req)
          .then((res) => {
            if (res && res.ok) cache.put(req, res.clone());
            return res;
          })
          .catch(() => cache.match('./index.html'));
      }),
    ),
  );
});
