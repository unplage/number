// @ts-nocheck -- service worker global scope typing differs from the DOM lib.
// Service worker: precache the full app (including the ~1 MB numbering engine)
// so the tool works completely offline after the first visit.

const CACHE = 'antibody-number-v2';

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
  './js/number-worker.js',
  './js/fasta.js',
  './js/kabat-light.js',
  './js/vhh.js',
  './js/export.js',
  './lib/anarcism/index.js',
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
      Promise.all(keys.filter((k) => k !== CACHE).map((k) => caches.delete(k))),
    ).then(() => self.clients.claim()),
  );
});

self.addEventListener('fetch', (event) => {
  const req = event.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url);
  if (url.origin !== self.location.origin) return;

  event.respondWith(
    caches.match(req).then((cached) => {
      if (cached) return cached;
      return fetch(req)
        .then((res) => {
          if (res && res.ok) {
            const copy = res.clone();
            caches.open(CACHE).then((cache) => cache.put(req, copy));
          }
          return res;
        })
        .catch(() => caches.match('./index.html'));
    }),
  );
});
