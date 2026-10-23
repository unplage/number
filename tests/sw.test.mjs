import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import vm from 'node:vm';

const ORIGIN = 'https://example.github.io';
const SCOPE = `${ORIGIN}/number/`;

async function loadWorker({ keys = [], fetchRejects = false } = {}) {
  const source = await readFile(new URL('../sw.js', import.meta.url), 'utf8');
  const handlers = new Map();
  const deleted = [];
  const openCalls = [];
  const fetches = [];
  const stores = new Map();

  const cacheFor = (name) => {
    if (!stores.has(name)) {
      stores.set(name, {
        name,
        entries: new Map(),
        puts: [],
        async match(req) {
          const key = typeof req === 'string' ? new URL(req, `${ORIGIN}/number/sw.js`).href : req.url;
          return this.entries.get(key);
        },
        async put(req, res) {
          this.puts.push(req.url);
          this.entries.set(req.url, res);
        },
        async addAll() {},
      });
    }
    return stores.get(name);
  };

  const context = vm.createContext({
    URL,
    Promise,
    console,
    self: {
      location: { pathname: '/number/sw.js', origin: ORIGIN },
      registration: { scope: SCOPE },
      clients: { claim: async () => {} },
      skipWaiting: async () => {},
      addEventListener: (type, fn) => handlers.set(type, fn),
    },
    caches: {
      keys: async () => keys,
      delete: async (name) => {
        deleted.push(name);
        stores.delete(name);
        return true;
      },
      open: async (name) => {
        openCalls.push(name);
        return cacheFor(name);
      },
    },
    fetch: async (req) => {
      fetches.push(req.url);
      if (fetchRejects) throw new TypeError('offline');
      return { ok: true, url: req.url, clone() { return this; } };
    },
  });
  vm.runInContext(`${source}\n;globalThis.__cacheName = CACHE;`, context);

  return { handlers, deleted, openCalls, fetches, stores, cacheFor, context };
}

test('activate deletes only this project\'s old caches', async () => {
  const { handlers, deleted } = await loadWorker({
    keys: [
      'pwa-cache-other-project--v1',
      'pwa-cache-number--v7',
      'pwa-cache-number--v8',
      'pwa-cache-number--v9',
      'antibody-number-v7',
      'unrelated-cache',
    ],
  });

  let promise;
  handlers.get('activate')({ waitUntil: (p) => { promise = p; } });
  await promise;

  assert.deepEqual(deleted.sort(), ['antibody-number-v7', 'pwa-cache-number--v7', 'pwa-cache-number--v8']);
  assert.ok(!deleted.includes('pwa-cache-other-project--v1'));
  assert.ok(!deleted.includes('pwa-cache-number--v9'));
  assert.ok(!deleted.includes('unrelated-cache'));
});

test('fetch only serves from this project\'s cache and only inside scope', async () => {
  const { handlers, stores, fetches, cacheFor, openCalls } = await loadWorker();
  const fetchHandler = handlers.get('fetch');
  const current = cacheFor('pwa-cache-number--v9');
  const other = cacheFor('pwa-cache-other-project--v1');
  other.entries.set(`${SCOPE}app.js`, 'other-pwa-copy');

  let responded;
  fetchHandler({
    request: { method: 'GET', url: `${SCOPE}app.js` },
    respondWith: (p) => { responded = p; },
  });
  assert.ok(responded, 'in-scope requests must be handled');
  const res = await responded;
  assert.notEqual(res, 'other-pwa-copy', 'must not read another PWA\'s cache');
  assert.deepEqual(fetches, [`${SCOPE}app.js`]);
  assert.ok(current.puts.includes(`${SCOPE}app.js`), 'network response cached locally');
  assert.ok(openCalls.every((name) => name === 'pwa-cache-number--v9'));

  let outside;
  fetchHandler({
    request: { method: 'GET', url: `${ORIGIN}/other/index.html` },
    respondWith: (p) => { outside = p; },
  });
  assert.equal(outside, undefined, 'out-of-scope requests must not be intercepted');
});

test('offline navigation falls back to this project\'s cached index', async () => {
  const { handlers, cacheFor } = await loadWorker({ fetchRejects: true });
  const fetchHandler = handlers.get('fetch');
  const current = cacheFor('pwa-cache-number--v9');
  current.entries.set(`${SCOPE}index.html`, 'cached-shell');

  let responded;
  fetchHandler({
    request: { method: 'GET', url: `${SCOPE}missing.html` },
    respondWith: (p) => { responded = p; },
  });
  assert.equal(await responded, 'cached-shell');
});
