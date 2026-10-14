// Web Worker: runs anarcism numbering off the main thread.

import { Anarcism } from '../lib/anarcism/index.js';

let enginePromise = null;

self.onmessage = async (event) => {
  const { seqs, options } = event.data || {};
  try {
    if (!enginePromise) enginePromise = Anarcism.create();
    const engine = await enginePromise;
    const results = engine.numberSequences(seqs, options);
    self.postMessage({ ok: true, results });
  } catch (err) {
    self.postMessage({ ok: false, error: err && err.message ? err.message : String(err) });
  }
};
