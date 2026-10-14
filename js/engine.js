// Analysis orchestrator: chain typing + IMGT via anarcism; other schemes via
// the ANARCI port (js/schemes.js). Kabat light chains use the dedicated anchor
// engine (js/kabat-light.js) to avoid ANARCI/AbRSA VL CDR3 errors.

import { Anarcism } from '../lib/anarcism/index.js';
import { schemeNumberingFromDomain } from './schemes.js';
import { numberLight } from './kabat-light.js';
import { isVhh } from './vhh.js';
import { validateSequence } from './validate.js';
import { cdrLengths, findLiabilities } from './analysis.js';

export const SCHEMES = ['imgt', 'kabat', 'chothia', 'martin', 'aho', 'wolfguy'];
export const ALLOWED_CHAINS = ['H', 'K', 'L', 'A', 'B', 'G', 'D'];
const IG_CHAINS = new Set(['H', 'K', 'L']);

const CHAIN_LABELS = {
  H: 'Heavy (IGH)',
  K: 'Kappa (IGK)',
  L: 'Lambda (IGL)',
  A: 'TCR alpha',
  B: 'TCR beta',
  G: 'TCR gamma',
  D: 'TCR delta',
};

// CDR (FR) boundaries per scheme and chain letter. Values are inclusive
// position-number ranges; regions not listed are framework.
const CDR_RANGES = {
  kabat: {
    H: { CDR1: [31, 35], CDR2: [50, 65], CDR3: [95, 102], FR1: [1, 30], FR2: [36, 49], FR3: [66, 94], FR4: [103, 999] },
    L: { CDR1: [24, 34], CDR2: [50, 56], CDR3: [89, 97], FR1: [1, 23], FR2: [35, 49], FR3: [57, 88], FR4: [98, 999] },
  },
  chothia: {
    H: { CDR1: [26, 32], CDR2: [52, 56], CDR3: [95, 102], FR1: [1, 25], FR2: [33, 51], FR3: [57, 94], FR4: [103, 999] },
    L: { CDR1: [24, 34], CDR2: [50, 56], CDR3: [89, 97], FR1: [1, 23], FR2: [35, 49], FR3: [57, 88], FR4: [98, 999] },
  },
  aho: {
    H: { CDR1: [25, 42], CDR2: [58, 77], CDR3: [107, 138], FR1: [1, 24], FR2: [43, 57], FR3: [78, 106], FR4: [139, 999] },
    L: { CDR1: [25, 42], CDR2: [58, 77], CDR3: [107, 138], FR1: [1, 24], FR2: [43, 57], FR3: [78, 106], FR4: [139, 999] },
  },
  wolfguy: {
    H: { CDR1: [151, 199], CDR2: [251, 299], CDR3: [351, 399], FR1: [101, 150], FR2: [200, 250], FR3: [300, 350], FR4: [400, 999] },
    L: { CDR1: [551, 599], CDR2: [651, 699], CDR3: [751, 799], FR1: [501, 550], FR2: [600, 650], FR3: [700, 750], FR4: [800, 999] },
  },
};
CDR_RANGES.martin = CDR_RANGES.chothia;

function lightLetter(chainType) {
  return 'L';
}

function regionForPosition(scheme, chainType, n) {
  const chain = chainType === 'H' ? 'H' : 'L';
  const table = (CDR_RANGES[scheme] || CDR_RANGES.kabat)[chain];
  if (!table) return '';
  for (const [name, [a, b]] of Object.entries(table)) if (n >= a && n <= b) return name;
  return '';
}

function confidenceFromScore(bitScore) {
  if (typeof bitScore !== 'number') return 'unknown';
  if (bitScore >= 200) return 'high';
  if (bitScore >= 100) return 'medium';
  return 'low';
}

function extractCdrs(residues) {
  const out = {};
  for (const name of ['CDR1', 'CDR2', 'CDR3']) out[name] = residues.filter((r) => r.region === name).map((r) => r.aa).join('');
  return out;
}

function posLabel(chainType, n, ins) {
  const prefix = chainType === 'H' ? 'H' : 'L';
  return `${prefix}${n}${ins && ins !== ' ' ? ins : ''}`;
}

// ---------------------------------------------------------------------------
// Engine loading
// ---------------------------------------------------------------------------

let enginePromise = null;
let engineVersion = null;

export function loadEngine(wasmSource) {
  if (!enginePromise) {
    enginePromise = Anarcism.create(wasmSource ? { source: wasmSource } : undefined).then((e) => {
      engineVersion = e.version;
      return e;
    });
  }
  return enginePromise;
}

export function getEngineVersion() {
  return engineVersion;
}

// ---------------------------------------------------------------------------
// Numbering (worker with sync fallback)
// ---------------------------------------------------------------------------

let worker = null;
function getWorker() {
  if (worker !== null) return worker;
  if (typeof Worker === 'undefined') {
    worker = false;
    return worker;
  }
  try {
    worker = new Worker(new URL('./number-worker.js', import.meta.url), { type: 'module' });
  } catch {
    worker = false;
  }
  return worker;
}

const NUMBER_OPTIONS = { allowedChains: ALLOWED_CHAINS, assignGermline: true };

function numberSync(engine, seqs) {
  return engine.numberSequences(seqs, NUMBER_OPTIONS);
}

async function numberViaWorker(w, seqs) {
  return new Promise((resolve, reject) => {
    const onMessage = (e) => {
      w.removeEventListener('message', onMessage);
      w.removeEventListener('error', onError);
      if (e.data && e.data.ok) resolve(e.data.results);
      else reject(new Error(e.data && e.data.error ? e.data.error : 'worker error'));
    };
    const onError = (err) => {
      w.removeEventListener('message', onMessage);
      w.removeEventListener('error', onError);
      reject(err);
    };
    w.addEventListener('message', onMessage);
    w.addEventListener('error', onError);
    w.postMessage({ seqs, options: NUMBER_OPTIONS });
  });
}

// ---------------------------------------------------------------------------
// Domain processing
// ---------------------------------------------------------------------------

function processDomain(record, domain, scheme, engine) {
  const chainType = domain.chainType;
  const out = {
    chainType,
    chainLabel: CHAIN_LABELS[chainType] || chainType,
    species: domain.species,
    start: domain.start,
    end: domain.end,
    bitScore: domain.bitScore,
    eValue: domain.eValue,
    confidence: confidenceFromScore(domain.bitScore),
    germline: domain.germline || null,
    isVhh: false,
    vhhPositions: [],
    residues: [],
    cdrs: {},
    cdrLengths: {},
    liabilities: [],
    warnings: [],
    errors: [],
  };

  if (scheme === 'imgt') {
    out.residues = domain.numbering.map((r) => ({
      pos: `${chainType === 'H' ? 'H' : 'L'}${r.position}${r.insertionCode || ''}`,
      num: r.position,
      ins: r.insertionCode || '',
      aa: r.aminoAcid,
      region: r.region,
    }));
  } else if (!IG_CHAINS.has(chainType)) {
    out.errors.push(`Scheme "${scheme}" is only defined for antibody heavy/kappa/lambda chains (got ${chainType}). Use IMGT for TCR.`);
    return out;
  } else if (scheme === 'kabat' && chainType !== 'H') {
    // Kabat light chains: dedicated anchor engine.
    const normalized = domain.normalizedSequence || record.sequence;
    const seq = normalized.slice(domain.start, domain.end);
    try {
      const { residues, warnings } = numberLight(seq);
      out.residues = residues.map((r) => ({ pos: r.pos, aa: r.aa, region: r.region }));
      out.warnings.push(...warnings);
    } catch (e) {
      out.errors.push(`Kabat light numbering failed: ${e.message}`);
      return out;
    }
  } else {
    try {
      const numbered = schemeNumberingFromDomain(scheme, chainType, domain);
      out.residues = numbered.map((r) => ({
        pos: posLabel(chainType, r.n, r.ins),
        aa: r.aa,
        region: regionForPosition(scheme, chainType, r.n),
      }));
    } catch (e) {
      out.errors.push(`${scheme} numbering failed: ${e.message}`);
      return out;
    }
  }

  if (out.residues.length === 0) {
    out.errors.push('No residues could be numbered.');
    return out;
  }

  out.cdrs = extractCdrs(out.residues);
  out.cdrLengths = cdrLengths(out.cdrs);
  out.liabilities = findLiabilities(out.residues, chainType);

  if (chainType === 'H') {
    try {
      const kabatResidues = schemeNumberingFromDomain('kabat', 'H', domain).map((r) => ({
        pos: posLabel('H', r.n, r.ins),
        aa: r.aa,
      }));
      const vhh = isVhh(kabatResidues);
      out.isVhh = vhh.isVhh;
      out.vhhPositions = vhh.positions;
    } catch {
      /* VHH detection is best-effort */
    }
  }

  return out;
}

// ---------------------------------------------------------------------------
// Public analysis API
// ---------------------------------------------------------------------------

/**
 * Analyze FASTA records.
 * @param {{id:string, description?:string, sequence:string}[]} records
 * @param {string} scheme
 * @param {object} [opts]
 * @param {(done:number,total:number)=>void} [opts.onProgress]
 * @param {AbortSignal} [opts.signal]
 * @param {ArrayBuffer|Uint8Array} [opts.wasmSource]
 * @returns {Promise<{results:object[], engineVersion:string|null}>}
 */
export async function analyze(records, scheme, opts = {}) {
  const { onProgress, signal, wasmSource } = opts;
  const engine = await loadEngine(wasmSource);
  const prepared = records.map((rec) => ({ rec, validation: validateSequence(rec.sequence) }));

  const indexes = [];
  const seqs = [];
  prepared.forEach((p, i) => {
    if (p.validation.ok) {
      indexes.push(i);
      seqs.push({ id: p.rec.id, sequence: p.validation.sequence });
    }
  });

  const total = seqs.length;
  let done = 0;
  if (onProgress) onProgress(0, total);

  const batchSize = 50;
  const resultByIndex = new Map();
  const w = total > batchSize ? getWorker() : false;

  for (let start = 0; start < seqs.length; start += batchSize) {
    if (signal && signal.aborted) throw signal.reason || new Error('aborted');
    const chunk = seqs.slice(start, start + batchSize);
    let numbered;
    try {
      numbered = w ? await numberViaWorker(w, chunk) : numberSync(engine, chunk);
    } catch (e) {
      // Fall back to the main-thread engine if the worker fails.
      numbered = numberSync(engine, chunk);
      void e;
    }
    chunk.forEach((s, j) => resultByIndex.set(start + j, numbered[j]));
    done += chunk.length;
    if (onProgress) onProgress(done, total);
    await new Promise((r) => setTimeout(r, 0));
  }

  const results = prepared.map((p, i) => {
    const base = {
      id: p.rec.id,
      description: p.rec.description || '',
      sequence: p.validation.sequence,
      scheme,
      errors: [],
      warnings: [...p.validation.warnings],
      domains: [],
    };
    if (!p.validation.ok) {
      base.errors = p.validation.errors.map((x) => x.message);
      return base;
    }
    const localIdx = indexes.indexOf(i);
    const imgt = resultByIndex.get(localIdx);
    if (!imgt || !imgt.domains || imgt.domains.length === 0) {
      base.errors.push('No antibody variable domain detected (chain typing failed). Check the sequence or try IMGT.');
      return base;
    }
    base.warnings.push(...(imgt.warnings || []));
    for (const domain of imgt.domains) {
      domain.normalizedSequence = imgt.normalizedSequence;
      base.domains.push(processDomain({ ...base, sequence: p.validation.sequence }, domain, scheme, engine));
    }
    return base;
  });

  return { results, engineVersion };
}
