// Analysis orchestrator: chain typing + IMGT via anarcism, Kabat via ports.

import { Anarcism } from '../lib/anarcism/index.js';
import { kabatHeavyFromDomain } from './anarci-port.js';
import { numberLight } from './kabat-light.js';
import { isVhh } from './vhh.js';

let enginePromise = null;

export function loadEngine(wasmSource) {
  if (!enginePromise) {
    enginePromise = Anarcism.create(wasmSource ? { source: wasmSource } : undefined);
  }
  return enginePromise;
}

const CHAIN_LABELS = {
  H: 'Heavy (IGH)',
  K: 'Kappa (IGK)',
  L: 'Lambda (IGL)',
  A: 'TCR alpha',
  B: 'TCR beta',
  G: 'TCR gamma',
  D: 'TCR delta',
};

const KABAT_CDR_RANGES = {
  H: { CDR1: [31, 35], CDR2: [50, 65], CDR3: [95, 102] },
  L: { CDR1: [24, 34], CDR2: [50, 56], CDR3: [89, 97] },
  K: { CDR1: [24, 34], CDR2: [50, 56], CDR3: [89, 97] },
};

const IMGT_CDR_RANGES = {
  H: { CDR1: [27, 38], CDR2: [56, 65], CDR3: [105, 117] },
  K: { CDR1: [27, 38], CDR2: [56, 65], CDR3: [105, 117] },
  L: { CDR1: [27, 38], CDR2: [56, 65], CDR3: [105, 117] },
};

function regionForPosition(pos, scheme, chainType) {
  const ranges = (scheme === 'imgt' ? IMGT_CDR_RANGES : KABAT_CDR_RANGES)[chainType];
  if (!ranges) return '';
  const m = /^([A-Z]+)(\d+)/.exec(pos);
  if (!m) return '';
  const n = Number(m[2]);
  for (const [name, [a, b]] of Object.entries(ranges)) {
    if (n >= a && n <= b) return name;
  }
  // FR assignment
  if (scheme === 'imgt') {
    if (n <= 26) return 'FR1';
    if (n <= 55) return 'FR2';
    if (n <= 104) return 'FR3';
    return 'FR4';
  }
  if (chainType === 'H') {
    if (n <= 30) return 'FR1';
    if (n <= 49) return 'FR2';
    if (n <= 94) return 'FR3';
    return 'FR4';
  }
  if (n <= 23) return 'FR1';
  if (n <= 49) return 'FR2';
  if (n <= 88) return 'FR3';
  return 'FR4';
}

function extractCdrs(residues) {
  const out = {};
  for (const name of ['CDR1', 'CDR2', 'CDR3']) {
    out[name] = residues.filter((r) => r.region === name).map((r) => r.aa).join('');
  }
  return out;
}

/**
 * Analyze sequences with the chosen scheme.
 * @param {{ id: string, sequence: string }[]} records
 * @param {'imgt'|'kabat'} scheme
 * @returns {Promise<object[]>}
 */
export async function analyze(records, scheme, wasmSource) {
  const engine = await loadEngine(wasmSource);
  const imgtResults = engine.numberSequences(records, { allowedChains: ['H', 'K', 'L'] });

  const results = [];
  for (let i = 0; i < records.length; i++) {
    const record = records[i];
    const imgt = imgtResults[i];
    const base = {
      id: record.id,
      description: record.description,
      sequence: record.sequence,
      scheme,
      warnings: [...(imgt.warnings || [])],
      errors: [],
    };

    if (!imgt.domains || imgt.domains.length === 0) {
      base.errors.push('No antibody variable domain detected (chain typing failed).');
      results.push(base);
      continue;
    }

    base.domains = [];
    for (const domain of imgt.domains) {
      const chainType = domain.chainType;
      const domainResult = {
        chainType,
        chainLabel: CHAIN_LABELS[chainType] || chainType,
        species: domain.species,
        start: domain.start,
        end: domain.end,
        bitScore: domain.bitScore,
        cdrs: {},
        residues: [],
        isVhh: false,
        vhhPositions: [],
        warnings: [],
      };

      if (scheme === 'imgt') {
        domainResult.residues = domain.numbering.map((r) => ({
          pos: `${chainType === 'H' ? 'H' : 'L'}${r.position}${r.insertionCode || ''}`,
          imgtPos: r.position,
          imgtIns: r.insertionCode || '',
          aa: r.aminoAcid,
          region: r.region,
        }));
        domainResult.cdrs = extractCdrs(domainResult.residues);
      } else {
        // Kabat
        if (chainType === 'H') {
          try {
            domainResult.residues = kabatHeavyFromDomain(domain).map((r) => ({
              pos: r.pos,
              aa: r.aa,
              region: regionForPosition(r.pos, 'kabat', 'H'),
            }));
          } catch (e) {
            domainResult.warnings.push(`Kabat heavy numbering failed: ${e.message}`);
          }
        } else {
          // Light chains: use the anchor-based engine on the raw domain sequence.
          const normalized = imgt.normalizedSequence || record.sequence;
          const domainSeq = normalized.slice(domain.start, domain.end);
          try {
            const { residues, warnings } = numberLight(domainSeq);
            domainResult.residues = residues.map((r) => ({
              pos: r.pos,
              aa: r.aa,
              region: r.region,
            }));
            domainResult.warnings.push(...warnings);
          } catch (e) {
            domainResult.warnings.push(`Kabat light numbering failed: ${e.message}`);
          }
        }
        if (domainResult.residues.length) {
          domainResult.cdrs = extractCdrs(domainResult.residues);
        }
      }

      if (chainType === 'H' && domainResult.residues.length) {
        const vhh = isVhh(domainResult.residues);
        domainResult.isVhh = vhh.isVhh;
        domainResult.vhhPositions = vhh.positions;
      }

      base.domains.push(domainResult);
    }
    results.push(base);
  }
  return results;
}
