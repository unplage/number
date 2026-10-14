import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { analyze, loadEngine, getEngineVersion, ALLOWED_CHAINS, SCHEMES } from '../js/engine.js';

const wasm = readFileSync(new URL('../lib/anarcism/anarcism.wasm', import.meta.url));
await loadEngine(wasm);

const VH = 'EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS';
const VL = 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK';

test('exposes the six schemes and all chain types', () => {
  assert.deepEqual(SCHEMES, ['imgt', 'kabat', 'chothia', 'martin', 'aho', 'wolfguy']);
  assert.ok(ALLOWED_CHAINS.includes('H') && ALLOWED_CHAINS.includes('A'));
});

test('IMGT: detects VH and annotates CDRs, germline and confidence', async () => {
  const { results, engineVersion } = await analyze([{ id: 'vh', sequence: VH }], 'imgt', { wasmSource: wasm });
  const r = results[0];
  assert.equal(r.errors.length, 0);
  const d = r.domains[0];
  assert.equal(d.chainType, 'H');
  assert.ok(d.residues.length > 100);
  assert.ok(d.cdrs.CDR3.length > 0);
  assert.ok(['high', 'medium', 'low'].includes(d.confidence));
  assert.ok(d.germline && d.germline.vGene);
  assert.ok(engineVersion);
});

test('IMGT: kappa light chain and lambda light chain', async () => {
  const [r] = (await analyze([{ id: 'vl', sequence: VL }], 'imgt', { wasmSource: wasm })).results;
  assert.equal(r.domains[0].chainType, 'K');
  assert.ok(r.domains[0].cdrs.CDR1.length > 0);
});

test('Kabat: VH uses ANARCI rules, VL uses anchor engine', async () => {
  const { results } = await analyze(
    [
      { id: 'vh', sequence: VH },
      { id: 'vl', sequence: VL },
    ],
    'kabat',
    { wasmSource: wasm },
  );
  assert.equal(results[0].domains[0].chainType, 'H');
  assert.equal(results[1].domains[0].chainType, 'K');
  assert.equal(results[1].domains[0].cdrs.CDR1, 'RASQDVNTAVA');
  assert.equal(results[1].domains[0].cdrs.CDR3, 'QQHYTTPPT');
  assert.ok(results[1].domains[0].residues.some((x) => x.pos === 'L23'));
});

test('all non-IMGT schemes produce numbering for VH', async () => {
  for (const scheme of ['kabat', 'chothia', 'martin', 'aho', 'wolfguy']) {
    const { results } = await analyze([{ id: 'vh', sequence: VH }], scheme, { wasmSource: wasm });
    const d = results[0].domains[0];
    assert.equal(d.errors.length, 0, scheme);
    assert.ok(d.residues.length > 100, scheme);
    assert.ok(d.cdrs.CDR1 && d.cdrs.CDR2 && d.cdrs.CDR3, scheme);
  }
});

test('computes liabilities and CDR lengths', async () => {
  const { results } = await analyze([{ id: 'vl', sequence: VL }], 'kabat', { wasmSource: wasm });
  const d = results[0].domains[0];
  assert.deepEqual(d.cdrLengths, { CDR1: 11, CDR2: 7, CDR3: 9 });
  assert.ok(Array.isArray(d.liabilities));
});

test('rejects nucleotide input', async () => {
  const { results } = await analyze([{ id: 'dna', sequence: 'ATGGCTGATATTGCCGTGTTAAAGCCGTAA' }], 'imgt', { wasmSource: wasm });
  assert.match(results[0].errors.join(' '), /nucleotide/i);
});

test('reports non-antibody input without throwing', async () => {
  const { results } = await analyze([{ id: 'x', sequence: 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ' }], 'imgt', { wasmSource: wasm });
  assert.ok(results[0].errors.length > 0 || results[0].domains.length === 0);
});

test('handles multi-domain scFv', async () => {
  const scfv = VH + 'GGGGSGGGGSGGGGS' + VL;
  const { results } = await analyze([{ id: 'scfv', sequence: scfv }], 'kabat', { wasmSource: wasm });
  const chains = results[0].domains.map((d) => d.chainType);
  assert.deepEqual(chains, ['H', 'K']);
});

test('progress callback is invoked', async () => {
  const seen = [];
  await analyze([{ id: 'vh', sequence: VH }], 'imgt', { wasmSource: wasm, onProgress: (d, t) => seen.push([d, t]) });
  assert.ok(seen.length >= 2);
  assert.equal(seen[seen.length - 1][0], seen[seen.length - 1][1]);
});

test('supports cancellation via AbortSignal', async () => {
  const controller = new AbortController();
  controller.abort(new Error('cancelled'));
  await assert.rejects(
    analyze([{ id: 'vh', sequence: VH }], 'imgt', { wasmSource: wasm, signal: controller.signal }),
    /cancelled/,
  );
});

test('getEngineVersion returns a string', () => {
  assert.equal(typeof getEngineVersion(), 'string');
});
