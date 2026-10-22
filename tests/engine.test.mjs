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

test('allowedSpecies restricts chain typing without breaking results', async () => {
  const { results } = await analyze([{ id: 'vh', sequence: VH }], 'imgt', {
    wasmSource: wasm,
    allowedSpecies: ['human'],
  });
  const d = results[0].domains[0];
  assert.equal(d.chainType, 'H');
  assert.equal(d.species, 'human');
  assert.ok(d.cdrs.CDR3.length > 0);
});

const VL_LAMBDA = 'QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL';

test('Kabat VL: the lambda FR1 gap is preserved and the conserved Cys is L23', async () => {
  const { results } = await analyze([{ id: 'lambda', sequence: VL_LAMBDA }], 'kabat', { wasmSource: wasm });
  const d = results[0].domains[0];
  assert.equal(d.errors.length, 0);
  const byPos = new Map(d.residues.map((r) => [r.pos, r.aa]));
  // A lambda FR1 is one residue shorter than the Kabat reference: L10 is a gap
  // and the conserved cysteine keeps its canonical L23.
  assert.equal(byPos.get('L23'), 'C');
  assert.equal(byPos.has('L10'), false);
  assert.equal(byPos.get('L11'), 'V');
  assert.equal(byPos.get('L98'), 'F', 'FR4 starts at the J region');
  // The CDR assignment is unchanged.
  assert.equal(d.cdrs.CDR1, 'SGSSSNIGNNYVS');
  assert.equal(d.cdrs.CDR2, 'DNNKRPS');
  assert.equal(d.cdrs.CDR3, 'GTWDSSLSAVV');
  const labels = d.residues.map((r) => r.pos);
  assert.equal(new Set(labels).size, labels.length, 'no duplicated position labels');
});

test('Kabat VL: a long CDR3 keeps the J region', async () => {
  const seq =
    'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYC' +
    'QQWYSSPPPTAAAA' +
    'FGQGTKVEIK';
  const { results } = await analyze([{ id: 'vl_long', sequence: seq }], 'kabat', { wasmSource: wasm });
  const d = results[0].domains[0];
  assert.equal(d.errors.length, 0);
  assert.equal(d.cdrs.CDR3, 'QQWYSSPPPTAAAA');
  const fr4 = d.residues.filter((r) => r.region === 'FR4');
  assert.equal(fr4.length, 10, 'FR4 must not be swallowed by CDR3');
  assert.equal(fr4[0].pos, 'L98');
  assert.equal(fr4[0].aa, 'F');
});

test('Wolfguy: the heavy CDR3 covers the 331/332 positions', async () => {
  // The Wolfguy CDRH3 is {331, 332} plus 351-399: the first two positions are
  // the N-terminal end of the IMGT CDR3 and must not fall into FR3.
  const { results } = await analyze([{ id: 'vh', sequence: VH }], 'wolfguy', { wasmSource: wasm });
  const d = results[0].domains[0];
  assert.equal(d.errors.length, 0);
  const cdr3 = d.residues.filter((r) => r.region === 'CDR3');
  assert.equal(cdr3[0].pos, 'H331');
  assert.equal(cdr3[1].pos, 'H332');
  assert.equal(d.cdrs.CDR3, 'ARMDY');
  // The FR3 must stop at the conserved cysteine (H330).
  const fr3 = d.residues.filter((r) => r.region === 'FR3');
  assert.equal(fr3.at(-1).pos, 'H330');
});

test('the CDR3 content is the same in the IMGT and Wolfguy numbering', async () => {
  const [imgt, wolfguy] = await Promise.all([
    analyze([{ id: 'vh', sequence: VH }], 'imgt', { wasmSource: wasm }),
    analyze([{ id: 'vh', sequence: VH }], 'wolfguy', { wasmSource: wasm }),
  ]);
  assert.equal(imgt.results[0].domains[0].cdrs.CDR3, wolfguy.results[0].domains[0].cdrs.CDR3);
});
