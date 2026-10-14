import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { loadEngine, analyze } from '../js/engine.js';

const wasm = readFileSync(new URL('../lib/anarcism/anarcism.wasm', import.meta.url));
const engine = await loadEngine(wasm);

const VH = 'EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS';
const VL = 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK';

test('detects VH chain and produces IMGT numbering', async () => {
  const [r] = await analyze([{ id: 'vh', sequence: VH }], 'imgt', wasm);
  assert.equal(r.errors.length, 0);
  assert.equal(r.domains[0].chainType, 'H');
  assert.ok(r.domains[0].residues.length > 100);
  assert.ok(r.domains[0].cdrs.CDR3.length > 0);
});

test('detects kappa light chain', async () => {
  const [r] = await analyze([{ id: 'vl', sequence: VL }], 'imgt', wasm);
  assert.equal(r.domains[0].chainType, 'K');
  assert.ok(r.domains[0].cdrs.CDR1.length > 0);
  assert.ok(r.domains[0].cdrs.CDR3.length > 0);
});

test('kabat mode numbers VH and VL', async () => {
  const results = await analyze(
    [
      { id: 'vh', sequence: VH },
      { id: 'vl', sequence: VL },
    ],
    'kabat',
    wasm,
  );
  assert.equal(results[0].domains[0].chainType, 'H');
  assert.equal(results[1].domains[0].chainType, 'K');
  assert.ok(results[1].domains[0].residues.some((r) => r.pos === 'L23'));
  assert.equal(results[1].domains[0].cdrs.CDR3, 'QQHYTTPPT');
});

test('reports error for non-antibody input', async () => {
  const [r] = await analyze([{ id: 'x', sequence: 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ' }], 'imgt', wasm);
  assert.ok(r.errors.length > 0 || r.domains.length === 0);
});
