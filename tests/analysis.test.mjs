import { test } from 'node:test';
import assert from 'node:assert/strict';
import { cdrLengths, findLiabilities } from '../js/analysis.js';

const residues = (pairs) => pairs.map(([pos, aa, region]) => ({ pos, aa, region }));

test('cdrLengths counts residues per CDR', () => {
  assert.deepEqual(cdrLengths({ CDR1: 'ABC', CDR2: 'DE', CDR3: '' }), { CDR1: 3, CDR2: 2, CDR3: 0 });
});

test('detects an N-glycosylation sequon', () => {
  const r = residues([
    ['H50', 'N', 'CDR2'],
    ['H51', 'S', 'CDR2'],
    ['H52', 'T', 'CDR2'],
  ]);
  const liabs = findLiabilities(r, 'H');
  assert.ok(liabs.some((l) => l.type === 'N-glycosylation sequon' && l.motif === 'NST'));
});

test('detects deamidation and oxidation', () => {
  const r = residues([
    ['H30', 'N', 'FR1'],
    ['H31', 'G', 'FR1'],
    ['H40', 'M', 'FR2'],
    ['H41', 'W', 'FR2'],
  ]);
  const types = findLiabilities(r, 'H').map((l) => l.type);
  assert.ok(types.includes('Deamidation hotspot'));
  assert.ok(types.includes('Oxidation risk'));
});

test('flags free cysteine in a CDR', () => {
  const r = residues([
    ['H95', 'C', 'CDR3'],
    ['H96', 'A', 'CDR3'],
  ]);
  const liabs = findLiabilities(r, 'H');
  assert.ok(liabs.some((l) => l.type === 'Free cysteine (CDR)'));
});
