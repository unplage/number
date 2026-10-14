import { test } from 'node:test';
import assert from 'node:assert/strict';
import { validateSequence, ERROR_CODES } from '../js/validate.js';

test('accepts a normal protein sequence', () => {
  const v = validateSequence('EVQLVESGGGLVQ');
  assert.equal(v.ok, true);
  assert.equal(v.sequence, 'EVQLVESGGGLVQ');
});

test('normalizes whitespace and case', () => {
  const v = validateSequence('evql vesg\n gglvq');
  assert.equal(v.sequence, 'EVQLVESGGGLVQ');
});

test('flags nucleotide input', () => {
  const v = validateSequence('ATGGCTGATATTGCCGTGTTAAAGCCGTAA');
  assert.equal(v.ok, false);
  assert.equal(v.errors[0].code, ERROR_CODES.DNA);
});

test('flags invalid characters', () => {
  const v = validateSequence('EVQLVESGG$%');
  assert.equal(v.ok, false);
  assert.equal(v.errors[0].code, ERROR_CODES.INVALID_CHARS);
});

test('warns on non-standard residues', () => {
  const v = validateSequence('EVQLVESGGXZB');
  assert.equal(v.ok, true);
  assert.match(v.warnings.join(' '), /Non-standard/);
});

test('warns on poly-His tag', () => {
  const v = validateSequence('HHHHHHEVQLVESGGGLVQPGGSLRLSCAA');
  assert.match(v.warnings.join(' '), /His/);
});

test('empty input errors', () => {
  const v = validateSequence('   ');
  assert.equal(v.ok, false);
  assert.equal(v.errors[0].code, ERROR_CODES.EMPTY);
});
