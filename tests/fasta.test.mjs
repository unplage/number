import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseFasta } from '../js/fasta.js';

test('parses multiple records with wrapped lines', () => {
  const text = '>seq1 description here\nEVQLV\nESGGG\n>seq2\nDIQMT\n';
  const recs = parseFasta(text);
  assert.equal(recs.length, 2);
  assert.equal(recs[0].id, 'seq1');
  assert.equal(recs[0].description, 'description here');
  assert.equal(recs[0].sequence, 'EVQLVESGGG');
  assert.equal(recs[1].id, 'seq2');
  assert.equal(recs[1].sequence, 'DIQMT');
});

test('handles CRLF and blank lines', () => {
  const text = '>a\r\nEVQL\r\n\r\nVESG\r\n';
  const recs = parseFasta(text);
  assert.equal(recs.length, 1);
  assert.equal(recs[0].sequence, 'EVQLVESG');
});

test('treats headerless input as a single record', () => {
  const recs = parseFasta('EVQLVESGGG');
  assert.equal(recs.length, 1);
  assert.equal(recs[0].sequence, 'EVQLVESGGG');
});

test('drops empty records', () => {
  const recs = parseFasta('>empty\n>full\nEVQL\n');
  assert.equal(recs.length, 1);
  assert.equal(recs[0].id, 'full');
});
