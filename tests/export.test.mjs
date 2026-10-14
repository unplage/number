import { test } from 'node:test';
import assert from 'node:assert/strict';
import { toTsv, toCsv, toJson, toAligned, toAnarci, toFastaWithCdr, toCdrSummary } from '../js/export.js';

const results = [
  {
    id: 'vh',
    scheme: 'kabat',
    errors: [],
    warnings: [],
    domains: [
      {
        chainType: 'H',
        chainLabel: 'Heavy (IGH)',
        species: 'human',
        start: 0,
        end: 3,
        bitScore: 210.5,
        confidence: 'high',
        germline: { vGene: 'IGHV1-2*01', vIdentity: 0.9, jGene: 'IGHJ4*01' },
        residues: [
          { pos: 'H1', aa: 'E', region: 'FR1' },
          { pos: 'H2', aa: 'V', region: 'FR1' },
          { pos: 'H31', aa: 'D', region: 'CDR1' },
        ],
        cdrs: { CDR1: 'D', CDR2: '', CDR3: '' },
        cdrLengths: { CDR1: 1, CDR2: 0, CDR3: 0 },
        liabilities: [{ type: 'Oxidation risk', motif: 'M', region: 'FR1', positions: 'H1', note: 'x', severity: 'low' }],
      },
    ],
  },
];

test('TSV has header and one row per residue', () => {
  const tsv = toTsv(results);
  const lines = tsv.trim().split('\n');
  assert.match(lines[0], /seq_id\tchain\tspecies/);
  assert.equal(lines.length, 4); // header + 3 residues
  assert.match(lines[1], /^vh\tH\thuman\tIGHV1-2\*01\t0\.900/);
});

test('CSV quotes when needed and parses shape', () => {
  const csv = toCsv(results);
  const lines = csv.trim().split('\n');
  assert.equal(lines.length, 4);
  assert.match(lines[0], /^seq_id,chain,species/);
});

test('JSON round-trips and includes metadata', () => {
  const obj = JSON.parse(toJson(results, { engineVersion: '1.0.1' }));
  assert.equal(obj.results.length, 1);
  assert.equal(obj.engineVersion, '1.0.1');
  assert.equal(obj.results[0].domains[0].residues.length, 3);
});

test('aligned text contains positions and CDR summary', () => {
  const text = toAligned(results);
  assert.match(text, /H1/);
  assert.match(text, /CDR1=D/);
});

test('ANARCI-style output uses // separators', () => {
  const text = toAnarci(results);
  assert.match(text, /# vh \| chain=H/);
  assert.match(text, /H1\tE\tFR1/);
  assert.match(text, /\/\//);
});

test('annotated FASTA carries CDR header', () => {
  const text = toFastaWithCdr(results);
  assert.match(text, /^>vh_H \| scheme=kabat/);
  assert.match(text, /CDR1=D/);
  assert.match(text, /\nEVD/);
});

test('CDR summary has id + CDR1/2/3 columns', () => {
  const lines = toCdrSummary(results).split('\n');
  assert.equal(lines[0], 'id\tCDR1\tCDR2\tCDR3');
  assert.equal(lines[1], 'vh\tD\t\t');
});

test('CDR summary labels multi-domain records with the chain', () => {
  const multi = [
    {
      id: 'scfv',
      scheme: 'kabat',
      errors: [],
      domains: [
        { chainType: 'H', cdrs: { CDR1: 'A', CDR2: 'B', CDR3: 'C' } },
        { chainType: 'K', cdrs: { CDR1: 'D', CDR2: 'E', CDR3: 'F' } },
      ],
    },
  ];
  const lines = toCdrSummary(multi).trim().split('\n');
  assert.equal(lines[1], 'scfv (H)\tA\tB\tC');
  assert.equal(lines[2], 'scfv (K)\tD\tE\tF');
});
