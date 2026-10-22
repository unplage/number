import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Anarcism } from '../lib/anarcism/index.js';
import { domainToStateVector, schemeNumberingFromDomain } from '../js/schemes.js';

const wasm = readFileSync(new URL('../lib/anarcism/anarcism.wasm', import.meta.url));
const engine = await Anarcism.create({ source: wasm });

function domainOf(seq) {
  const res = engine.numberSequence(seq, { allowedChains: ['H', 'K', 'L'] });
  assert.ok(res.domains.length, 'expected a domain to be recognised');
  return res.domains[0];
}

const SEQS = {
  // A VH with an FR3 insertion (IMGT 87A) followed by a deleted state
  // (IMGT 111). The padded IMGT alignment is 129 characters long here and
  // shifted by the insertion, so deciding "deleted or not" from that string
  // duplicates one residue and drops another.
  vh_fr3_insertion:
    'EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLNYLQMNSLRAEDTAVYYCARDRGGSYAMDYWGQGTLVTVSS',
  // A kappa chain with an extra residue in FR2 ("WYQQK" -> "WYQEQK"), i.e. an
  // IMGT insertion in the middle of the framework.
  vl_fr2_insertion:
    'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQEQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK',
  // A VH with a long CDR3 (17 IMGT insertions on 111/112).
  vh_long_cdr3:
    'QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARRRRRRRRRRRSSSSSSSSSYYYYYYYYYWGQGTLVTVSS',
  // Plain chains.
  vh: 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS',
  vl: 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK',
  lambda:
    'QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL',
};

test('the reconstructed state vector keeps every residue, in order', () => {
  for (const [name, seq] of Object.entries(SEQS)) {
    const d = domainOf(seq);
    const { stateVector, sequence } = domainToStateVector(d);

    assert.equal(
      sequence,
      d.numbering.map((r) => r.aminoAcid).join(''),
      `${name}: the reconstructed sequence must follow the anarcism residue order`,
    );
    assert.equal(
      stateVector.filter(([, si]) => si !== null).length,
      d.numbering.length,
      `${name}: no residue may be dropped or duplicated`,
    );

    // A state that carries a residue must be matched exactly once, every other
    // state must be a deletion.
    const withResidue = new Set(
      d.numbering.filter((r) => !r.insertionCode).map((r) => r.position),
    );
    for (let stateId = 1; stateId <= 128; stateId += 1) {
      const entries = stateVector.filter(([st]) => st[0] === stateId);
      const matches = entries.filter(([[, type]]) => type === 'm');
      const insertions = entries.filter(([[, type]]) => type === 'i');
      if (withResidue.has(stateId)) {
        assert.equal(matches.length, 1, `${name}: state ${stateId} must be matched once`);
        assert.notEqual(matches[0][1], null, `${name}: state ${stateId} must point at a residue`);
      } else {
        assert.equal(matches.length, 0, `${name}: state ${stateId} must not be matched`);
        assert.equal(insertions.length, 0, `${name}: state ${stateId} must not carry insertions`);
      }
    }
  }
});

test('an insertion before a deletion does not shift the numbering (VH)', () => {
  const d = domainOf(SEQS.vh_fr3_insertion);
  const num = schemeNumberingFromDomain('kabat', d.chainType, d);
  const byPos = new Map(num.map((r) => [r.n, r.aa]));

  // The residue of IMGT 112 belongs to Kabat 99 - the Gly of IMGT 110 must not
  // be pushed into its place.
  assert.equal(byPos.get(99), 'S');
  assert.equal(
    num
      .filter((r) => r.n >= 95 && r.n <= 102)
      .map((r) => r.aa)
      .join(''),
    'DRGGSYAMDY',
  );
  assert.equal(num.length, d.numbering.length, 'no residue may be lost');
});

test('an insertion before a deletion does not shift the numbering (VL)', () => {
  const d = domainOf(SEQS.vl_fr2_insertion);
  const num = schemeNumberingFromDomain('chothia', d.chainType, d);
  const byPos = new Map(num.map((r) => [r.n, r.aa]));

  // These three positions were lost when the deletion states were read from
  // the shifted padded alignment.
  assert.equal(byPos.get(52), 'S');
  assert.equal(byPos.get(60), 'S');
  assert.equal(byPos.get(67), 'S');
  assert.equal(num.length, d.numbering.length, 'no residue may be lost');
});

test('the CDR3 insertions stay inside the CDR3 region', () => {
  const d = domainOf(SEQS.vh_long_cdr3);
  const { stateVector } = domainToStateVector(d);

  // Every IMGT 111x/112x insertion has to be attached inside the CDR3 block
  // (states 105-117), not at the end of the vector where it would fall into
  // FR4.
  const cdr3Insertions = d.numbering.filter((r) => r.insertionCode);
  assert.ok(cdr3Insertions.length > 10, 'the test sequence must have CDR3 insertions');
  const states = stateVector.map(([s]) => s[0]);
  const lastCdr3State = states.lastIndexOf(117);
  for (const r of cdr3Insertions) {
    const idx = stateVector.findIndex(([s]) => s[1] === 'i' && s[0] === 111);
    assert.ok(idx >= 0 && idx < lastCdr3State, 'CDR3 insertions must stay before FR4');
  }

  const num = schemeNumberingFromDomain('kabat', d.chainType, d);
  const cdr3 = num.filter((r) => r.n >= 95 && r.n <= 102);
  assert.equal(cdr3.map((r) => r.aa).join(''), 'RRRRRRRRRRSSSSSSSSSYYYYYYYYY');
});
