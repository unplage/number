import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { numberLight } from '../js/kabat-light.js';

const fixture = JSON.parse(readFileSync(new URL('./fixtures/vl_reference.json', import.meta.url)));
const ladder = JSON.parse(readFileSync(new URL('./fixtures/vl_ladder.json', import.meta.url)));

const VL = {
  vl1: 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK',
  vl2: 'DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK',
  vl_insert: 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQWYSSPPTFGQGTKVEIK',
  lambda: 'QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL',
};

// Reference outputs generated with the original unplage/humanization number_light (MIT).
for (const [name, seq] of Object.entries(VL)) {
  test(`kabat light matches Python reference: ${name}`, () => {
    const { residues } = numberLight(seq);
    const mine = residues.map((r) => [r.pos, r.aa, r.region]);
    assert.deepEqual(mine, fixture[name]);
  });
}

test('VL CDR3 insertion ladder matches Python reference', () => {
  // The upstream reference looks for the J anchor in a 16-residue window, which
  // truncates the "FG.G" motif once the CDR3 reaches 14 residues: from
  // cdr3_5 on the reference loses the J region and turns the whole tail into
  // CDR3. The fixed engine keeps FR4 there - see the dedicated test below.
  const base = VL.vl1;
  for (const [name, ref] of Object.entries(ladder)) {
    if (name === 'cdr3_5') continue;
    const extra = Number(name.split('_')[1]);
    const cdr3 = 'QQHYTTP' + 'A'.repeat(extra) + 'PT';
    const seq = base.replace('QQHYTTPPT', cdr3);
    const { residues } = numberLight(seq);
    assert.deepEqual(residues.map((r) => [r.pos, r.aa]), ref, name);
  }
});

test('the 14-residue CDR3 step keeps the J region (upstream reference loses it)', () => {
  const extra = 5;
  const cdr3 = 'QQHYTTP' + 'A'.repeat(extra) + 'PT';
  const seq = VL.vl1.replace('QQHYTTPPT', cdr3);
  const { residues } = numberLight(seq);

  const cdr3Residues = residues.filter((r) => r.region === 'CDR3');
  const fr4 = residues.filter((r) => r.region === 'FR4');
  assert.equal(cdr3Residues.length, 14, 'CDR3 must be 14 residues');
  assert.deepEqual(
    cdr3Residues.map((r) => r.pos),
    ['L89', 'L90', 'L91', 'L92', 'L93', 'L94', 'L95', 'L95A', 'L95B', 'L95C', 'L95D', 'L95E', 'L96', 'L97'],
  );
  assert.deepEqual(
    fr4.map((r) => r.pos),
    ['L98', 'L99', 'L100', 'L101', 'L102', 'L103', 'L104', 'L105', 'L106', 'L107'],
  );

  // The upstream golden output for this step swallowed the J region: the
  // CDR3 ran from L95 up to L95O instead of ending at L97.
  const upstream = ladder.cdr3_5;
  assert.equal(upstream.filter(([p]) => p.startsWith('L95')).length, 16);
  assert.notDeepEqual(residues.map((r) => [r.pos, r.aa]), upstream);
});

test('conserved anchors are numbered at Kabat positions', () => {
  const { residues } = numberLight(VL.vl1);
  const byPos = new Map(residues.map((r) => [r.pos, r.aa]));
  assert.equal(byPos.get('L23'), 'C');
  assert.equal(byPos.get('L35'), 'W');
  assert.equal(byPos.get('L88'), 'C');
});

// ---------------------------------------------------------------------------
// Regressions for the anchor engine
// ---------------------------------------------------------------------------

// IMGT alignment of `seqLen` residues starting at state `from`; the states in
// `gaps` are deleted (i.e. the framework gap of a lambda light chain).
function imgtAlignment(seqLen, { from = 1, gaps = [] } = {}) {
  const gapSet = new Set(gaps);
  const map = [];
  for (let state = from; map.length < seqLen; state += 1) {
    if (!gapSet.has(state)) map.push({ pos: state, ins: '' });
  }
  return map;
}

const VL_PREFIX = 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYC';

test('FR1 follows the alignment: the conserved Cys lands on L23 for lambda chains', () => {
  // A lambda light chain is one FR1 residue shorter than the Kabat reference:
  // IMGT position 10 is a gap and the conserved Cys is the 22nd residue.
  const imgt = imgtAlignment(VL.lambda.length, { gaps: [10] });
  const { residues, warnings } = numberLight(VL.lambda, { imgt });
  const byPos = new Map(residues.map((r) => [r.pos, r.aa]));

  assert.equal(byPos.get('L23'), 'C', 'conserved Cys must be numbered L23');
  assert.equal(byPos.has('L10'), false, 'the FR1 gap at L10 must be preserved');
  assert.equal(byPos.get('L11'), 'V');
  assert.equal(byPos.get('L22'), 'S');
  assert.equal(byPos.get('L24'), 'S', 'CDR1 starts at L24');

  const labels = residues.map((r) => r.pos);
  assert.equal(new Set(labels).size, labels.length, 'no duplicated position labels');
  assert.ok(
    !warnings.some((w) => /Cys23 anchor/.test(w)),
    'no anchor warning when the alignment is used',
  );
});

test('a shifted Cys anchor cannot assign the same position twice', () => {
  // One extra N-terminal residue: the Cys sits at position 24, so the FR1 must
  // be capped at 23 residues instead of running into the CDR1 start.
  const seq = 'X' + VL.vl1;
  const { residues, warnings } = numberLight(seq);
  const labels = residues.map((r) => r.pos);
  assert.equal(new Set(labels).size, labels.length, 'no duplicated position labels');
  assert.equal(labels.filter((p) => p === 'L24').length, 1);
  assert.ok(
    warnings.some((w) => /FR1 capped to 23 residues/.test(w)),
    'the capped FR1 must be reported',
  );
});

test('a long CDR3 does not swallow the J region', () => {
  // 14-residue CDR3: the reference search window (16 residues) truncated the
  // "FGQG" motif, so the anchor was lost and the whole J region became CDR3.
  const cdr3 = 'QYTYTYTYTYTYTY';
  const seq = VL_PREFIX + cdr3 + 'FGQGTKVEIK';
  const { residues } = numberLight(seq);

  const cdr3Residues = residues.filter((r) => r.region === 'CDR3');
  const fr4 = residues.filter((r) => r.region === 'FR4');
  assert.equal(cdr3Residues.length, cdr3.length, 'CDR3 must not swallow the J region');
  assert.equal(fr4.length, 10, 'FR4 must be complete');
  assert.equal(fr4[0].pos, 'L98');
  assert.equal(fr4[0].aa, 'F');
  assert.equal(residues.at(-1).pos, 'L107');
});

test('a broken J region is reported instead of silently swallowing the J region', () => {
  // The conserved Phe of the J region is mutated, so no anchor can be found -
  // but the alignment still places residues in the FR4 region.
  const seq = VL_PREFIX + 'QQHYTTPPT' + 'YGQGTKVEIK';
  const imgt = imgtAlignment(seq.length - 10).concat(
    Array.from({ length: 10 }, (_, i) => ({ pos: 118 + i, ins: '' })),
  );
  assert.throws(() => numberLight(seq, { imgt }), /J-region/);
});

test('the J anchor is taken from the J region, not from an FG motif inside CDR3', () => {
  // The CDR3 contains an "FG" pair; the real J region follows it.
  const seq = VL_PREFIX + 'QQFGHTTPPT' + 'FGQGTKVEIK';
  const { residues } = numberLight(seq);
  const fr4 = residues.filter((r) => r.region === 'FR4');
  assert.equal(fr4.length, 10, 'FR4 must start at the real J region');
  assert.equal(fr4[0].pos, 'L98');
  assert.equal(fr4[0].aa, 'F');
  assert.equal(residues.filter((r) => r.region === 'CDR3').length, 10);
});
