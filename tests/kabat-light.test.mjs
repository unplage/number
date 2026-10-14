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
  const base = VL.vl1;
  for (const [name, ref] of Object.entries(ladder)) {
    const extra = Number(name.split('_')[1]);
    const cdr3 = 'QQHYTTP' + 'A'.repeat(extra) + 'PT';
    const seq = base.replace('QQHYTTPPT', cdr3);
    const { residues } = numberLight(seq);
    assert.deepEqual(residues.map((r) => [r.pos, r.aa]), ref, name);
  }
});

test('conserved anchors are numbered at Kabat positions', () => {
  const { residues } = numberLight(VL.vl1);
  const byPos = new Map(residues.map((r) => [r.pos, r.aa]));
  assert.equal(byPos.get('L23'), 'C');
  assert.equal(byPos.get('L35'), 'W');
  assert.equal(byPos.get('L88'), 'C');
});
