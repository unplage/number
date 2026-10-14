import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Anarcism } from '../lib/anarcism/index.js';
import { schemeNumberingFromDomain } from '../js/schemes.js';

const wasm = readFileSync(new URL('../lib/anarcism/anarcism.wasm', import.meta.url));
const engine = await Anarcism.create({ source: wasm });
const fixture = JSON.parse(readFileSync(new URL('./fixtures/anarci_all_schemes.json', import.meta.url)));

// Sequences used to generate the ANARCI reference (see fixtures note).
const SEQS = {
  trast_vh: 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS',
  tm_vh: 'EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS',
  long_vh: 'QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARRRRRRRRRRRSSSSSSSSSYYYYYYYYYWGQGTLVTVSS',
  vhh: 'QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAAGDYWGQGTLVTVSS',
  vh2: 'QVQLQQSGAELARPGASVKLSCKASGYTFTSYWMQWVKQRPGQGLEWIGAIYPGDGDTRYTQKFKGKATLTADKSSSTAYMQLSSLASEDSAVYYCARGGTGYFDYWGQGTTLTVSS',
  vh_short: 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYGISWVRQAPGQGLEWMGWISAYNGNTNYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARDLGYSYYYYGMDVWGQGTTVTVSS',
  trast_vl: 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK',
  vk2: 'DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK',
  vl_lambda: 'QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL',
  vl_long_cdr3: 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQWYSSPPPTFGQGTKVEIK',
};

const SCHEMES = ['kabat', 'chothia', 'martin', 'aho', 'wolfguy'];

function normalizeRefNumbering(numbering) {
  const map = new Map();
  for (const [p, aa] of numbering) {
    if (aa === '-') continue;
    map.set(('H' + p).replace(/\s+/g, ''), aa);
  }
  return map;
}

for (const [name, seq] of Object.entries(SEQS)) {
  const dom = engine.numberSequence(seq, { allowedChains: ['H', 'K', 'L'] }).domains[0];
  for (const scheme of SCHEMES) {
    const key = `${name}::${scheme}`;
    const ref = fixture[key];
    if (!ref || ref.error) continue;
    test(`${name} ${scheme} matches ANARCI`, () => {
      const mine = new Map(
        schemeNumberingFromDomain(scheme, dom.chainType, dom).map((r) => [
          'H' + r.n + (r.ins && r.ins !== ' ' ? r.ins : ''),
          r.aa,
        ]),
      );
      const refMap = normalizeRefNumbering(ref.numbering);
      assert.equal(mine.size, refMap.size, 'residue count');
      for (const [pos, aa] of refMap) assert.equal(mine.get(pos), aa, `position ${pos}`);
    });
  }
}
