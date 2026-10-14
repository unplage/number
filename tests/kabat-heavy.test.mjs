import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { Anarcism } from '../lib/anarcism/index.js';
import { kabatHeavyFromDomain } from '../js/anarci-port.js';

const wasm = readFileSync(new URL('../lib/anarcism/anarcism.wasm', import.meta.url));
const engine = await Anarcism.create({ source: wasm });
const fixture = JSON.parse(readFileSync(new URL('./fixtures/anarci_kabat_reference.json', import.meta.url)));

const VH = {
  trast_vh: 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS',
  tm_vh: 'EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS',
  long_vh: 'QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARRRRRRRRRRRSSSSSSSSSYYYYYYYYYWGQGTLVTVSS',
  vhh: 'QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAAGDYWGQGTLVTVSS',
  vh2: 'QVQLQQSGAELARPGASVKLSCKASGYTFTSYWMQWVKQRPGQGLEWIGAIYPGDGDTRYTQKFKGKATLTADKSSSTAYMQLSSLASEDSAVYYCARGGTGYFDYWGQGTTLTVSS',
};

// Reference outputs generated with real ANARCI (oxpig/ANARCI), scheme=kabat.
for (const [name, seq] of Object.entries(VH)) {
  test(`kabat heavy matches ANARCI: ${name}`, () => {
    const domain = engine.numberSequence(seq, { allowedChains: ['H'] }).domains[0];
    const mine = new Map(kabatHeavyFromDomain(domain).filter((x) => x.aa && x.aa !== '-').map((x) => [x.pos, x.aa]));
    const ref = new Map();
    for (const [p, aa] of fixture[name].numbering) {
      if (aa === '-') continue;
      ref.set(('H' + p).replace(/\s+/g, ''), aa);
    }
    assert.equal(mine.size, ref.size, 'residue count');
    for (const [pos, aa] of ref) {
      assert.equal(mine.get(pos), aa, `position ${pos}`);
    }
  });
}
