# Attributions

This project is distributed under the **GNU General Public License v3.0** (see `LICENSE`),
because it incorporates rules derived from ANARCI (GPLv3).

## @revanta/anarcism (vendored under `lib/anarcism/`)

- License: Apache-2.0
- Source: https://www.npmjs.com/package/@revanta/anarcism
- Used for chain-type classification and IMGT numbering (WebAssembly).
- See `lib/anarcism/LICENSE` and `lib/anarcism/THIRD_PARTY_NOTICES.md`.

## ANARCI (`js/anarci-port.js`)

- License: GPLv3
- Source: https://github.com/oxpig/ANARCI — `lib/python/anarci/schemes.py`
- `number_kabat_heavy`, `_number_regions`, `smooth_insertions`, `get_cdr3_annotations`,
  and `gap_missing` were ported to JavaScript to apply the official Kabat heavy-chain
  numbering rules to anarcism's IMGT alignment. Cite:
  Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor classification.
  Bioinformatics 32(2):298–300 (2016).

## unplage/humanization (`js/kabat-light.js`, `js/vhh.js`)

- License: MIT
- Source: https://github.com/unplage/humanization — `scripts/humanize/numbering.py`
- The anchor-based Kabat light-chain numbering engine (`number_light` and helpers) and
  the VHH hallmark detector (`is_vhh`) were ported to JavaScript. Used instead of
  ANARCI/AbRSA for light chains because those tools mis-number VL CDR3.

## Reference sources

`docs/reference/` contains the upstream files used to author the ports, kept for
traceability:
- `anarci_schemes.py` — oxpig/ANARCI, GPLv3
- `humanization_numbering.py` — unplage/humanization, MIT

## Golden test data

`tests/fixtures/anarci_kabat_reference.json` was generated with the official ANARCI
package. `tests/fixtures/vl_reference.json` and `vl_ladder.json` were generated with the
original unplage/humanization `number_light`.
