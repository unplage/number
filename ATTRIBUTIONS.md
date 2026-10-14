# Attributions

This project is distributed under the **GNU General Public License v3.0** (see `LICENSE`).
All bundled components are permissively licensed (BSD-3-Clause, Apache-2.0, MIT) and are
compatible with the GPLv3.

## @revanta/anarcism (vendored under `lib/anarcism/`)

- License: Apache-2.0
- Source: https://www.npmjs.com/package/@revanta/anarcism
- Used for chain-type classification, germline assignment and IMGT numbering (WebAssembly).
- See `lib/anarcism/LICENSE` and `lib/anarcism/THIRD_PARTY_NOTICES.md`.

## ANARCI (`js/schemes.js`, `js/data/blosum62.js`)

- License: **BSD 3-Clause** (`Copyright (C) 2016 Oxford Protein Informatics Group (OPIG)`)
- Source: https://github.com/oxpig/ANARCI — `lib/python/anarci/schemes.py`
- The Kabat, Chothia, Martin (extended Chothia), AHo and Wolfguy numbering schemes
  (`_number_regions`, `smooth_insertions`, `get_cdr3_annotations`, `gap_missing` and the
  per-scheme renumbering functions) were ported to JavaScript to apply ANARCI's official
  rules to anarcism's IMGT alignment. The Blosum62 subset used by the Wolfguy L1
  assignment was extracted from the same file. Cite:
  Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor classification.
  Bioinformatics 32(2):298–300 (2016).

## unplage/humanization (`js/kabat-light.js`, `js/vhh.js`)

- License: MIT
- Source: https://github.com/unplage/humanization — `scripts/humanize/numbering.py`
- The anchor-based Kabat light-chain numbering engine (`number_light` and helpers) and the
  VHH hallmark detector (`is_vhh`) were ported to JavaScript. Used instead of ANARCI/AbRSA
  for light chains because those tools mis-number VL CDR3.

## Reference sources

`docs/reference/` contains the upstream files used to author the ports, kept for
traceability:
- `anarci_schemes.py` — oxpig/ANARCI, BSD-3-Clause
- `humanization_numbering.py` — unplage/humanization, MIT

## Golden test data

Generated from the real upstream tools, and reproducible via the scripts in `tools/`:
- `tests/fixtures/anarci_all_schemes.json` — real ANARCI, all schemes (`tools/gen_anarci_fixtures.py`).
- `tests/fixtures/vl_reference.json`, `vl_ladder.json` — original `number_light` (`tools/gen_vl_fixtures.py`).
