# Antibody Numbering

A fully offline, backend-free web app (PWA) for numbering antibody variable-region
sequences and annotating their CDRs.

- **Six numbering schemes**: IMGT, Kabat, Chothia, Martin (extended Chothia), AHo, Wolfguy.
- Chain typing (H / K / L, plus VHH hallmark detection), germline V/J assignment, numbering
  confidence, CDR loop lengths and sequence-liability (PTM) flags.
- Paste FASTA, upload a file, or drag & drop. Multiple sequences and multi-domain
  constructs (e.g. scFv) are supported.
- Everything runs **locally in the browser** — no server, no API, no data leaves the device.
  After the first visit the whole app (including the ~1 MB engine) is cached and works
  **completely offline**.
- Results start with a **CDR overview table** (one row per sequence/domain: `id  CDR1  CDR2  CDR3`)
  followed by the per-sequence detail, and export to TSV, CSV, JSON, aligned text,
  ANARCI-style output, FASTA + CDRs, and the CDR summary.

## How it works

| Concern | Engine |
| --- | --- |
| Chain typing, germline assignment, IMGT numbering | [`@revanta/anarcism`](https://www.npmjs.com/package/@revanta/anarcism) (Rust→WASM, IMGT-identical to ANARCI) |
| Kabat, Chothia, Martin, AHo, Wolfguy (VH) | Rules ported from [`oxpig/ANARCI`](https://github.com/oxpig/ANARCI) `schemes.py` (BSD-3-Clause) |
| Kabat light chain (VL) | Anchor-based engine ported from [`unplage/humanization`](https://github.com/unplage/humanization) (MIT). Used because ANARCI/AbRSA mis-number VL CDR3. |
| VHH detection | `is_vhh` hallmarks ported from `unplage/humanization` (MIT) |

The Kabat/Chothia/Martin/AHo/Wolfguy conversions are validated to be **identical to real
ANARCI** across a set of diverse VH and VL sequences. See `tests/`.

## Usage

Open the deployed app at <https://unplage.github.io/number/>, paste sequences, pick a
scheme, and analyze. Nothing to install.

## Development

There is **no build step** — the app is plain HTML/CSS/ES modules.

```sh
npm test          # node --test tests/   (Node 20+, zero runtime deps)
npm run typecheck # tsc -p tsconfig.json (optional; needs `npm install`)
npm run serve     # python3 -m http.server 5173  -> http://localhost:5173
```

## Layout

```
index.html            app shell
manifest.webmanifest  PWA manifest
sw.js                 service worker (precaches all assets for offline use)
css/style.css
js/app.js             UI controller
js/engine.js          orchestrator (scheme dispatch, germline, confidence, progress/cancel)
js/schemes.js         ANARCI scheme ports (Kabat/Chothia/Martin/AHo/Wolfguy) + state vector
js/kabat-light.js     Kabat light-chain anchor engine (MIT)
js/vhh.js             VHH hallmark detection (MIT)
js/validate.js        input validation (DNA, non-standard residues, tags)
js/analysis.js        CDR lengths and liability (PTM) motifs
js/fasta.js           FASTA parser
js/export.js          TSV/CSV/JSON/aligned/ANARCI/FASTA export
js/number-worker.js   Web Worker running numbering off the main thread
js/data/blosum62.js   Blosum62 subset (ANARCI, BSD-3-Clause)
lib/anarcism/         vendored anarcism WASM engine (Apache-2.0)
tests/                node:test suite + golden fixtures
tools/                scripts to regenerate the golden fixtures
docs/reference/       upstream source used for the ports (reference only)
```

## Deployment

GitHub Pages serves this repository directly from the `main` branch root (no build step).
The `unplage/number` token originally used had no `workflow` scope, so there is no GitHub
Actions workflow. Pushing to `main` updates the live site.

## License

GNU GPLv3 (all bundled components are permissive and GPL-compatible). See
`ATTRIBUTIONS.md`. Not for clinical use.
