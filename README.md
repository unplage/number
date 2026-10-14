# Antibody Numbering

A fully offline, backend-free web app (PWA) for numbering antibody variable-region
sequences and annotating their CDRs with the **IMGT** or **Kabat** schemes.

- Paste FASTA or upload a `.fa` / `.fasta` / `.txt` file.
- Chain typing (H / K / L, plus VHH hallmark detection) and numbering run **locally** in
  the browser via WebAssembly — no server, no API, no data leaves the device.
- After the first visit the whole app (including the ~1 MB numbering engine) is cached
  by a service worker and works **completely offline**.
- Results can be downloaded as TSV, JSON, or aligned numbered text.

## How it works

| Concern | Engine |
| --- | --- |
| Chain typing + IMGT numbering | [`@revanta/anarcism`](https://www.npmjs.com/package/@revanta/anarcism) (Rust→WASM, IMGT-identical to ANARCI) |
| Kabat — heavy chain (VH) | Rules ported from [`oxpig/ANARCI`](https://github.com/oxpig/ANARCI) `number_kabat_heavy` (GPLv3) |
| Kabat — light chain (VL) | Anchor-based engine ported from [`unplage/humanization`](https://github.com/unplage/humanization) (MIT). Used because ANARCI/AbRSA are known to mis-number VL CDR3. |
| VHH detection | `is_vhh` hallmarks ported from `unplage/humanization` (MIT) |

The Kabat VH conversion is validated to be **byte-identical to real ANARCI** on a set of
diverse VH sequences (standard, long CDR3 with insertions, VHH). See `tests/`.

## Usage

Open the deployed app, paste sequences, pick a scheme, and analyze. Nothing to install.

## Development

There is **no build step** — the app is plain HTML/CSS/ES modules.

```sh
# run the test suite (Node 20+, no dependencies)
npm test          # or: node --test tests/

# serve locally (service worker works on http://localhost)
npm run serve     # python3 -m http.server 5173
```

Then open <http://localhost:5173>.

## Layout

```
index.html            app shell
manifest.webmanifest  PWA manifest
sw.js                 service worker (precaches all assets for offline use)
css/style.css
js/app.js             UI controller
js/engine.js          orchestrator (analyze)
js/anarci-port.js     ANARCI Kabat heavy-chain conversion (GPLv3)
js/kabat-light.js     Kabat light-chain anchor engine (MIT)
js/vhh.js             VHH hallmark detection (MIT)
js/fasta.js           FASTA parser
js/export.js          TSV / JSON / aligned-text export
lib/anarcism/         vendored anarcism WASM engine (Apache-2.0)
tests/                node:test suite + golden fixtures
docs/reference/       upstream source used for the ports (reference only)
```

## Deployment

Pushing to `main` triggers `.github/workflows/deploy.yml`, which runs the tests and
publishes the static site to GitHub Pages.

## License

GNU GPLv3. The Kabat VH conversion is derived from ANARCI (GPLv3), so the combined work
is GPLv3. See `ATTRIBUTIONS.md`. Not for clinical use.
