# AGENTS.md

## What this is

Offline-first PWA that numbers antibody variable regions and annotates CDRs with six
schemes (IMGT, Kabat, Chothia, Martin, AHo, Wolfguy). Pure static HTML/CSS/ES modules —
**no build step, no bundler, no backend**. Deploys to GitHub Pages from `main`.

## Commands

```sh
npm test          # node --test tests/  (Node 20+, zero runtime deps)
npm run typecheck # tsc -p tsconfig.json  (optional; run `npm install` first)
npm run serve     # python3 -m http.server 5173  -> http://localhost:5173
```

No `dev`/`build`/`lint` step exists. Do not add a framework or bundler without a strong
reason; the whole app is served as-is.

## Architecture (read before changing numbering output)

`js/engine.js` is the orchestrator. All schemes except IMGT are produced by converting
anarcism's IMGT output into the 128-state HMM vector and applying a numbering-scheme port:

- **IMGT** — `@revanta/anarcism` (vendored WASM in `lib/anarcism/`), residue-identical to
  ANARCI. Region labels come straight from anarcism.
- **Kabat, Chothia, Martin, AHo, Wolfguy** — `js/schemes.js`, a port of ANARCI
  `schemes.py`. `domainToStateVector()` rebuilds the state vector from the anarcism domain
  (`paddedImgtAlignment` + per-residue `position`/`insertionCode`); `convertStateVector()`
  dispatches to the scheme.
- **Kabat light chain is special** — `js/kabat-light.js` is a standalone anchor engine
  (Cys23→Trp35→Cys88→Phe98) that ignores anarcism entirely and runs on the raw domain
  sequence. This is deliberate: ANARCI/AbRSA mis-number Kabat VL CDR3.

`js/vhh.js` flags camelid VHH sequences (Kabat FR2 positions 37/44/45/47). `js/validate.js`
rejects nucleotide / malformed input and flags tags. `js/analysis.js` computes CDR lengths
and PTM/liability motifs.

**Numbering dispatch** (`js/engine.js`): batches of `>= POOL_THRESHOLD` (100) sequences are
numbered by anarcism's multicore `AnarcismWorkerPool` (vendored `lib/anarcism/worker-pool.js`
+ `worker.js`); smaller batches and environments without `Worker` use the in-process WASM
engine (`loadEngine`). An `allowedSpecies` option (UI: species filter) restricts chain typing
and roughly halves the time. The pool is created lazily and any pool failure falls back to the
synchronous engine.

## Critical conventions

- **Do not "fix" the Kabat-VL split.** Do not route Kabat light chains through
  `schemes.js`. `convertStateVector('kabat','L',…)` exists only for scheme comparison in
  tests; the app intentionally uses `kabat-light.js`.
- **Do not change the anarcism output → state-vector reconstruction casually.** Parity with
  real ANARCI for every scheme depends on it. Keep `tests/schemes.test.mjs` green.
- anarcism's IMGT CDR3 insertions (labels `111x`/`112x`) are attached **inside the state
  loop at `s === 111`** in `domainToStateVector`. Attaching them at the end of the vector
  (a past bug) pushes them into FR4 and breaks long-CDR3 numbering.
- Python-vs-JS ordering matters in scheme ports: several ANARCI blocks **slice before
  sorting** (e.g. AHo FW1). Preserve the exact order of operations.

## Testing

`npm test` (89 tests). Golden fixtures in `tests/fixtures/` were generated from the real
upstream tools and can be regenerated with the committed scripts:

- `anarci_all_schemes.json` — real ANARCI for Kabat/Chothia/Martin/AHo/Wolfguy, VH and VL
  (`python3 tools/gen_anarci_fixtures.py`, needs `pip install anarci` + HMMER).
- `vl_reference.json` / `vl_ladder.json` — original unplage/humanization `number_light`
  (`python3 tools/gen_vl_fixtures.py`).

`tests/schemes.test.mjs` asserts the JS ports match ANARCI residue-for-residue. When
touching `js/schemes.js`, regenerate and re-diff. Reference upstream sources are in
`docs/reference/`.

## Licensing (important)

The project is **GPLv3**. `js/schemes.js` derives from ANARCI, which is **BSD-3-Clause**
(not GPLv3); the humanization ports are MIT and anarcism is Apache-2.0. All are
GPL-compatible. Keep `ATTRIBUTIONS.md` in sync with any newly ported code.

## Deployment / PWA gotchas

- GitHub Pages serves the repo directly from the `main` branch root (`.nojekyll` present),
  so there is no CI/build step. The `unplage/number` token in this project's history lacked
  the `workflow` scope, which is why there is no GitHub Actions workflow.
- `sw.js` precaches every asset by URL. **When you add, rename, or remove a file, update the
  `ASSETS` list and bump the `CACHE` version string**, or returning users get stale files.
- All paths are relative (`./`) so the app works from the Pages subpath `/number/`. Keep it
  that way — do not switch to absolute `/…` URLs.
- anarcism requires WebAssembly SIMD128 (all evergreen browsers since ~2021–22).
- The service worker only activates over `https://` or `http://localhost`.
