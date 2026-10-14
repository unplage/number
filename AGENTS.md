# AGENTS.md

## What this is

Offline-first PWA that numbers antibody variable regions and annotates CDRs with IMGT or
Kabat. Pure static HTML/CSS/ES modules — **no build step, no bundler, no backend**. Deploys
to GitHub Pages.

## Commands

```sh
npm test        # node --test tests/  (Node 20+, zero deps)
npm run serve   # python3 -m http.server 5173  -> http://localhost:5173
```

There is no `dev`/`build`/`lint` step. Do not add a framework or bundler without a strong
reason; the whole app is served as-is.

## Architecture (read before changing numbering output)

Two independent numbering paths, chosen at runtime in `js/engine.js`:

- **IMGT** — every chain is numbered by `@revanta/anarcism` (vendored WASM in `lib/anarcism/`),
  which is residue-identical to ANARCI's IMGT. Region labels come straight from anarcism.
- **Kabat** — split by chain:
  - VH: `js/anarci-port.js` reconstructs the 128-state IMGT state vector from the anarcism
    domain (`paddedImgtAlignment` + per-residue `position`/`insertionCode`) and applies
    `number_kabat_heavy` ported from ANARCI `schemes.py`.
  - VL: `js/kabat-light.js` is a standalone anchor engine (Cys23→Trp35→Cys88→Phe98) that
    ignores anarcism entirely and runs on the raw domain sequence.

`js/vhh.js` flags camelid VHH sequences (Kabat FR2 positions 37/44/45/47).

## Critical conventions

- **Do not "fix" the split.** ANARCI/AbRSA are known to mis-number Kabat VL CDR3; that is
  exactly why VL uses the anchor engine. Do not route VL through `anarci-port.js`.
- **Do not change the anarcism output → state-vector reconstruction casually.** Kabat VH
  parity with real ANARCI depends on it. Any change must keep `tests/kabat-heavy.test.mjs`
  green against the ANARCI golden fixtures.
- anarcism's IMGT CDR3 insertions (labels `111x`/`112x`) are attached after match state 111
  in `domainToStateVector`; Kabat CDR3 renumbering is purely length-based, so the attach
  point only needs to stay inside IMGT states 105–117.

## Testing

`npm test`. Golden fixtures in `tests/fixtures/` were generated from the real upstream tools:
- `anarci_kabat_reference.json` — real ANARCI (`scheme=kabat`).
- `vl_reference.json` / `vl_ladder.json` — original unplage/humanization `number_light`.

When touching `js/anarci-port.js`, ideally regenerate the ANARCI fixtures (see below) and
re-diff. Reference upstream sources are in `docs/reference/`.

## Licensing (important)

The project is **GPLv3** because `js/anarci-port.js` derives from ANARCI (GPLv3). The
humanization ports are MIT and anarcism is Apache-2.0 (both compatible). Any newly ported
code must keep the combined work GPLv3-compatible; keep `ATTRIBUTIONS.md` in sync.

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

## Regenerating the ANARCI VH fixtures (optional, needs hmmer)

```sh
pip install anarci           # also needs hmmscan/hmmbuild on PATH
# script used: see git history / docs; writes tests/fixtures/anarci_kabat_reference.json
```
