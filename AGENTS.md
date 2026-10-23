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
  (Cys23→Trp35→Cys88→Phe98) that runs on the raw domain sequence. Kabat light chains
  are deliberately not routed through `schemes.js`; the anchor engine keeps the VL
  handling in one place. The engine is a port of `unplage/humanization` `number_light`
  (MIT) with these deliberate deviations, all covered by
  `tests/kabat-light.test.mjs`:
  1. The J-anchor (Phe98) search window is 45 residues instead of 16. With the upstream
     window the `FG.G` motif is truncated for CDR3s of 14 residues or more, the anchor is
     lost and the whole J region is swallowed by CDR3 (FR4 silently empty).
  2. `findJAnchor` evaluates all three motifs before choosing (as the Python reference
     does) instead of stopping after the first hit, and its J-consensus tail is fixed to
     match the canonical light-chain J tails (`QGTKVE`, `GGTKLT`).
  3. When the IMGT alignment is supplied (`numberLight(seq, {imgt})`, wired up in
     `js/engine.js`), the FR1 is numbered from the alignment: the conserved Cys is
     placed on L23 even for light chains with an FR1 deletion (lambda, IGLV2/3/…) and
     FR1 gaps (e.g. the missing L10) are preserved instead of being closed up.
  4. The FR1 is capped at 23 residues only when the Cys anchor comes from the fallback
     scan (no alignment). An aligned Cys (including an FR1 insertion before it) stays in
     FR1 at L23 instead of being pushed into CDR1.
  5. Short CDR1/CDR3 loops use ANARCI's ordered deletion tables. For CDR3 this keeps
     L96/L97 for loops of 8 residues or less; upstream `number_light` numbered them
     forwards and dropped those positions.
  6. With an IMGT alignment, Phe98 is the first FR4 residue (IMGT >= 118), so a
     CDR3-internal `FG`/`FXXG` motif cannot become the J anchor, and FR2/CDR2/FR3 labels
     are derived from the IMGT states: CDR2 insertions get L52A…, FR3 insertions keep
     their actual position, and deletions leave gaps instead of shifting CDR2. The
     aligned Cys23 is only trusted when the residue is actually Cys.
  Without the alignment the engine falls back to the anchor heuristics (CDR1/CDR3 still
  use the ANARCI deletion order). `tests/fixtures/vl_reference.json` (generated from the
  upstream Python) still matches residue-for-residue because it contains no short loops.
  The other intentional divergence is `vl_ladder.json` `cdr3_5` (14-residue CDR3), where
  upstream loses the J region — see the test comments.

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
  tests; the app intentionally uses `kabat-light.js`. Short CDR3 loops are intentionally
  numbered with ANARCI's deletion order (L96/L97 preserved), which differs from upstream
  `number_light`; keep those tests green.
- **Do not change the anarcism output → state-vector reconstruction casually.** Parity with
  real ANARCI for every scheme depends on it. Keep `tests/schemes.test.mjs` green.
- **Never derive the deletion states from `domain.paddedImgtAlignment`.** That string is the
  *sequence* with gaps: HMM insertion residues are inlined in it, so its length is
  128 + (number of insertions) (145 for a long-CDR3 VH, 127 for a light chain) and its
  characters are shifted by the number of preceding insertions. Reading "is this state
  deleted?" from it mis-classifies every state after an insertion that is followed by a
  deletion — one residue is duplicated and another dropped, so the numbered output gets
  the wrong amino acid. The state vector must be rebuilt from `domain.numbering`
  (`matchByPos`): a state is a deletion iff no residue carries that IMGT position.
  `tests/state-vector.test.mjs` guards this.
- anarcism's IMGT CDR3 insertions (labels `111x`/`112x`) are attached **inside the state
  loop at `s === 111`** in `domainToStateVector`. Attaching them at the end of the vector
  (a past bug) pushes them into FR4 and breaks long-CDR3 numbering.
- Python-vs-JS ordering matters in scheme ports: several ANARCI blocks **slice before
  sorting** (e.g. AHo FW1). Preserve the exact order of operations.
- **The CDR boundaries in `js/engine.js` (`CDR_RANGES`) are per scheme, not per IMGT.**
  Kabat/Chothia/Martin use the published narrow definitions (Kabat H CDR1 31-35, …);
  AHo uses the AHo scheme's own CDR *regions* (25-42 / 58-77 / 107-138) because the
  AHo numbering (1-149) is *defined* by those regions — the narrow IMGT-style ranges are
  not expressible in AHo coordinates (AHo 27-38 would give "GFNIKD" for trastuzumab).
  So an AHo CDR length is wider than the same loop in Kabat and is not comparable with
  it. The Wolfguy CDRH3 is `{331, 332}` plus 351-399, so its range starts at 331.

## Testing

`npm test` (121 tests). Golden fixtures in `tests/fixtures/` were generated from the real
upstream tools and can be regenerated with the committed scripts:

- `anarci_all_schemes.json` — real ANARCI for Kabat/Chothia/Martin/AHo/Wolfguy, VH and VL
  (`python3 tools/gen_anarci_fixtures.py`, needs `pip install anarci` + HMMER).
- `vl_reference.json` / `vl_ladder.json` — original unplage/humanization `number_light`
  (`python3 tools/gen_vl_fixtures.py`).

`tests/schemes.test.mjs` asserts the JS ports match ANARCI residue-for-residue. When
touching `js/schemes.js`, regenerate and re-diff. Reference upstream sources are in
`docs/reference/`. `tests/kabat-light.test.mjs` additionally covers short CDR1/CDR3,
CDR2 insertions, FR2/FR3 indels, and the IMGT FR4/Cys23 anchors. `tests/engine.test.mjs`
exercises the same cases end-to-end through anarcism (including lambda and long CDR3).

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
