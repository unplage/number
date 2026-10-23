// Kabat light-chain numbering. The anchor-based engine is a port of
// unplage/humanization number_light (MIT License); short CDR1/CDR3 labels follow
// ANARCI's ordered deletion tables and, when the IMGT alignment is available,
// FR2/CDR2/FR3 are labelled from the IMGT states.
// Anchor-based engine: Cys23 -> Trp35 -> Cys88 -> Phe98.
//
// Deliberate deviations from the Python reference (documented in AGENTS.md):
//   1. The J-anchor (Phe98) search window is widened (16 -> 45 residues). With
//      the reference window the "FG.G" motif is truncated for CDR3s longer than
//      ~12 residues, so the anchor is lost and the whole J region is swallowed
//      by CDR3 (FR4 silently empty).
//   2. findJAnchor evaluates all three motifs before choosing, exactly like the
//      reference does. The original JS port stopped after the first hit.
//   3. When the IMGT alignment is supplied (opts.imgt) the FR1 is numbered from
//      the alignment instead of sequentially. The conserved Cys therefore lands
//      on L23 for light chains with an FR1 deletion (lambda, IGLV2/3/...) and
//      FR1 gaps are preserved instead of being closed up.
//   4. The FR1 is capped at 23 residues only when the Cys anchor comes from the
//      fallback scan; an aligned Cys (or an FR1 insertion before it) is never
//      pushed into CDR1.
//   5. Short CDR1/CDR3 labels use ANARCI's ordered deletion tables. The upstream
//      reference numbered short loops forwards, which dropped L96/L97 for short
//      CDR3 loops.
//   6. With an IMGT alignment, Phe98 is the first FR4 residue (IMGT >= 118),
//      FR2/CDR2/FR3 labels are derived from the IMGT states (CDR2 insertions
//      land on L52A.., FR3 indels keep their position), and an IMGT Cys23 entry
//      is only trusted when the residue is actually Cys.

function letters(n) {
  const out = [];
  for (let i = 0; i < n; i++) {
    const a = Math.floor(i / 26);
    const b = i % 26;
    out.push(a === 0 ? String.fromCharCode(65 + b) : String.fromCharCode(64 + a) + String.fromCharCode(65 + b));
  }
  return out;
}

function insLetter(i) {
  return i < 0 ? '' : letters(i + 1).at(-1);
}

// ANARCI's ordered deletions for light-chain CDR1 (Kabat 24-34).
function cdr1Positions(length) {
  const out = [];
  for (let n = 24; n <= 27 && out.length < length; n++) out.push({ n, ins: '' });
  for (let i = 0; i < Math.max(0, length - 11); i++) out.push({ n: 27, ins: insLetter(i) });
  const tail = [28, 29, 30, 31, 32, 33, 34];
  for (let i = Math.abs(Math.min(0, length - 11)); i < tail.length && out.length < length; i++) {
    out.push({ n: tail[i], ins: '' });
  }
  return out.map((p) => `L${p.n}${p.ins}`);
}

// ANARCI's ordered deletions for light-chain CDR3 (Kabat 89-97). Short loops
// delete 95, 94, ... first, so L96/L97 stay at the C-terminal end.
function cdr3Positions(length) {
  if (length > 35) throw new Error('[VL] too many CDR3 residues for Kabat numbering (>35)');
  const insertions = Math.max(0, length - 9);
  const ordered = [95, 94, 93, 92, 91, 96, 97, 90, 89];
  const items = [];
  for (let i = Math.max(0, 9 - length); i < ordered.length; i++) items.push({ n: ordered[i], ins: '' });
  for (let i = 0; i < insertions; i++) items.push({ n: 95, ins: insLetter(i) });
  items.sort((a, b) => (a.n - b.n) || (a.ins === '' ? -1 : b.ins === '' ? 1 : a.ins < b.ins ? -1 : 1));
  return items.map((p) => `L${p.n}${p.ins}`);
}

// Kabat labels for the W35..C88 block, derived from the IMGT states 41..104.
// This mirrors ANARCI's number_kabat_light state machine: 41-56 are FR2
// (state - 6), 57-67 are the CDR2 block (insertions at 52), and 68-104 are FR3
// (state - 13, minus deleted insertion states). When the IMGT alignment is not
// usable this returns null and the caller falls back to the anchor heuristics.
function kabatSpanLabels(imgt, fromIdx, toIdx) {
  const byState = new Map();
  for (let i = fromIdx; i <= toIdx; i++) {
    const m = imgt[i];
    if (!m || !m.pos || m.pos < 41 || m.pos > 104) return null;
    if (!byState.has(m.pos)) byState.set(m.pos, []);
    byState.get(m.pos).push({ idx: i, ins: m.ins || '' });
  }
  if (!byState.has(41) || !byState.has(104)) return null;
  for (const entries of byState.values()) {
    entries.sort((a, b) => {
      if (a.ins === b.ins) return a.idx - b.idx;
      if (a.ins === '') return -1;
      if (b.ins === '') return 1;
      return a.ins < b.ins ? -1 : 1;
    });
  }

  const isInsertionState = (state) =>
    (state >= 58 && state <= 64) || state === 73 || state === 81 || state === 82;
  const labels = new Map();
  let rels = -6;
  let insertion = -1;
  let previousStateId = 1;
  let previousType = 'd';
  for (let state = 41; state <= 104; state++) {
    if (state === 57) {
      rels = -6;
      insertion = -1;
    } else if (state === 68) {
      rels = -13;
      insertion = -1;
    }
    const entries = byState.get(state);
    if (!entries) {
      previousType = 'd';
      if (isInsertionState(state)) rels -= 1;
      continue;
    }
    for (const entry of entries) {
      if (entry.ins === '') {
        const insState = isInsertionState(state);
        if (insState) {
          if (previousType !== 'd') insertion += 1;
          rels -= 1;
        } else {
          insertion = -1;
        }
        labels.set(entry.idx, { n: state + rels, ins: insState ? insLetter(insertion) : '' });
        previousStateId = state;
        previousType = 'm';
      } else {
        insertion += 1;
        labels.set(entry.idx, { n: previousStateId + rels, ins: insLetter(insertion) });
        previousType = 'i';
      }
    }
  }

  const cdr2 = [];
  for (let state = 57; state <= 67; state++) {
    const entries = byState.get(state);
    if (entries) cdr2.push(...entries);
  }
  if (cdr2.length > 4) {
    const annotations = [{ n: 51, ins: '' }, { n: 52, ins: '' }];
    for (let i = 0; i < cdr2.length - 4; i++) annotations.push({ n: 52, ins: insLetter(i) });
    annotations.push({ n: 53, ins: '' }, { n: 54, ins: '' });
    cdr2.forEach((entry, i) => labels.set(entry.idx, annotations[i]));
  }

  for (let i = fromIdx; i <= toIdx; i++) {
    const label = labels.get(i);
    if (!label || label.n < 35 || label.n > 88) return null;
  }
  return labels;
}

function findFr2Trp(s, start, end) {
  const ctxPat = /^WY[QFL]/;
  const prefer = 11;
  let best = null;
  let bestScore = -1e9;
  const window = s.slice(start, end);
  let idx = window.indexOf('W');
  while (idx !== -1) {
    const i = start + idx;
    const cdr1Len = i - start;
    let score = 0;
    if (ctxPat.test(s.slice(i, i + 3))) score += 10;
    if (cdr1Len >= 6 && cdr1Len <= 17) score += 5 - Math.abs(cdr1Len - prefer);
    if (score > bestScore) {
      best = i;
      bestScore = score;
    }
    idx = window.indexOf('W', idx + 1);
  }
  return best;
}

function findC88(s, wIdx, n) {
  const re = /Y[YF]?C/g;
  re.lastIndex = wIdx;
  let best = null;
  let m;
  const limit = Math.min(wIdx + 90, n);
  const segment = s.slice(wIdx, limit);
  while ((m = re.exec(segment)) !== null) {
    const idx = wIdx + m.index + m[0].length - 1;
    const rel = idx - wIdx;
    if (rel >= 52 && rel <= 59) {
      if (best === null) best = idx;
      else if (Math.abs(rel - 54) < Math.abs(best - wIdx - 54)) best = idx;
    }
  }
  if (best !== null) return best;
  for (let i = wIdx; i < limit; i++) {
    if (s[i] === 'C' && i - wIdx >= 52 && i - wIdx <= 59) return i;
  }
  return null;
}

// J-region anchor: the F of the conserved "FG.G" motif. All three motifs are
// evaluated (as in the Python reference); candidates whose downstream context
// matches the J consensus win, otherwise the most 3' candidate is used.
function findJAnchor(s, start, end) {
  const pats = [/FG[GQSTPEA][GTKQE]/g, /FG[A-Z]/g, /F[GSTPEA][GQSTPEK]/g];
  const jTail = /^[GQ][GQ][GT][KQTV][LVT]/;
  let best = null;
  let bestTail = false;
  for (const pat of pats) {
    pat.lastIndex = 0;
    const window = s.slice(start, end);
    let m;
    while ((m = pat.exec(window)) !== null) {
      const idx = start + m.index;
      const tailOk = jTail.test(s.slice(idx + 2, idx + 8));
      if (best === null || (tailOk && !bestTail) || (tailOk === bestTail && idx > best)) {
        best = idx;
        bestTail = tailOk;
      }
    }
  }
  if (best !== null) return best;
  for (let i = start; i < Math.min(start + 4, end); i++) {
    if (s[i] === 'F') return i;
  }
  return null;
}

// Window used to look for the J-region anchor. The reference uses 16 residues,
// which truncates the "FG.G" motif for long CDR3s.
const J_ANCHOR_WINDOW = 45;
// Kabat VL framework boundaries used to sanity check the alignment.
const FR1_MAX = 23;
const IMGT_CYS = 23;
const IMGT_FR4_START = 118;

/**
 * Number a light chain (kappa/lambda) with Kabat numbering.
 * @param {string} seq
 * @param {{imgt?: Array<{pos:number, ins:string}|null>}} [opts]
 *   `imgt` is the IMGT alignment of `seq` (one entry per residue, `null` when
 *   the residue has no IMGT position). When given, FR1 is numbered from the
 *   alignment (so lambda's L10 gap is preserved), FR2/CDR2/FR3 labels come from
 *   the IMGT states (CDR2 insertions on L52A..), and Phe98 is taken from the
 *   first FR4 residue. The Cys23 entry is only used when it is actually Cys.
 * @returns {{ residues: { pos: string, aa: string, region: string }[], warnings: string[] }}
 */
export function numberLight(seq, opts = {}) {
  const s = seq.toUpperCase().trim();
  const n = s.length;
  const warnings = [];
  const res = [];
  const imgt = Array.isArray(opts.imgt) ? opts.imgt : null;

  // ---- FR1 (1-23), Cys23 anchor ----
  let cIdx = null;
  let cysFromAlignment = false;
  let badAlignedCys = false;
  if (imgt) {
    for (let i = 0; i < Math.min(imgt.length, n); i++) {
      const m = imgt[i];
      if (m && m.pos === IMGT_CYS && !m.ins) {
        if (s[i] === 'C') {
          cIdx = i;
          cysFromAlignment = true;
          break;
        }
        badAlignedCys = true;
      }
    }
  }
  if (cIdx === null) {
    for (let i = 20; i < Math.min(27, n); i++) {
      if (s[i] === 'C') {
        cIdx = i;
        break;
      }
    }
  }
  if (cIdx === null) {
    throw new Error(`[VL] no conserved Cys23 anchor found in first 27 residues: ${s.slice(0, 27)}`);
  }
  if (badAlignedCys && !cysFromAlignment) {
    warnings.push('[VL] IMGT position 23 is not a Cys; using the anchor scan instead');
  }
  if (!cysFromAlignment && cIdx !== 22) {
    warnings.push(`[VL] Cys23 anchor at index ${cIdx + 1} (expected 23); FR1 length adjusted`);
  }
  let fr1Len = cIdx + 1;
  if (fr1Len > FR1_MAX && !(imgt && cysFromAlignment)) {
    // Fallback Cys anchor shifted (e.g. an extra N-terminal residue): keep the
    // standard 23-residue FR1 and let CDR1 absorb the extra residue. Capping
    // avoids double-assigning L24 (FR1 tail vs CDR1 start), which would
    // silently drop one residue. An aligned Cys (including an FR1 insertion
    // before it) is kept in FR1.
    warnings.push(
      `[VL] Cys23 anchor shifted (C at ${cIdx + 1}); FR1 capped to ${FR1_MAX} residues, CDR1 absorbs the extra`,
    );
    fr1Len = FR1_MAX;
  }
  if (imgt && cysFromAlignment) {
    // Number FR1 from the alignment: IMGT state numbers are the Kabat VL
    // positions for states 1-23, so framework gaps (e.g. the missing L10 of
    // lambda light chains) are preserved.
    for (let i = 0; i < fr1Len; i++) {
      const m = imgt[i];
      const pos = m && m.pos ? `L${m.pos}${m.ins || ''}` : `L${i + 1}`;
      res.push({ pos, aa: s[i], region: 'FR1' });
    }
  } else {
    for (let i = 0; i < fr1Len; i++) res.push({ pos: `L${i + 1}`, aa: s[i], region: 'FR1' });
  }

  // ---- CDR1 (24-34 + 27A..) and FR2 anchor Trp35 ----
  const wIdx = findFr2Trp(s, fr1Len, Math.min(fr1Len + 24, n));
  if (wIdx === null) throw new Error('[VL] no Trp35 (FR2 start) anchor found');
  const cdr1Seq = s.slice(fr1Len, wIdx);
  if (!(cdr1Seq.length >= 6 && cdr1Seq.length <= 17)) {
    warnings.push(`[VL] unusual CDR1 length ${cdr1Seq.length} (expected 11-17)`);
  }
  const cdr1Labels = cdr1Positions(cdr1Seq.length);
  for (let j = 0; j < cdr1Seq.length; j++) {
    res.push({ pos: cdr1Labels[j], aa: cdr1Seq[j], region: 'CDR1' });
  }

  // ---- FR2 (35-49, 15) / CDR2 (50-56, 7) / FR3 (57-88, 32) ----
  const c88Idx = findC88(s, wIdx, n);
  if (c88Idx === null) throw new Error('[VL] no Cys88 (FR3 end) anchor found');
  const seg = s.slice(wIdx, c88Idx + 1);
  const segLen = seg.length;
  const diff = segLen - 54;
  if (diff > 4 || diff < -2) warnings.push(`[VL] FR2+CDR2+FR3 span ${segLen} vs 54 expected (diff ${diff > 0 ? '+' : ''}${diff})`);
  const spanLabels = imgt ? kabatSpanLabels(imgt, wIdx, c88Idx) : null;
  if (imgt && !spanLabels) {
    warnings.push('[VL] IMGT alignment does not cover the W35-C88 span; using anchor heuristics');
  }
  if (spanLabels) {
    for (let i = wIdx; i <= c88Idx; i++) {
      const label = spanLabels.get(i);
      const region = label.n <= 49 ? 'FR2' : label.n <= 56 ? 'CDR2' : 'FR3';
      res.push({ pos: `L${label.n}${label.ins || ''}`, aa: s[i], region });
    }
  } else {
    let fr2Len = 15;
    if (diff < 0) fr2Len = Math.max(0, 15 + diff);
    const rem = segLen - fr2Len;
    const cdr2Len = Math.min(7, rem);
    const fr3Len = rem - cdr2Len;
    const fr3Ins = Math.max(0, fr3Len - 32);
    if (fr3Ins > 3) warnings.push(`[VL] ${fr3Ins} FR3 insertions (82A..) - unusual`);
    for (let j = 0; j < fr2Len; j++) res.push({ pos: `L${35 + j}`, aa: seg[j], region: 'FR2' });
    for (let j = 0; j < cdr2Len; j++) res.push({ pos: `L${50 + j}`, aa: seg[fr2Len + j], region: 'CDR2' });
    for (let j = 0; j < fr3Len; j++) {
      let pos;
      if (j < 26) pos = `L${57 + j}`;
      else if (j < 26 + fr3Ins) pos = `L82${letters(j - 25).slice(-1)[0]}`;
      else pos = `L${57 + j - fr3Ins}`;
      res.push({ pos, aa: seg[fr2Len + cdr2Len + j], region: 'FR3' });
    }
  }

  // ---- CDR3 (89-97, inserts on 95) and FR4 (98-107), Phe98 anchor ----
  const jStart = c88Idx + 1;
  let f98Idx = null;
  let fr4StartIdx = null;
  if (imgt) {
    for (let i = jStart; i < n; i++) {
      const m = imgt[i];
      if (m && m.pos >= IMGT_FR4_START) {
        fr4StartIdx = i;
        break;
      }
    }
  }
  if (fr4StartIdx !== null) {
    if (s[fr4StartIdx] !== 'F') {
      throw new Error('[VL] the IMGT FR4 region does not start with Phe98 - the J-region is broken');
    }
    f98Idx = fr4StartIdx;
  } else {
    f98Idx = findJAnchor(s, jStart, Math.min(jStart + J_ANCHOR_WINDOW, n));
  }
  if (f98Idx === null) {
    // The alignment says the domain continues into the FR4 region but no J
    // anchor could be found: the J region is broken (e.g. the conserved
    // Trp/Phe mutated), so numbering FR4 would silently mis-assign it.
    const fr4Residues = imgt
      ? imgt.filter((m) => m && m.pos >= IMGT_FR4_START).length
      : 0;
    if (fr4Residues > 0) {
      throw new Error(
        '[VL] no Phe98 (J-region) anchor found although the alignment places ' +
          `${fr4Residues} residue(s) in the FR4 region - the J region is broken`,
      );
    }
    // V-region-only sequence (no J region): treat the tail as CDR3.
    f98Idx = n;
    warnings.push('[VL] no J-region anchor found (V region only?) - the tail was assigned to CDR3, FR4 is empty');
  }
  const cdr3Seq = s.slice(jStart, f98Idx);
  const cdr3Labels = cdr3Positions(cdr3Seq.length);
  for (let j = 0; j < cdr3Seq.length; j++) {
    res.push({ pos: cdr3Labels[j], aa: cdr3Seq[j], region: 'CDR3' });
  }
  const fr4 = s.slice(f98Idx, f98Idx + 10);
  if (fr4.length < 10) warnings.push(`[VL] FR4 shorter than expected (${fr4.length}/10) - truncated J region?`);
  for (let j = 0; j < fr4.length; j++) res.push({ pos: `L${98 + j}`, aa: fr4[j], region: 'FR4' });
  const tailEnd = f98Idx + 10;
  if (tailEnd < n) warnings.push(`[VL] ${n - tailEnd} trailing residue(s) after FR4 ignored: ${s.slice(tailEnd)}`);

  return { residues: res, warnings };
}
