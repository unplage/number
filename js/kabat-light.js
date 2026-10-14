// Kabat light-chain numbering — port of number_light from
// unplage/humanization scripts/humanize/numbering.py (MIT License).
// Anchor-based engine: Cys23 -> Trp35 -> Cys88 -> Phe98.
// Used because ANARCI/AbRSA are known to mis-number VL CDR3.

function letters(n) {
  const out = [];
  for (let i = 0; i < n; i++) {
    const a = Math.floor(i / 26);
    const b = i % 26;
    out.push(a === 0 ? String.fromCharCode(65 + b) : String.fromCharCode(64 + a) + String.fromCharCode(65 + b));
  }
  return out;
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

function findJAnchor(s, start, end) {
  const pats = [/FG[GQSTPEA][GTKQE]/g, /FG[A-Z]/g, /F[GSTPEA][GQSTPEK]/g];
  const jTail = /^[GQ][GT][QKE][GT][KLTV]/;
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
    if (best !== null) break;
  }
  if (best !== null) return best;
  for (let i = start; i < Math.min(start + 4, end); i++) {
    if (s[i] === 'F') return i;
  }
  return null;
}

/**
 * Number a light chain (kappa/lambda) with Kabat numbering.
 * @param {string} seq
 * @returns {{ residues: { pos: string, aa: string, region: string }[], warnings: string[] }}
 */
export function numberLight(seq) {
  const s = seq.toUpperCase().trim();
  const n = s.length;
  const warnings = [];
  const res = [];

  // FR1 (1-23), Cys23 anchor
  let cIdx = null;
  for (let i = 20; i < Math.min(27, n); i++) {
    if (s[i] === 'C') {
      cIdx = i;
      break;
    }
  }
  if (cIdx === null) {
    throw new Error(`[VL] no conserved Cys23 anchor found in first 27 residues: ${s.slice(0, 27)}`);
  }
  if (cIdx !== 22) warnings.push(`[VL] Cys23 anchor at index ${cIdx + 1} (expected 23); FR1 length adjusted`);
  const fr1Len = cIdx + 1;
  for (let i = 0; i < fr1Len; i++) res.push({ pos: `L${i + 1}`, aa: s[i], region: 'FR1' });

  // CDR1 (24-34 + 27A..) and FR2 anchor Trp35
  const wIdx = findFr2Trp(s, fr1Len, Math.min(fr1Len + 24, n));
  if (wIdx === null) throw new Error('[VL] no Trp35 (FR2 start) anchor found');
  const cdr1Seq = s.slice(fr1Len, wIdx);
  if (!(cdr1Seq.length >= 6 && cdr1Seq.length <= 17)) {
    warnings.push(`[VL] unusual CDR1 length ${cdr1Seq.length} (expected 11-17)`);
  }
  const cdr1Ins = Math.max(0, cdr1Seq.length - 11);
  for (let j = 0; j < cdr1Seq.length; j++) {
    let pos;
    if (j < 4) pos = `L${24 + j}`;
    else if (j < 4 + cdr1Ins) pos = `L27${letters(j - 3).slice(-1)[0]}`;
    else pos = `L${28 + j - 4 - cdr1Ins}`;
    res.push({ pos, aa: cdr1Seq[j], region: 'CDR1' });
  }

  // FR2 (35-49, 15) / CDR2 (50-56, 7) / FR3 (57-88, 32)
  const c88Idx = findC88(s, wIdx, n);
  if (c88Idx === null) throw new Error('[VL] no Cys88 (FR3 end) anchor found');
  const seg = s.slice(wIdx, c88Idx + 1);
  const segLen = seg.length;
  const diff = segLen - 54;
  if (diff > 4 || diff < -2) warnings.push(`[VL] FR2+CDR2+FR3 span ${segLen} vs 54 expected (diff ${diff > 0 ? '+' : ''}${diff})`);
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

  // CDR3 (89-97, inserts on 95) and FR4 (98-107), Phe98 anchor
  const jStart = c88Idx + 1;
  let f98Idx = findJAnchor(s, jStart, Math.min(jStart + 16, n));
  if (f98Idx === null) f98Idx = n;
  const cdr3Seq = s.slice(jStart, f98Idx);
  const cdr3Ins = Math.max(0, cdr3Seq.length - 9);
  for (let j = 0; j < cdr3Seq.length; j++) {
    let pos;
    if (j < 7) pos = `L${89 + j}`;
    else if (j < 7 + cdr3Ins) pos = `L95${letters(j - 6).slice(-1)[0]}`;
    else pos = `L${89 + j - cdr3Ins}`;
    res.push({ pos, aa: cdr3Seq[j], region: 'CDR3' });
  }
  const fr4 = s.slice(f98Idx, f98Idx + 10);
  if (fr4.length < 10) warnings.push(`[VL] FR4 shorter than expected (${fr4.length}/10) - truncated J region?`);
  for (let j = 0; j < fr4.length; j++) res.push({ pos: `L${98 + j}`, aa: fr4[j], region: 'FR4' });
  const tailEnd = f98Idx + 10;
  if (tailEnd < n) warnings.push(`[VL] ${n - tailEnd} trailing residue(s) after FR4 ignored: ${s.slice(tailEnd)}`);

  return { residues: res, warnings };
}
