// Post-numbering sequence analysis: CDR lengths and liability (PTM) motifs.

/** @typedef {{pos:string, aa:string, region:string}} Residue */

/**
 * CDR loop lengths.
 * @param {Record<string,string>} cdrs
 * @returns {{CDR1:number, CDR2:number, CDR3:number}}
 */
export function cdrLengths(cdrs) {
  return {
    CDR1: (cdrs.CDR1 || '').length,
    CDR2: (cdrs.CDR2 || '').length,
    CDR3: (cdrs.CDR3 || '').length,
  };
}

// motif -> { type, note, severity }
const MOTIFS = [
  { re: /N[^P][ST]/g, type: 'N-glycosylation sequon', note: 'N-X-S/T (X≠P)', severity: 'high' },
  { re: /N[GST]/g, type: 'Deamidation hotspot', note: 'Asn followed by Gly/Ser/Thr', severity: 'medium' },
  { re: /D[GST]/g, type: 'Isomerization hotspot', note: 'Asp followed by Gly/Ser/Thr', severity: 'medium' },
  { re: /DP/g, type: 'Acid-labile DP', note: 'Asp-Pro fragmentation', severity: 'medium' },
  { re: /RGD|KGD|NGR/g, type: 'Integrin-binding motif', note: 'RGD/KGD/NGR', severity: 'low' },
];

/**
 * Find positional liabilities across the numbered residues.
 * @param {Residue[]} residues
 * @param {string} [chainType]
 * @returns {Array<{type:string, motif:string, region:string, positions:string, note:string, severity:string}>}
 */
export function findLiabilities(residues, chainType) {
  const out = [];
  if (!residues || residues.length === 0) return out;
  const seq = residues.map((r) => r.aa).join('');

  // Collect CDR spans once.
  const inCdr = residues.map((r) => r.region === 'CDR1' || r.region === 'CDR2' || r.region === 'CDR3');

  const seen = new Set();
  for (const m of MOTIFS) {
    m.re.lastIndex = 0;
    let match;
    while ((match = m.re.exec(seq)) !== null) {
      const start = match.index;
      const end = start + match[0].length - 1;
      const region = residues[start].region;
      const positions = `${residues[start].pos}–${residues[end].pos}`;
      const key = `${m.type}:${positions}`;
      if (seen.has(key)) continue;
      seen.add(key);
      out.push({
        type: m.type,
        motif: match[0],
        region,
        positions,
        note: m.note,
        severity: m.severity,
      });
    }
  }

  // Oxidation-prone residues (Met/Trp), reported individually.
  for (let i = 0; i < residues.length; i++) {
    const aa = residues[i].aa;
    if (aa === 'M' || aa === 'W') {
      out.push({
        type: 'Oxidation risk',
        motif: aa,
        region: residues[i].region,
        positions: residues[i].pos,
        note: `${aa} is susceptible to oxidation`,
        severity: 'low',
      });
    }
  }

  // Free / unpaired cysteines: flag Cys inside CDRs and any beyond the two
  // conserved framework disulfide cysteines.
  const cysIdx = residues.map((r, i) => (r.aa === 'C' ? i : -1)).filter((i) => i >= 0);
  for (const i of cysIdx) {
    if (inCdr[i]) {
      out.push({
        type: 'Free cysteine (CDR)',
        motif: 'C',
        region: residues[i].region,
        positions: residues[i].pos,
        note: 'CDR cysteine is not expected to be disulfide-paired',
        severity: 'medium',
      });
    }
  }
  if (cysIdx.length > 2) {
    const extra = cysIdx.slice(2).map((i) => residues[i].pos).join(', ');
    out.push({
      type: 'Extra cysteines',
      motif: `${cysIdx.length}×C`,
      region: residues[2] ? residues[2].region : '',
      positions: extra,
      note: 'More than the two conserved framework cysteines',
      severity: 'medium',
    });
  }

  return out;
}
