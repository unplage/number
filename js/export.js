// Result export: TSV, CSV, JSON, aligned text, ANARCI-style, annotated FASTA.

function rows(results) {
  const out = [];
  for (const r of results) {
    if (r.errors && r.errors.length) {
      out.push({ seq_id: r.id, chain: '', species: '', scheme: r.scheme, position: '', amino_acid: '', region: '', note: r.errors.join('; ') });
      continue;
    }
    for (const d of r.domains || []) {
      const g = d.germline || {};
      for (const res of d.residues) {
        out.push({
          seq_id: r.id,
          chain: d.chainType,
          species: d.species || '',
          v_gene: g.vGene || '',
          v_identity: g.vIdentity != null ? g.vIdentity.toFixed(3) : '',
          j_gene: g.jGene || '',
          bit_score: d.bitScore != null ? d.bitScore.toFixed(1) : '',
          confidence: d.confidence || '',
          scheme: r.scheme,
          position: res.pos,
          amino_acid: res.aa,
          region: res.region,
          note: '',
        });
      }
    }
  }
  return out;
}

const COLUMNS = [
  'seq_id', 'chain', 'species', 'v_gene', 'v_identity', 'j_gene',
  'bit_score', 'confidence', 'scheme', 'position', 'amino_acid', 'region', 'note',
];

export function toTsv(results) {
  const lines = [COLUMNS.join('\t')];
  for (const row of rows(results)) lines.push(COLUMNS.map((c) => row[c] ?? '').join('\t'));
  return lines.join('\n') + '\n';
}

function csvCell(v) {
  const s = String(v ?? '');
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

export function toCsv(results) {
  const lines = [COLUMNS.join(',')];
  for (const row of rows(results)) lines.push(COLUMNS.map((c) => csvCell(row[c])).join(','));
  return lines.join('\n') + '\n';
}

export function toJson(results, meta = {}) {
  return JSON.stringify({ generated: new Date().toISOString(), tool: 'antibody-number', ...meta, results }, null, 2) + '\n';
}

export function toAligned(results) {
  const blocks = [];
  for (const r of results) {
    blocks.push(`>${r.id}  scheme=${r.scheme}`);
    if (r.errors && r.errors.length) {
      blocks.push(`# ERROR: ${r.errors.join('; ')}`);
      continue;
    }
    for (const d of r.domains || []) {
      const tags = [
        `chain=${d.chainType}`, `species=${d.species}`, `range=${d.start}-${d.end}`,
        d.confidence ? `confidence=${d.confidence}` : '',
        d.isVhh ? 'VHH' : '',
      ].filter(Boolean).join(' ');
      blocks.push(`# ${tags}`);
      if (d.germline && d.germline.vGene) {
        blocks.push(`# germline V=${d.germline.vGene}${d.germline.vIdentity != null ? ` (${(d.germline.vIdentity * 100).toFixed(0)}%)` : ''}${d.germline.jGene ? ` J=${d.germline.jGene}` : ''}`);
      }
      const res = d.residues;
      const W = 7;
      let numLine = [];
      let aaLine = [];
      for (let i = 0; i < res.length; i++) {
        numLine.push(res[i].pos.padEnd(W));
        aaLine.push(res[i].aa.padEnd(W));
        if ((i + 1) % 10 === 0 || i === res.length - 1) {
          blocks.push(numLine.join(' ').trimEnd());
          blocks.push(aaLine.join(' ').trimEnd());
          blocks.push('');
          numLine = [];
          aaLine = [];
        }
      }
      const cdrs = Object.entries(d.cdrs || {}).filter(([, v]) => v).map(([k, v]) => `${k}=${v}`).join('  ');
      if (cdrs) blocks.push(`# ${cdrs}`);
      if (d.liabilities && d.liabilities.length) {
        blocks.push(`# liabilities: ${d.liabilities.map((l) => `${l.type}(${l.positions})`).join('; ')}`);
      }
    }
    blocks.push('');
  }
  return blocks.join('\n') + '\n';
}

/** ANARCI-style numbered output ("//" separated domain records). */
export function toAnarci(results) {
  const blocks = [];
  for (const r of results) {
    if (r.errors && r.errors.length) continue;
    for (const d of r.domains || []) {
      const g = d.germline || {};
      blocks.push(`# ${r.id} | chain=${d.chainType} | species=${d.species} | scheme=${r.scheme}${g.vGene ? ` | v_gene=${g.vGene}` : ''}`);
      blocks.push('# pos\taa\tregion');
      for (const res of d.residues) blocks.push(`${res.pos}\t${res.aa}\t${res.region}`);
      blocks.push('//');
    }
  }
  return blocks.join('\n') + '\n';
}

/** FASTA with CDR annotations in the header, one record per variable domain. */
export function toFastaWithCdr(results) {
  const blocks = [];
  for (const r of results) {
    if (r.errors && r.errors.length) continue;
    for (let di = 0; di < (r.domains || []).length; di++) {
      const d = r.domains[di];
      const g = d.germline || {};
      const header = [
        `${r.id}_${d.chainType}${di > 0 ? `_${di + 1}` : ''}`,
        `scheme=${r.scheme}`,
        `chain=${d.chainType}`,
        g.vGene ? `v_gene=${g.vGene}` : '',
        `CDR1=${d.cdrs.CDR1 || ''}`,
        `CDR2=${d.cdrs.CDR2 || ''}`,
        `CDR3=${d.cdrs.CDR3 || ''}`,
      ].filter(Boolean).join(' | ');
      blocks.push(`>${header}`);
      blocks.push(d.residues.map((x) => x.aa).join(''));
    }
  }
  return blocks.join('\n') + '\n';
}

/** Compact CDR overview: one row per domain (id, CDR1, CDR2, CDR3), tab-separated. */
export function toCdrSummary(results) {
  const lines = ['id\tCDR1\tCDR2\tCDR3'];
  for (const r of results) {
    if (r.errors && r.errors.length) {
      lines.push(`${r.id}\t\t\t`);
      continue;
    }
    const domains = r.domains || [];
    const multi = domains.length > 1;
    for (const d of domains) {
      const id = multi ? `${r.id} (${d.chainType})` : r.id;
      lines.push([id, d.cdrs.CDR1 || '', d.cdrs.CDR2 || '', d.cdrs.CDR3 || ''].join('\t'));
    }
  }
  return lines.join('\n') + '\n';
}

export function download(filename, text, mime) {
  const blob = new Blob([text], { type: mime || 'text/plain;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
