// Result export: TSV, JSON, aligned numbered text.

function flattenRows(results) {
  const rows = [];
  for (const r of results) {
    if (r.errors && r.errors.length) {
      rows.push({
        seq_id: r.id,
        domain: '',
        chain: '',
        species: '',
        scheme: r.scheme,
        position: '',
        amino_acid: '',
        region: '',
        note: r.errors.join('; '),
      });
      continue;
    }
    for (const d of r.domains || []) {
      for (const res of d.residues) {
        rows.push({
          seq_id: r.id,
          domain: String(d.start) + '-' + String(d.end),
          chain: d.chainType,
          species: d.species,
          scheme: r.scheme,
          position: res.pos,
          amino_acid: res.aa,
          region: res.region,
          note: '',
        });
      }
    }
  }
  return rows;
}

export function toTsv(results) {
  const rows = flattenRows(results);
  const cols = ['seq_id', 'domain', 'chain', 'species', 'scheme', 'position', 'amino_acid', 'region', 'note'];
  const lines = [cols.join('\t')];
  for (const row of rows) lines.push(cols.map((c) => row[c]).join('\t'));
  return lines.join('\n') + '\n';
}

export function toJson(results) {
  return JSON.stringify(
    {
      generated: new Date().toISOString(),
      tool: 'antibody-number (anarcism + ANARCI-Kabat port)',
      results,
    },
    null,
    2,
  ) + '\n';
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
      blocks.push(
        `# chain=${d.chainType} species=${d.species} range=${d.start}-${d.end}` +
          (d.isVhh ? ' VHH' : '') +
          (d.warnings && d.warnings.length ? ` warnings=${d.warnings.join('; ')}` : ''),
      );
      // Numbered blocks of 10
      const res = d.residues;
      const POS_W = 7;
      const numLine = [];
      const aaLine = [];
      for (let i = 0; i < res.length; i++) {
        numLine.push(res[i].pos.padEnd(POS_W));
        aaLine.push(res[i].aa.padEnd(POS_W));
        if ((i + 1) % 10 === 0 || i === res.length - 1) {
          blocks.push(numLine.join(' '));
          blocks.push(aaLine.join(' '));
          blocks.push('');
          numLine.length = 0;
          aaLine.length = 0;
        }
      }
      const cdrs = Object.entries(d.cdrs || {})
        .filter(([, v]) => v)
        .map(([k, v]) => `${k}=${v}`)
        .join('  ');
      if (cdrs) blocks.push(`# ${cdrs}`);
    }
    blocks.push('');
  }
  return blocks.join('\n') + '\n';
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
