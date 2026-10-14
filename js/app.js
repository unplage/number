// UI controller.

import { parseFasta } from './fasta.js';
import { analyze, getEngineVersion, SCHEMES } from './engine.js';
import { toTsv, toCsv, toJson, toAligned, toAnarci, toFastaWithCdr, download } from './export.js';

const $ = (sel) => document.querySelector(sel);
const el = (tag, cls, text) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (text != null) n.textContent = text;
  return n;
};

const VH = 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS';
const VL = 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK';

const EXAMPLES = {
  vh: `>trastuzumab_VH\n${VH}\n`,
  vl: `>trastuzumab_VL_kappa\n${VL}\n`,
  lambda: '>lambda_VL\nQSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL\n',
  vhh: '>camelid_VHH\nQVQLQESGGGLVQAGGSLRLSCAASGRTFSDYVIGWFRQAPGKEREFVAAITWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAAGLGTVVAYWGQGTQVTVSS\n',
  scfv: `>scFv_example\n${VH}GGGGSGGGGSGGGGS${VL}\n`,
  longcdr3: '>VH_long_CDR3\nQVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARRRRRRRRRRRSSSSSSSSSYYYYYYYYYWGQGTLVTVSS\n',
};

const SCHEME_DESC = {
  imgt: 'IMGT unique numbering (Lefranc 2003). ANARCI-identical, via the anarcism WebAssembly engine. Region labels come straight from the engine.',
  kabat: 'Kabat. VH uses the official ANARCI Kabat rules; VL uses a dedicated anchor engine (Cys23→Trp35→Cys88→Phe98) that avoids ANARCI/AbRSA VL CDR3 errors.',
  chothia: 'Chothia structural numbering, ported from ANARCI.',
  martin: 'Martin (extended Chothia) numbering, ported from ANARCI.',
  aho: 'AHo numbering (Honegger & Plückthun), ported from ANARCI.',
  wolfguy: 'Wolfguy numbering (Bujotzek et al.), ported from ANARCI. Positions encode chain and region.',
};

const state = { results: null, records: null, scheme: 'imgt', engineVersion: null, controller: null, compare: null };

function setStatus(html, isError) {
  const s = $('#status');
  s.innerHTML = html;
  s.classList.toggle('err', Boolean(isError));
}

function currentScheme() {
  return $('#scheme-select').value;
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

function chip(text, cls) {
  return el('span', `chip ${cls || ''}`, text);
}

function renderResidueStrip(residues) {
  const strip = el('div', 'residue-strip');
  for (const r of residues) {
    const cell = el('div', `res ${r.region || ''}`);
    cell.append(el('span', 'pos', r.pos), el('span', 'aa', r.aa));
    strip.appendChild(cell);
  }
  return strip;
}

function renderTable(residues) {
  const table = el('table', 'numb');
  const thead = el('thead');
  const hr = el('tr');
  for (const h of ['Position', 'AA', 'Region']) hr.appendChild(el('th', null, h));
  thead.appendChild(hr);
  table.appendChild(thead);
  const tbody = el('tbody');
  for (const r of residues) {
    const tr = el('tr');
    const td = el('td', null, r.pos);
    const td2 = el('td', null, r.aa);
    const td3 = el('td', (r.region || '').toLowerCase(), r.region || '');
    tr.append(td, td2, td3);
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);
  const scroll = el('div');
  scroll.style.maxHeight = '420px';
  scroll.style.overflow = 'auto';
  scroll.appendChild(table);
  return scroll;
}

function renderAligned(residues) {
  const lines = [];
  const W = 7;
  let num = [];
  let aa = [];
  for (let i = 0; i < residues.length; i++) {
    num.push(residues[i].pos.padEnd(W));
    aa.push(residues[i].aa.padEnd(W));
    if ((i + 1) % 10 === 0 || i === residues.length - 1) {
      lines.push(num.join(' ').trimEnd());
      lines.push(aa.join(' ').trimEnd());
      lines.push('');
      num = [];
      aa = [];
    }
  }
  return el('pre', 'aligned', lines.join('\n'));
}

function renderCdrBoxes(d) {
  const wrap = el('div', 'cdr-boxes');
  for (const name of ['CDR1', 'CDR2', 'CDR3']) {
    const box = el('div', 'cdr-box');
    const len = d.cdrLengths ? d.cdrLengths[name] : (d.cdrs[name] || '').length;
    box.append(el('div', 'label', `${name} · ${len} aa`), el('div', 'val', d.cdrs[name] || '—'));
    wrap.appendChild(box);
  }
  return wrap;
}

function renderLiabilities(d) {
  if (!d.liabilities || d.liabilities.length === 0) return null;
  const wrap = el('div');
  wrap.appendChild(el('div', 'label', 'Sequence liabilities'));
  const row = el('div');
  for (const l of d.liabilities) {
    const pill = el('span', 'liability');
    pill.append(el('span', `sev sev-${l.severity}`), el('span', null, `${l.type} · ${l.positions} (${l.motif})`));
    row.appendChild(pill);
  }
  wrap.appendChild(row);
  return wrap;
}

function domainCopyText(d, seqId, scheme) {
  const lines = [`>${seqId} ${d.chainType} scheme=${scheme}`];
  lines.push(`# CDR1=${d.cdrs.CDR1} CDR2=${d.cdrs.CDR2} CDR3=${d.cdrs.CDR3}`);
  lines.push(d.residues.map((r) => `${r.pos}\t${r.aa}\t${r.region}`).join('\n'));
  return lines.join('\n') + '\n';
}

function renderDomain(seq, d, seqIndex, domIndex) {
  const wrap = el('div', 'domain');

  const head = el('div', 'domain-head');
  head.appendChild(chip(d.chainLabel, 'chain'));
  head.appendChild(chip(d.species || 'unknown', 'species'));
  head.appendChild(chip(`${d.confidence} confidence`, `conf-${d.confidence}`));
  if (typeof d.bitScore === 'number') head.appendChild(chip(`bit ${d.bitScore.toFixed(1)}`, ''));
  if (d.isVhh) head.appendChild(chip(`VHH (${d.vhhPositions.join(', ')})`, 'vhh'));
  if (d.germline && d.germline.vGene) {
    const id = d.germline.vIdentity != null ? ` ${(d.germline.vIdentity * 100).toFixed(0)}%` : '';
    head.appendChild(chip(`${d.germline.vGene}${id}${d.germline.jGene ? ` / ${d.germline.jGene}` : ''}`, 'germline'));
  }
  wrap.appendChild(head);

  if (d.errors && d.errors.length) {
    const n = el('div', 'notice err', d.errors.join(' '));
    wrap.appendChild(n);
  }

  if (d.residues && d.residues.length) {
    const tabs = el('div', 'view-tabs');
    const body = el('div');
    const views = {
      Strip: renderResidueStrip(d.residues),
      Table: renderTable(d.residues),
      Aligned: renderAligned(d.residues),
    };
    let active = 'Strip';
    const buttons = {};
    const show = (name) => {
      active = name;
      for (const [k, v] of Object.entries(views)) v.style.display = k === name ? '' : 'none';
      for (const [k, b] of Object.entries(buttons)) b.classList.toggle('active', k === name);
    };
    for (const name of Object.keys(views)) {
      const b = el('button', null, name);
      b.type = 'button';
      b.addEventListener('click', () => show(name));
      buttons[name] = b;
      tabs.appendChild(b);
    }
    for (const v of Object.values(views)) body.appendChild(v);
    wrap.append(tabs, body);
    show(active);

    wrap.appendChild(renderCdrBoxes(d));

    const legend = el('div', 'legend');
    legend.innerHTML = '<span class="l-fr">Framework</span><span class="l-cdr1">CDR1</span><span class="l-cdr2">CDR2</span><span class="l-cdr3">CDR3</span>';
    wrap.appendChild(legend);

    const liab = renderLiabilities(d);
    if (liab) wrap.appendChild(liab);

    const actions = el('div', 'mini-actions');
    const copyBtn = el('button', 'ghost-btn', 'Copy');
    copyBtn.type = 'button';
    copyBtn.addEventListener('click', async () => {
      try {
        await navigator.clipboard.writeText(domainCopyText(d, seq.id, seq.scheme));
        copyBtn.textContent = 'Copied';
        setTimeout(() => (copyBtn.textContent = 'Copy'), 1200);
      } catch {
        setStatus('Clipboard unavailable.', true);
      }
    });
    const dlBtn = el('button', 'ghost-btn', 'Download domain');
    dlBtn.type = 'button';
    dlBtn.addEventListener('click', () => {
      download(`${seq.id}_${d.chainType}.tsv`, domainCopyText(d, seq.id, seq.scheme).replace(/\t/g, '\t'), 'text/tab-separated-values;charset=utf-8');
    });
    actions.append(copyBtn, dlBtn);
    wrap.appendChild(actions);
  }

  if (d.warnings && d.warnings.length) {
    wrap.appendChild(el('div', 'notice warn', 'Warnings: ' + d.warnings.join(' ')));
  }
  return wrap;
}

function renderResults(results) {
  const container = $('#results');
  container.innerHTML = '';
  container.hidden = false;
  $('#downloads').hidden = false;

  for (const seq of results) {
    const card = el('article', 'seq-card');
    const head = el('div', 'seq-head');
    head.appendChild(el('span', 'sid', '>' + seq.id));
    head.appendChild(el('span', 'meta', `${seq.sequence.length} aa · ${seq.scheme.toUpperCase()}${seq.description ? ` · ${seq.description}` : ''}`));
    card.appendChild(head);

    if (seq.errors && seq.errors.length) {
      const w = el('div', 'domain');
      w.appendChild(el('div', 'notice err', seq.errors.join(' ')));
      card.appendChild(w);
      container.appendChild(card);
      continue;
    }
    seq.domains.forEach((d, di) => card.appendChild(renderDomain(seq, d, 0, di)));
    container.appendChild(card);
  }
  if (state.engineVersion) {
    $('#engine-info').textContent = `${results.length} sequence(s) · engine anarcism v${state.engineVersion} · generated ${new Date().toLocaleString()}`;
  }
}

// ---------------------------------------------------------------------------
// Compare IMGT vs Kabat
// ---------------------------------------------------------------------------

async function renderCompare() {
  const text = $('#fasta-input').value.trim();
  if (!text) {
    setStatus('Please enter sequences first.', true);
    return;
  }
  const records = parseFasta(text);
  if (records.length === 0) {
    setStatus('No valid sequences found.', true);
    return;
  }
  setBusy(true);
  setStatus('<span class="spinner"></span>Numbering with IMGT and Kabat for comparison…');
  try {
    const [imgt, kabat] = await Promise.all([
      analyze(records, 'imgt'),
      analyze(records, 'kabat'),
    ]);
    state.engineVersion = imgt.engineVersion;
    renderCompareResults(imgt.results, kabat.results);
    setStatus('Comparison ready.');
  } catch (e) {
    setStatus('Comparison failed: ' + e.message, true);
  } finally {
    setBusy(false);
  }
}

function renderCompareResults(imgtResults, kabatResults) {
  const container = $('#results');
  container.innerHTML = '';
  container.hidden = false;
  for (let i = 0; i < imgtResults.length; i++) {
    const a = imgtResults[i];
    const b = kabatResults[i];
    const card = el('article', 'seq-card');
    const head = el('div', 'seq-head');
    head.appendChild(el('span', 'sid', '>' + a.id));
    head.appendChild(el('span', 'meta', 'IMGT vs Kabat'));
    card.appendChild(head);
    const body = el('div', 'domain');
    const da = a.domains && a.domains[0];
    const db = b.domains && b.domains[0];
    if (!da || !db) {
      body.appendChild(el('div', 'notice err', (a.errors || []).concat(b.errors || []).join(' ') || 'Could not compare.'));
    } else {
      const table = el('table', 'numb');
      const thead = el('thead');
      const hr = el('tr');
      for (const h of ['#', 'IMGT', 'Kabat']) hr.appendChild(el('th', null, h));
      thead.appendChild(hr);
      table.appendChild(thead);
      const tbody = el('tbody');
      const n = Math.max(da.residues.length, db.residues.length);
      for (let j = 0; j < n; j++) {
        const ra = da.residues[j];
        const rb = db.residues[j];
        const tr = el('tr');
        tr.appendChild(el('td', null, String(j + 1)));
        tr.appendChild(el('td', (ra && ra.region ? ra.region.toLowerCase() : ''), ra ? `${ra.pos} ${ra.aa}` : '—'));
        tr.appendChild(el('td', (rb && rb.region ? rb.region.toLowerCase() : ''), rb ? `${rb.pos} ${rb.aa}` : '—'));
        tbody.appendChild(tr);
      }
      table.appendChild(tbody);
      const scroll = el('div');
      scroll.style.maxHeight = '480px';
      scroll.style.overflow = 'auto';
      scroll.appendChild(table);
      body.appendChild(scroll);
    }
    card.appendChild(body);
    container.appendChild(card);
  }
  $('#results').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

// ---------------------------------------------------------------------------
// Analyze
// ---------------------------------------------------------------------------

function setBusy(busy) {
  $('#analyze-btn').disabled = busy;
  $('#compare-btn').disabled = busy;
  $('#cancel-btn').hidden = !busy;
  if (!busy) $('#progress').hidden = true;
}

async function run() {
  const text = $('#fasta-input').value.trim();
  if (!text) {
    setStatus('Please paste or upload sequences first.', true);
    return;
  }
  const records = parseFasta(text);
  if (records.length === 0) {
    setStatus('No valid sequences found in the input.', true);
    return;
  }

  state.records = records;
  state.scheme = currentScheme();
  const controller = new AbortController();
  state.controller = controller;
  setBusy(true);
  $('#progress').hidden = false;
  $('#progress-bar').style.width = '0%';
  setStatus('<span class="spinner"></span>Loading engine and numbering…');

  try {
    const { results, engineVersion } = await analyze(records, state.scheme, {
      signal: controller.signal,
      onProgress: (done, total) => {
        $('#progress-bar').style.width = total ? `${Math.round((done / total) * 100)}%` : '0%';
      },
    });
    state.results = results;
    state.engineVersion = engineVersion;
    renderResults(results);
    const failed = results.filter((r) => r.errors && r.errors.length).length;
    setStatus(`Done. ${results.length - failed}/${results.length} sequence(s) numbered.`);
    $('#results').scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (e) {
    if (controller.signal.aborted) setStatus('Cancelled.');
    else {
      console.error(e);
      setStatus('Analysis failed: ' + (e && e.message ? e.message : String(e)), true);
    }
  } finally {
    setBusy(false);
    state.controller = null;
  }
}

// ---------------------------------------------------------------------------
// Event wiring
// ---------------------------------------------------------------------------

$('#file-input').addEventListener('change', async (ev) => {
  const file = ev.target.files && ev.target.files[0];
  if (!file) return;
  $('#fasta-input').value = await file.text();
  setStatus(`Loaded ${file.name}.`);
});

const dz = $('#drop-zone');
dz.addEventListener('dragover', (e) => {
  e.preventDefault();
  dz.classList.add('dragover');
});
dz.addEventListener('dragleave', () => dz.classList.remove('dragover'));
dz.addEventListener('drop', async (e) => {
  e.preventDefault();
  dz.classList.remove('dragover');
  const file = e.dataTransfer.files && e.dataTransfer.files[0];
  if (file) {
    $('#fasta-input').value = await file.text();
    setStatus(`Loaded ${file.name}.`);
  }
});

$('#example-select').addEventListener('change', (e) => {
  const key = e.target.value;
  if (key && EXAMPLES[key]) {
    $('#fasta-input').value = EXAMPLES[key];
    setStatus('Example loaded. Click Analyze.');
  }
  e.target.value = '';
});

$('#clear-btn').addEventListener('click', () => {
  $('#fasta-input').value = '';
  $('#results').innerHTML = '';
  $('#results').hidden = true;
  $('#downloads').hidden = true;
  state.results = null;
  setStatus('');
});

$('#scheme-select').addEventListener('change', () => {
  $('#scheme-desc').textContent = SCHEME_DESC[currentScheme()];
});

$('#analyze-btn').addEventListener('click', run);
$('#compare-btn').addEventListener('click', renderCompare);
$('#cancel-btn').addEventListener('click', () => {
  if (state.controller) state.controller.abort();
});

const DL = { tsv: toTsv, csv: toCsv, json: toJson, aligned: toAligned, anarci: toAnarci, fasta: toFastaWithCdr };
const MIME = { json: 'application/json;charset=utf-8', csv: 'text/csv;charset=utf-8', tsv: 'text/tab-separated-values;charset=utf-8' };
for (const btn of document.querySelectorAll('.dl')) {
  btn.addEventListener('click', () => {
    if (!state.results) return;
    const fmt = /** @type {HTMLElement} */ (btn).dataset.fmt;
    const ext = { aligned: 'txt', anarci: 'txt', fasta: 'fasta' }[fmt] || fmt;
    download(`antibody_numbering.${ext}`, DL[fmt](state.results), MIME[fmt] || 'text/plain;charset=utf-8');
  });
}

$('#theme-btn').addEventListener('click', () => {
  const cur = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark';
  document.documentElement.dataset.theme = cur;
});

// Service worker + offline badge.
if ('serviceWorker' in navigator) {
  navigator.serviceWorker.register('sw.js').then(() => {
    if (navigator.serviceWorker.controller || navigator.onLine) $('#offline-badge').hidden = false;
  }).catch(() => {});
}
window.addEventListener('online', () => { $('#offline-badge').hidden = false; });

$('#scheme-desc').textContent = SCHEME_DESC.imgt;
