// UI controller for the antibody numbering PWA.

import { parseFasta } from './fasta.js';
import { analyze } from './engine.js';
import { toTsv, toJson, toAligned, download } from './export.js';

const $ = (sel) => document.querySelector(sel);

const els = {
  fasta: $('#fasta-input'),
  file: $('#file-input'),
  example: $('#example-btn'),
  clear: $('#clear-btn'),
  analyze: $('#analyze-btn'),
  status: $('#status'),
  results: $('#results'),
  downloads: $('#downloads'),
  dlTsv: $('#dl-tsv'),
  dlJson: $('#dl-json'),
  dlAligned: $('#dl-aligned'),
  offlineBadge: $('#offline-badge'),
};

const EXAMPLE = `>trastuzumab_VH
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
>trastuzumab_VL
DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK`;

let lastResults = null;

function setStatus(html, isError) {
  els.status.innerHTML = html;
  els.status.classList.toggle('err', Boolean(isError));
}

function currentScheme() {
  return document.querySelector('input[name="scheme"]:checked').value;
}

function renderResults(results) {
  els.results.innerHTML = '';
  els.results.hidden = false;
  els.downloads.hidden = false;

  for (const seq of results) {
    const card = document.createElement('article');
    card.className = 'seq-card';

    const head = document.createElement('div');
    head.className = 'seq-head';
    const sid = document.createElement('span');
    sid.className = 'sid';
    sid.textContent = '>' + seq.id;
    const meta = document.createElement('span');
    meta.className = 'meta';
    meta.textContent = `${seq.sequence.length} aa · ${seq.scheme.toUpperCase()}` + (seq.description ? ` · ${seq.description}` : '');
    head.append(sid, meta);
    card.appendChild(head);

    if (seq.errors && seq.errors.length) {
      const n = document.createElement('div');
      n.className = 'notice err';
      n.textContent = seq.errors.join(' ');
      const wrap = document.createElement('div');
      wrap.className = 'domain';
      wrap.appendChild(n);
      card.appendChild(wrap);
      els.results.appendChild(card);
      continue;
    }

    for (const d of seq.domains) {
      card.appendChild(renderDomain(d));
    }
    els.results.appendChild(card);
  }
}

function renderDomain(d) {
  const wrap = document.createElement('div');
  wrap.className = 'domain';

  const head = document.createElement('div');
  head.className = 'domain-head';
  const chain = document.createElement('span');
  chain.className = 'chip chain';
  chain.textContent = d.chainLabel;
  const species = document.createElement('span');
  species.className = 'chip species';
  species.textContent = d.species;
  head.append(chain, species);
  if (d.isVhh) {
    const vhh = document.createElement('span');
    vhh.className = 'chip vhh';
    vhh.textContent = `VHH hallmark (${d.vhhPositions.join(', ')})`;
    head.appendChild(vhh);
  }
  if (typeof d.bitScore === 'number') {
    const score = document.createElement('span');
    score.className = 'chip score';
    score.textContent = `bit score ${d.bitScore.toFixed(1)}`;
    head.appendChild(score);
  }
  wrap.appendChild(head);

  if (!d.residues || d.residues.length === 0) {
    const n = document.createElement('div');
    n.className = 'notice warn';
    n.textContent = 'No residues numbered for this domain. ' + (d.warnings || []).join(' ');
    wrap.appendChild(n);
    return wrap;
  }

  const strip = document.createElement('div');
  strip.className = 'residue-strip';
  for (const r of d.residues) {
    const cell = document.createElement('div');
    cell.className = 'res ' + (r.region || '');
    const pos = document.createElement('span');
    pos.className = 'pos';
    pos.textContent = r.pos;
    const aa = document.createElement('span');
    aa.className = 'aa';
    aa.textContent = r.aa;
    cell.append(pos, aa);
    strip.appendChild(cell);
  }
  wrap.appendChild(strip);

  const cdrBoxes = document.createElement('div');
  cdrBoxes.className = 'cdr-boxes';
  for (const [name, val] of Object.entries(d.cdrs || {})) {
    const box = document.createElement('div');
    box.className = 'cdr-box';
    const label = document.createElement('div');
    label.className = 'label';
    label.textContent = name;
    const value = document.createElement('div');
    value.className = 'val';
    value.textContent = val || '—';
    box.append(label, value);
    cdrBoxes.appendChild(box);
  }
  wrap.appendChild(cdrBoxes);

  const legend = document.createElement('div');
  legend.className = 'legend';
  legend.innerHTML = '<span class="l-fr">Framework</span><span class="l-cdr1">CDR1</span><span class="l-cdr2">CDR2</span><span class="l-cdr3">CDR3</span>';
  wrap.appendChild(legend);

  if (d.warnings && d.warnings.length) {
    const n = document.createElement('div');
    n.className = 'notice warn';
    n.textContent = 'Warnings: ' + d.warnings.join(' ');
    wrap.appendChild(n);
  }
  return wrap;
}

async function run() {
  const text = els.fasta.value.trim();
  if (!text) {
    setStatus('Please paste or upload sequences first.', true);
    return;
  }
  const records = parseFasta(text);
  if (records.length === 0) {
    setStatus('No valid sequences found in the input.', true);
    return;
  }

  els.analyze.disabled = true;
  setStatus(`<span class="spinner"></span>Loading numbering engine and analyzing ${records.length} sequence(s)…`);

  try {
    const results = await analyze(records, currentScheme());
    lastResults = results;
    renderResults(results);
    const failed = results.filter((r) => r.errors && r.errors.length).length;
    setStatus(`Done. ${results.length - failed}/${results.length} sequence(s) numbered.`);
    els.results.scrollIntoView({ behavior: 'smooth', block: 'start' });
  } catch (e) {
    console.error(e);
    setStatus('Analysis failed: ' + (e && e.message ? e.message : String(e)), true);
  } finally {
    els.analyze.disabled = false;
  }
}

els.file.addEventListener('change', async (ev) => {
  const file = ev.target.files && ev.target.files[0];
  if (!file) return;
  const content = await file.text();
  els.fasta.value = content;
  setStatus(`Loaded ${file.name}.`);
});

els.example.addEventListener('click', () => {
  els.fasta.value = EXAMPLE;
  setStatus('Example loaded. Click Analyze.');
});

els.clear.addEventListener('click', () => {
  els.fasta.value = '';
  els.results.innerHTML = '';
  els.results.hidden = true;
  els.downloads.hidden = true;
  lastResults = null;
  setStatus('');
});

els.analyze.addEventListener('click', run);

els.dlTsv.addEventListener('click', () => {
  if (lastResults) download('antibody_numbering.tsv', toTsv(lastResults), 'text/tab-separated-values;charset=utf-8');
});
els.dlJson.addEventListener('click', () => {
  if (lastResults) download('antibody_numbering.json', toJson(lastResults), 'application/json;charset=utf-8');
});
els.dlAligned.addEventListener('click', () => {
  if (lastResults) download('antibody_numbering_numbered.txt', toAligned(lastResults), 'text/plain;charset=utf-8');
});

// Service worker + offline badge.
if ('serviceWorker' in navigator) {
  navigator.serviceWorker.register('sw.js').then(() => {
    if (navigator.serviceWorker.controller || navigator.onLine) {
      els.offlineBadge.hidden = false;
    }
  }).catch(() => {});
}

window.addEventListener('online', () => { els.offlineBadge.hidden = false; });
