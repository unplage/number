// Port of ANARCI's Kabat heavy-chain conversion (GPLv3).
// Source: oxpig/ANARCI lib/python/anarci/schemes.py
// number_kabat_heavy / _number_regions / smooth_insertions / get_cdr3_annotations / gap_missing
// Operates on a reconstructed 128-state IMGT HMM state vector.

const ALPHABET = [
  'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
  'AA','BB','CC','DD','EE','FF','GG','HH','II','JJ','KK','LL','MM','NN','OO','PP','QQ','RR','SS','TT','UU','VV','WW','XX','YY','ZZ',' '
];

// Python allows alphabet[-1] == ' '; mirror that for insertion index -1.
function alphaAt(i) {
  if (i < 0) return ' ';
  return ALPHABET[i];
}

const HEAVY_STATE_STRING =
  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX';
const HEAVY_REGION_STRING =
  '11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888';
const HEAVY_REGION_INDEX = { '1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7 };
const HEAVY_RELS = { 0: 0, 1: -1, 2: -1, 3: -5, 4: -5, 5: -8, 6: -12, 7: -15 };
const HEAVY_N_REGIONS = 8;
const HEAVY_EXCLUDE_DELETIONS = [2, 4, 6];

// stateVector entry: [[{ id: number, type: 'm'|'i'|'d' }], seqIndex: number|null]

function smoothInsertions(stateVector) {
  const enforcedPatterns = [
    [[25, 'm'], [26, 'm'], [27, 'm'], [28, 'i']],
    [[38, 'i'], [38, 'm'], [39, 'm'], [40, 'm']],
    [[54, 'm'], [55, 'm'], [56, 'm'], [57, 'i']],
    [[65, 'i'], [65, 'm'], [66, 'm'], [67, 'm']],
    [[103, 'm'], [104, 'm'], [105, 'm'], [106, 'i']],
    [[117, 'i'], [117, 'm'], [118, 'm'], [119, 'm']],
  ];

  let buffer = [];
  let reg = null;
  const out = [];

  const flush = () => {
    if (buffer.length === 0) return;
    const nins = buffer.filter((s) => s[0][1] === 'i').length;
    if (nins > 0) {
      if (reg === -1) {
        let ntDels = buffer[0][0][0] - 1;
        for (const [st] of buffer) {
          if (st[1] === 'd') ntDels += 1;
          else break;
        }
        if (ntDels >= nins) {
          const matchStates = buffer.filter(([st]) => st[1] === 'm').map(([st]) => st);
          const first = matchStates[0][0];
          const cleaned = buffer.filter(([st]) => st[1] !== 'd');
          const add = cleaned.length - matchStates.length;
          const extended = [];
          for (let i = first - add; i < first; i++) extended.push([i, 'm']);
          extended.push(...matchStates);
          for (let i = 0; i < cleaned.length; i++) {
            out.push([extended[i], cleaned[i][1]]);
          }
        } else {
          out.push(...buffer);
        }
      } else {
        const cleaned = buffer.filter(([st]) => st[1] !== 'd');
        let newStates;
        if (reg % 2) {
          const headCount = Math.max(0, cleaned.length - 3);
          newStates = [];
          for (let i = 0; i < headCount; i++) newStates.push(enforcedPatterns[reg][0]);
          newStates.push(...enforcedPatterns[reg].slice(Math.max(4 - cleaned.length, 1)));
        } else {
          newStates = enforcedPatterns[reg].slice(0, 3);
          const tailCount = Math.max(0, cleaned.length - 3);
          for (let i = 0; i < tailCount; i++) newStates.push(enforcedPatterns[reg][2]);
        }
        for (let i = 0; i < cleaned.length; i++) {
          out.push([newStates[i], cleaned[i][1]]);
        }
      }
    } else {
      out.push(...buffer);
    }
    buffer = [];
  };

  for (const [[stateId, stateType], si] of stateVector) {
    if (stateId < 23) {
      buffer.push([[stateId, stateType], si]);
      reg = -1;
    } else if (stateId >= 25 && stateId < 28) {
      buffer.push([[stateId, stateType], si]);
      reg = 0;
    } else if (stateId > 37 && stateId <= 40) {
      buffer.push([[stateId, stateType], si]);
      reg = 1;
    } else if (stateId >= 54 && stateId < 57) {
      buffer.push([[stateId, stateType], si]);
      reg = 2;
    } else if (stateId > 64 && stateId <= 67) {
      buffer.push([[stateId, stateType], si]);
      reg = 3;
    } else if (stateId >= 103 && stateId < 106) {
      buffer.push([[stateId, stateType], si]);
      reg = 4;
    } else if (stateId > 116 && stateId <= 119) {
      buffer.push([[stateId, stateType], si]);
      reg = 5;
    } else if (buffer.length !== 0) {
      flush();
      out.push([[stateId, stateType], si]);
    } else {
      out.push([[stateId, stateType], si]);
    }
  }
  flush();
  return out;
}

function numberRegions(sequence, stateVector, stateString, regionString, regionIndexDict, relsIn, nRegions, excludeDeletions) {
  const sv = smoothInsertions(stateVector);
  const rels = { ...relsIn };
  const regions = Array.from({ length: nRegions }, () => []);
  let insertion = -1;
  let previousStateId = 1;
  let previousStateType = 'd';
  let startIndex = null;
  let endIndex = null;
  let region = null;

  for (const [[stateId, stateType], si] of sv) {
    if (stateType !== 'i' || region === null) {
      region = regionIndexDict[regionString[stateId - 1]];
    }

    if (stateType === 'm') {
      if (stateString[stateId - 1] === 'I') {
        if (previousStateType !== 'd') insertion += 1;
        rels[region] -= 1;
      } else {
        insertion = -1;
      }
      regions[region].push([{ n: stateId + rels[region], ins: alphaAt(insertion) }, sequence[si]]);
      previousStateId = stateId;
      if (startIndex === null) startIndex = si;
      endIndex = si;
      previousStateType = stateType;
    } else if (stateType === 'i') {
      insertion += 1;
      regions[region].push([{ n: previousStateId + rels[region], ins: alphaAt(insertion) }, sequence[si]]);
      if (startIndex === null) startIndex = si;
      endIndex = si;
      previousStateType = stateType;
    } else {
      previousStateType = stateType;
      if (stateString[stateId - 1] === 'I') {
        rels[region] -= 1;
        continue;
      }
      insertion = -1;
      previousStateId = stateId;
    }

    if (insertion >= 25 && excludeDeletions.includes(region)) {
      insertion = 0;
    }
  }

  return { regions, startIndex, endIndex };
}

function getCdr3AnnotationsHeavy(length) {
  const az = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  const insertions = Math.max(length - 10, 0);
  if (insertions >= 27) throw new Error('Too many CDR3 insertions for Kabat heavy');
  const orderedDeletions = [100, 99, 98, 97, 96, 95, 101, 102, 94, 93];
  const startDel = Math.max(0, 10 - length);
  const items = [];
  for (let i = startDel; i < orderedDeletions.length; i++) items.push({ n: orderedDeletions[i], ins: ' ' });
  for (let i = 0; i < insertions; i++) items.push({ n: 100, ins: az[i] });
  items.sort((a, b) => {
    if (a.n !== b.n) return a.n - b.n;
    if (a.ins === ' ') return -1;
    if (b.ins === ' ') return 1;
    return a.ins < b.ins ? -1 : 1;
  });
  return items;
}

function gapMissing(numbering) {
  const num = [{ n: 0, ins: ' ', aa: null }];
  for (const [p, a] of numbering.flat()) {
    if (!p) continue; // skip slots that were never filled
    if (p.n > num[num.length - 1].n + 1) {
      for (let i = num[num.length - 1].n + 1; i < p.n; i++) {
        num.push({ n: i, ins: ' ', aa: '-' });
      }
    }
    num.push({ n: p.n, ins: p.ins, aa: a });
  }
  return num.slice(1);
}

/**
 * Apply the Kabat numbering scheme to a heavy chain given an IMGT state vector.
 */
export function numberKabatHeavy(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence,
    stateVector,
    HEAVY_STATE_STRING,
    HEAVY_REGION_STRING,
    HEAVY_REGION_INDEX,
    HEAVY_RELS,
    HEAVY_N_REGIONS,
    HEAVY_EXCLUDE_DELETIONS,
  );

  const numbering = [[], regions[1], [], regions[3], [], regions[5], [], regions[7]];

  // FR1 zone (region index 0): insertions at Kabat 6.
  {
    const region = regions[0];
    const insertions = region.filter((r) => r[0].ins !== ' ').length;
    if (insertions) {
      const start = region[0][0].n;
      const annotations = [];
      for (let i = start; i < 7; i++) annotations.push({ n: i, ins: ' ' });
      for (let i = 0; i < insertions; i++) annotations.push({ n: 6, ins: alphaAt(i) });
      annotations.push({ n: 7, ins: ' ' }, { n: 8, ins: ' ' }, { n: 9, ins: ' ' });
      numbering[0] = region.map(([, aa], i) => [annotations[i], aa]);
    } else {
      numbering[0] = region;
    }
  }

  // CDR1 zone (region index 2): insertions onto 35.
  {
    const region = regions[2];
    const length = region.length;
    const insertions = Math.max(0, length - 13);
    const annotations = [];
    for (let i = 23; i < 36 && annotations.length < length; i++) annotations.push({ n: i, ins: ' ' });
    for (let i = 0; i < insertions; i++) annotations.push({ n: 35, ins: alphaAt(i) });
    numbering[2] = region.map(([, aa], i) => [annotations[i], aa]);
  }

  // CDR2 zone (region index 4): insertions onto 52.
  {
    const region = regions[4];
    const length = region.length;
    const insertions = Math.max(length - 8, 0);
    const annotations = [];
    const head = [{ n: 50, ins: ' ' }, { n: 51, ins: ' ' }, { n: 52, ins: ' ' }];
    annotations.push(...head.slice(0, Math.max(0, length - 5)));
    for (let i = 0; i < insertions; i++) annotations.push({ n: 52, ins: alphaAt(i) });
    const tail = [{ n: 53, ins: ' ' }, { n: 54, ins: ' ' }, { n: 55, ins: ' ' }, { n: 56, ins: ' ' }, { n: 57, ins: ' ' }];
    annotations.push(...tail.slice(Math.abs(Math.min(0, length - 5))));
    numbering[4] = region.map(([, aa], i) => [annotations[i], aa]);
  }

  // CDR3 zone (region index 6): renumber by length.
  {
    const region = regions[6];
    const length = region.length;
    if (length > 36) return { numbering: [], startIndex, endIndex };
    const annotations = getCdr3AnnotationsHeavy(length);
    numbering[6] = region.map(([, aa], i) => [annotations[i], aa]);
  }

  return { numbering: gapMissing(numbering), startIndex, endIndex };
}

/**
 * Reconstruct an IMGT HMM state vector from an anarcism domain result.
 */
export function domainToStateVector(domain) {
  const sv = [];
  const insertionsByPos = new Map();
  const cdr3InsertionResidues = [];
  for (const r of domain.numbering) {
    if (!r.insertionCode) continue;
    if (r.position >= 105 && r.position <= 117) {
      cdr3InsertionResidues.push(r);
    } else {
      if (!insertionsByPos.has(r.position)) insertionsByPos.set(r.position, []);
      insertionsByPos.get(r.position).push(r);
    }
  }

  const matchByPos = new Map();
  for (const r of domain.numbering) {
    if (!r.insertionCode) matchByPos.set(r.position, r);
  }

  const seqParts = [];
  let cursor = 0;

  for (let s = 1; s <= 128; s++) {
    const ch = domain.paddedImgtAlignment[s - 1];
    if (ch === '-') {
      sv.push([[s, 'd'], null]);
      continue;
    }
    const match = matchByPos.get(s);
    const aa = match ? match.aminoAcid : ch;
    seqParts.push(aa);
    sv.push([[s, 'm'], cursor]);
    cursor += 1;
    const ins = insertionsByPos.get(s);
    if (ins) {
      for (const r of ins) {
        sv.push([[s, 'i'], cursor]);
        seqParts.push(r.aminoAcid);
        cursor += 1;
      }
    }
    // CDR3 insertions (IMGT labels 111x/112x) are placed between match state
    // 111 and 112, preserving sequence order. Kabat CDR3 renumbering is purely
    // length-based, so the precise attach point within the loop does not affect
    // the resulting labels; placing them here keeps region assignment correct.
    if (s === 111) {
      for (const r of cdr3InsertionResidues) {
        sv.push([[111, 'i'], cursor]);
        seqParts.push(r.aminoAcid);
        cursor += 1;
      }
    }
  }

  return { stateVector: sv, sequence: seqParts.join('') };
}

/**
 * Convert an anarcism domain (IMGT) to Kabat heavy numbering.
 * Returns Array<{ pos: string, aa: string }> e.g. { pos: 'H100A', aa: 'Y' }.
 */
export function kabatHeavyFromDomain(domain) {
  const { stateVector, sequence } = domainToStateVector(domain);
  const { numbering } = numberKabatHeavy(stateVector, sequence);
  return numbering
    .filter((r) => r.aa && r.aa !== '-')
    .map((r) => ({ pos: 'H' + r.n + (r.ins && r.ins !== ' ' ? r.ins : ''), aa: r.aa }));
}
