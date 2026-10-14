// Port of ANARCI numbering schemes (BSD-3-Clause).
// Source: oxpig/ANARCI lib/python/anarci/schemes.py
//
// All schemes convert from the 128-state IMGT HMM state vector, exactly as ANARCI
// does. The state vector is reconstructed from anarcism's IMGT output by
// `domainToStateVector` (see below).
//
// Chain letters in the output positions: H = heavy, L = light (kappa/lambda).

import blosum62 from './data/blosum62.js';

const ALPHABET = [
  'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
  'AA','BB','CC','DD','EE','FF','GG','HH','II','JJ','KK','LL','MM','NN','OO','PP','QQ','RR','SS','TT','UU','VV','WW','XX','YY','ZZ',' ',
];

function alphaAt(i) {
  if (i < 0) return ' ';
  return ALPHABET[i];
}

// ---------------------------------------------------------------------------
// Generic alignment -> numbering machinery
// ---------------------------------------------------------------------------

// stateVector entries: [[stateId, stateType], seqIndex|null]
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
          for (let i = 0; i < cleaned.length; i++) out.push([extended[i], cleaned[i][1]]);
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
        for (let i = 0; i < cleaned.length; i++) out.push([newStates[i], cleaned[i][1]]);
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
    if (stateType !== 'i' || region === null) region = regionIndexDict[regionString[stateId - 1]];

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

    if (insertion >= 25 && excludeDeletions.includes(region)) insertion = 0;
  }

  return { regions, startIndex, endIndex };
}

function gapMissing(numbering) {
  const num = [{ n: 0, ins: ' ', aa: null }];
  for (const [p, a] of numbering.flat()) {
    if (!p) continue;
    if (p.n > num[num.length - 1].n + 1) {
      for (let i = num[num.length - 1].n + 1; i < p.n; i++) num.push({ n: i, ins: ' ', aa: '-' });
    }
    num.push({ n: p.n, ins: p.ins, aa: a });
  }
  return num.slice(1);
}

function getCdr3Annotations(length, chainType) {
  const az = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  if (chainType === 'heavy') {
    const insertions = Math.max(length - 10, 0);
    if (insertions >= 27) throw new Error('Too many CDR3 insertions');
    const ordered = [100, 99, 98, 97, 96, 95, 101, 102, 94, 93];
    const items = [];
    for (let i = Math.max(0, 10 - length); i < ordered.length; i++) items.push({ n: ordered[i], ins: ' ' });
    for (let i = 0; i < insertions; i++) items.push({ n: 100, ins: az[i] });
    items.sort((a, b) => (a.n - b.n) || (a.ins === ' ' ? -1 : b.ins === ' ' ? 1 : a.ins < b.ins ? -1 : 1));
    return items;
  }
  const insertions = Math.max(length - 9, 0);
  if (insertions >= 27) throw new Error('Too many CDR3 insertions');
  const ordered = [95, 94, 93, 92, 91, 96, 97, 90, 89];
  const items = [];
  for (let i = Math.max(0, 9 - length); i < ordered.length; i++) items.push({ n: ordered[i], ins: ' ' });
  for (let i = 0; i < insertions; i++) items.push({ n: 95, ins: az[i] });
  items.sort((a, b) => (a.n - b.n) || (a.ins === ' ' ? -1 : b.ins === ' ' ? 1 : a.ins < b.ins ? -1 : 1));
  return items;
}

// ---------------------------------------------------------------------------
// Kabat
// ---------------------------------------------------------------------------

const KABAT_H = {
  state: 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX',
  region: '11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888',
  index: { '1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7 },
  rels: { 0: 0, 1: -1, 2: -1, 3: -5, 4: -5, 5: -8, 6: -12, 7: -15 },
  n: 8,
  excl: [2, 4, 6],
};

function numberKabatHeavy(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, KABAT_H.state, KABAT_H.region, KABAT_H.index, KABAT_H.rels, KABAT_H.n, KABAT_H.excl,
  );
  const numbering = [[], regions[1], [], regions[3], [], regions[5], [], regions[7]];
  // FR1 zone
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
  // CDR1 -> 35
  {
    const region = regions[2];
    const length = region.length;
    const insertions = Math.max(0, length - 13);
    const annotations = [];
    for (let i = 23; i < 36 && annotations.length < length; i++) annotations.push({ n: i, ins: ' ' });
    for (let i = 0; i < insertions; i++) annotations.push({ n: 35, ins: alphaAt(i) });
    numbering[2] = region.map(([, aa], i) => [annotations[i], aa]);
  }
  // CDR2 -> 52
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
  // CDR3
  {
    const region = regions[6];
    const length = region.length;
    if (length > 36) return { numbering: [], startIndex, endIndex };
    const annotations = getCdr3Annotations(length, 'heavy');
    numbering[6] = region.map(([, aa], i) => [annotations[i], aa]);
  }
  return { numbering: gapMissing(numbering), startIndex, endIndex };
}

// ---------------------------------------------------------------------------
// Chothia (Martin light reuses Chothia light)
// ---------------------------------------------------------------------------

const CHOTHIA_H = {
  state: 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX',
  region: '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888',
  index: { '1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7 },
  rels: { 0: 0, 1: -1, 2: -1, 3: -5, 4: -5, 5: -8, 6: -12, 7: -15 },
  n: 8,
  excl: [0, 2, 4, 6],
};

function numberChothiaHeavy(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, CHOTHIA_H.state, CHOTHIA_H.region, CHOTHIA_H.index, CHOTHIA_H.rels, CHOTHIA_H.n, CHOTHIA_H.excl,
  );
  const numbering = [[], regions[1], [], regions[3], [], regions[5], [], regions[7]];
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
  {
    const region = regions[2];
    const length = region.length;
    const insertions = Math.max(length - 11, 0);
    let annotations;
    if (insertions) {
      annotations = [];
      for (let i = 23; i < 32; i++) annotations.push({ n: i, ins: ' ' });
      for (let i = 0; i < insertions; i++) annotations.push({ n: 31, ins: alphaAt(i) });
      annotations.push({ n: 32, ins: ' ' }, { n: 33, ins: ' ' });
    } else {
      annotations = [];
      for (let i = 23; i < 32; i++) annotations.push({ n: i, ins: ' ' });
      annotations = annotations.slice(0, Math.max(0, length - 2)).concat([{ n: 32, ins: ' ' }, { n: 33, ins: ' ' }].slice(0, length));
    }
    numbering[2] = region.map(([, aa], i) => [annotations[i], aa]);
  }
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
  {
    const region = regions[6];
    const length = region.length;
    if (length > 36) return { numbering: [], startIndex, endIndex };
    const annotations = getCdr3Annotations(length, 'heavy');
    numbering[6] = region.map(([, aa], i) => [annotations[i], aa]);
  }
  return { numbering: gapMissing(numbering), startIndex, endIndex };
}

const CHOTHIA_L = {
  state: 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX',
  region: '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777',
  index: { '1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6 },
  rels: { 0: 0, 1: 0, 2: -6, 3: -6, 4: -13, 5: -16, 6: -20 },
  n: 7,
  excl: [1, 3, 4, 5],
};

function numberKabatLight(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, CHOTHIA_L.state, CHOTHIA_L.region, CHOTHIA_L.index, CHOTHIA_L.rels, 7, [1, 3, 5],
  );
  const numbering = [regions[0], [], regions[2], [], regions[4], [], regions[6]];
  // CDR1 -> 27
  {
    const region = regions[1];
    const length = region.length;
    const insertions = Math.max(length - 11, 0);
    const base = [];
    for (let i = 24; i <= 27; i++) base.push({ n: i, ins: ' ' });
    const ann = base.slice(0, Math.max(0, length));
    for (let i = 0; i < insertions; i++) ann.push({ n: 27, ins: alphaAt(i) });
    const tail = [];
    for (let i = 28; i <= 34; i++) tail.push({ n: i, ins: ' ' });
    ann.push(...tail.slice(Math.abs(Math.min(0, length - 11))));
    numbering[1] = region.map(([, aa], i) => [ann[i], aa]);
  }
  // CDR2 -> 52
  {
    const region = regions[3];
    const length = region.length;
    const insertions = Math.max(length - 4, 0);
    if (insertions > 0) {
      const ann = [{ n: 51, ins: ' ' }, { n: 52, ins: ' ' }];
      for (let i = 0; i < insertions; i++) ann.push({ n: 52, ins: alphaAt(i) });
      ann.push({ n: 53, ins: ' ' }, { n: 54, ins: ' ' });
      numbering[3] = region.map(([, aa], i) => [ann[i], aa]);
    } else {
      numbering[3] = region;
    }
  }
  // CDR3
  {
    const region = regions[5];
    const length = region.length;
    if (length > 35) return { numbering: [], startIndex, endIndex };
    const ann = getCdr3Annotations(length, 'light');
    numbering[5] = region.map(([, aa], i) => [ann[i], aa]);
  }
  return { numbering: gapMissing(numbering), startIndex, endIndex };
}

function numberChothiaLight(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, CHOTHIA_L.state, CHOTHIA_L.region, CHOTHIA_L.index, CHOTHIA_L.rels, CHOTHIA_L.n, CHOTHIA_L.excl,
  );
  const numbering = [regions[0], [], regions[2], [], regions[4], [], regions[6]];
  // CDR1 -> 30
  {
    const region = regions[1];
    const length = region.length;
    const insertions = Math.max(length - 11, 0);
    const annotations = [];
    for (let i = 24; i <= 30; i++) annotations.push({ n: i, ins: ' ' });
    const capped = annotations.slice(0, Math.max(0, length));
    for (let i = 0; i < insertions; i++) capped.push({ n: 30, ins: alphaAt(i) });
    const tail = [{ n: 31, ins: ' ' }, { n: 32, ins: ' ' }, { n: 33, ins: ' ' }, { n: 34, ins: ' ' }];
    capped.push(...tail.slice(Math.abs(Math.min(0, length - 11))));
    numbering[1] = region.map(([, aa], i) => [capped[i], aa]);
  }
  // CDR2 -> 52
  {
    const region = regions[3];
    const length = region.length;
    const insertions = Math.max(length - 4, 0);
    if (insertions > 0) {
      const annotations = [{ n: 51, ins: ' ' }, { n: 52, ins: ' ' }];
      for (let i = 0; i < insertions; i++) annotations.push({ n: 52, ins: alphaAt(i) });
      annotations.push({ n: 53, ins: ' ' }, { n: 54, ins: ' ' });
      numbering[3] = region.map(([, aa], i) => [annotations[i], aa]);
    } else {
      numbering[3] = region;
    }
  }
  // FW3 -> 68
  {
    const region = regions[4];
    const length = region.length;
    const insertions = Math.max(length - 34, 0);
    if (insertions > 0) {
      const annotations = [];
      for (let i = 55; i < 69; i++) annotations.push({ n: i, ins: ' ' });
      for (let i = 0; i < insertions; i++) annotations.push({ n: 68, ins: alphaAt(i) });
      for (let i = 69; i < 89; i++) annotations.push({ n: i, ins: ' ' });
      numbering[4] = region.map(([, aa], i) => [annotations[i], aa]);
    } else if (length === 33) {
      const annotations = [];
      for (let i = 55; i < 68; i++) annotations.push({ n: i, ins: ' ' });
      for (let i = 69; i < 89; i++) annotations.push({ n: i, ins: ' ' });
      numbering[4] = region.map(([, aa], i) => [annotations[i], aa]);
    } else {
      numbering[4] = region;
    }
  }
  // CDR3
  {
    const region = regions[5];
    const length = region.length;
    if (length > 35) return { numbering: [], startIndex, endIndex };
    const annotations = getCdr3Annotations(length, 'light');
    numbering[5] = region.map(([, aa], i) => [annotations[i], aa]);
  }
  return { numbering: gapMissing(numbering), startIndex, endIndex };
}

// ---------------------------------------------------------------------------
// Martin (extended Chothia) — light == chothia light
// ---------------------------------------------------------------------------

const MARTIN_H = {
  state: 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX',
  region: '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888',
  index: { '1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6, '8': 7 },
  rels: { 0: 0, 1: -1, 2: -1, 3: -5, 4: -5, 5: -8, 6: -12, 7: -15 },
  n: 8,
  excl: [2, 4, 5, 6],
};

function numberMartinHeavy(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, MARTIN_H.state, MARTIN_H.region, MARTIN_H.index, MARTIN_H.rels, MARTIN_H.n, MARTIN_H.excl,
  );
  const numbering = [[], regions[1], [], regions[3], [], regions[5], [], regions[7]];
  {
    const region = regions[0];
    const insertions = region.filter((r) => r[0].ins !== ' ').length;
    if (insertions) {
      const start = region[0][0].n;
      const annotations = [];
      for (let i = start; i < 9; i++) annotations.push({ n: i, ins: ' ' });
      for (let i = 0; i < insertions; i++) annotations.push({ n: 8, ins: alphaAt(i) });
      annotations.push({ n: 9, ins: ' ' });
      numbering[0] = region.map(([, aa], i) => [annotations[i], aa]);
    } else {
      numbering[0] = region;
    }
  }
  {
    const region = regions[2];
    const length = region.length;
    const insertions = Math.max(length - 11, 0);
    let annotations;
    if (insertions) {
      annotations = [];
      for (let i = 23; i < 32; i++) annotations.push({ n: i, ins: ' ' });
      for (let i = 0; i < insertions; i++) annotations.push({ n: 31, ins: alphaAt(i) });
      annotations.push({ n: 32, ins: ' ' }, { n: 33, ins: ' ' });
    } else {
      annotations = [];
      for (let i = 23; i < 32; i++) annotations.push({ n: i, ins: ' ' });
      annotations = annotations.slice(0, Math.max(0, length - 2)).concat([{ n: 32, ins: ' ' }, { n: 33, ins: ' ' }].slice(0, length));
    }
    numbering[2] = region.map(([, aa], i) => [annotations[i], aa]);
  }
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
  // FW3 -> 72
  {
    const region = regions[5];
    const length = region.length;
    const insertions = Math.max(length - 35, 0);
    if (insertions > 0) {
      const annotations = [];
      for (let i = 58; i < 73; i++) annotations.push({ n: i, ins: ' ' });
      for (let i = 0; i < insertions; i++) annotations.push({ n: 72, ins: alphaAt(i) });
      for (let i = 73; i < 93; i++) annotations.push({ n: i, ins: ' ' });
      numbering[5] = region.map(([, aa], i) => [annotations[i], aa]);
    } else {
      numbering[5] = region;
    }
  }
  {
    const region = regions[6];
    const length = region.length;
    if (length > 36) return { numbering: [], startIndex, endIndex };
    const annotations = getCdr3Annotations(length, 'heavy');
    numbering[6] = region.map(([, aa], i) => [annotations[i], aa]);
  }
  return { numbering: gapMissing(numbering), startIndex, endIndex };
}

// ---------------------------------------------------------------------------
// AHo
// ---------------------------------------------------------------------------

const AHO_CONFIG = {
  state: 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
  region: 'BBBBBBBBBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFHHHHHHHHHHHHHHHHIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK',
  index: { A: 0, B: 1, C: 2, D: 3, E: 4, F: 5, G: 6, H: 7, I: 8, J: 9, K: 10 },
  rels: { 0: 0, 1: 0, 2: 0, 3: 0, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2, 10: 21 },
  n: 11,
  excl: [1, 3, 4, 5, 7, 9],
};

const AHO_CDR1 = {
  L: [28, 36, 35, 37, 34, 38, 27, 29, 33, 39, 32, 40, 26, 30, 25, 31, 41, 42],
  K: [28, 27, 36, 35, 37, 34, 38, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  H: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  A: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  B: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  D: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
  G: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
};

function numberAho(stateVector, sequence, chainType) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, AHO_CONFIG.state, AHO_CONFIG.region, AHO_CONFIG.index, AHO_CONFIG.rels, AHO_CONFIG.n, AHO_CONFIG.excl,
  );
  const numbering = [regions[0], regions[1], regions[2], [], regions[4], [], regions[6], [], regions[8], regions[9], regions[10]];

  // Move indel in FW1 onto 8
  {
    const region = regions[1];
    const length = region.length;
    if (length > 0) {
      const start = region[0][0].n;
      const stretchLen = 10 - (start - 1);
      let annotations;
      if (length > stretchLen) {
        annotations = [];
        for (let i = start; i < 9; i++) annotations.push({ n: i, ins: ' ' });
        for (let i = 0; i < length - stretchLen; i++) annotations.push({ n: 8, ins: alphaAt(i) });
        annotations.push({ n: 9, ins: ' ' }, { n: 10, ins: ' ' });
      } else {
        const od = [{ n: 8, ins: ' ' }];
        for (let i = start; i <= 10; i++) if (i !== 8) od.push({ n: i, ins: ' ' });
        annotations = od.slice(Math.max(stretchLen - length, 0));
        annotations.sort((a, b) => a.n - b.n);
      }
      numbering[1] = region.map(([, aa], i) => [annotations[i], aa]);
    }
  }

  // CDR1 (region D, index 3)
  {
    const region = regions[3];
    const length = region.length;
    const ordered = AHO_CDR1[chainType] || AHO_CDR1.L;
    const annotations = [];
    for (let i = Math.max(18 - length, 0); i < ordered.length; i++) annotations.push(ordered[i]);
    annotations.sort((a, b) => a - b);
    const annObj = annotations.map((n) => ({ n, ins: ' ' }));
    const insertions = Math.max(length - 18, 0);
    if (insertions > 26) return { numbering: [], startIndex, endIndex };
    if (insertions > 0) {
      const insertat = annObj.findIndex((a) => a.n === 36) + 1;
      const insArr = [];
      for (let i = 0; i < insertions; i++) insArr.push({ n: 36, ins: alphaAt(i) });
      annObj.splice(insertat, 0, ...insArr);
    }
    numbering[3] = region.map(([, aa], i) => [annObj[i], aa]);
  }

  // CDR2 (region F, index 5)
  {
    const region = regions[5];
    const length = region.length;
    const ordered = chainType === 'A'
      ? [74, 73, 63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 75, 76, 77]
      : [63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77];
    const annObj = ordered.slice(Math.max(20 - length, 0)).map((n) => ({ n, ins: ' ' })).sort((a, b) => a.n - b.n);
    const insertions = Math.max(length - 20, 0);
    if (insertions > 26) return { numbering: [], startIndex, endIndex };
    if (insertions > 0) {
      const insertat = annObj.findIndex((a) => a.n === 63) + 1;
      const insArr = [];
      for (let i = 0; i < insertions; i++) insArr.push({ n: 63, ins: alphaAt(i) });
      annObj.splice(insertat, 0, ...insArr);
    }
    numbering[5] = region.map(([, aa], i) => [annObj[i], aa]);
  }

  // FW3 (region H, index 7)
  {
    const region = regions[7];
    const length = region.length;
    const ordered = [86, 85, 87, 84, 88, 83, 89, 82, 90, 81, 91, 80, 92, 79, 93, 78];
    const annObj = ordered.slice(Math.max(16 - length, 0)).map((n) => ({ n, ins: ' ' })).sort((a, b) => a.n - b.n);
    const insertions = Math.max(length - 16, 0);
    if (insertions > 26) return { numbering: [], startIndex, endIndex };
    if (insertions > 0) {
      const insertat = annObj.findIndex((a) => a.n === 85) + 1;
      const insArr = [];
      for (let i = 0; i < insertions; i++) insArr.push({ n: 85, ins: alphaAt(i) });
      annObj.splice(insertat, 0, ...insArr);
    }
    numbering[7] = region.map(([, aa], i) => [annObj[i], aa]);
  }

  // CDR3 (region J, index 9)
  {
    const region = regions[9];
    const length = region.length;
    const ordered = [123, 124, 122, 125, 121, 126, 120, 127, 119, 128, 118, 129, 117, 130, 116, 131, 115, 132, 114, 133, 113, 134, 112, 135, 111, 136, 110, 137, 109, 138, 108, 107];
    const annObj = ordered.slice(Math.max(32 - length, 0)).map((n) => ({ n, ins: ' ' })).sort((a, b) => a.n - b.n);
    const insertions = Math.max(length - 32, 0);
    if (insertions > 26) return { numbering: [], startIndex, endIndex };
    if (insertions > 0) {
      const insertat = annObj.findIndex((a) => a.n === 123) + 1;
      const insArr = [];
      for (let i = 0; i < insertions; i++) insArr.push({ n: 123, ins: alphaAt(i) });
      annObj.splice(insertat, 0, ...insArr);
    }
    numbering[9] = region.map(([, aa], i) => [annObj[i], aa]);
  }

  let numberingOut = gapMissing(numbering);
  let end = endIndex;
  if (numberingOut.length > 0) {
    const last = numberingOut[numberingOut.length - 1];
    if (last.n === 148 && last.aa !== '-' && end + 1 < sequence.length) {
      numberingOut = numberingOut.concat([{ n: 149, ins: ' ', aa: sequence[end + 1] }]);
      end += 1;
    }
  }
  return { numbering: numberingOut, startIndex, endIndex: end };
}

// ---------------------------------------------------------------------------
// Wolfguy
// ---------------------------------------------------------------------------

const WOLFGUY_H = {
  state: 'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
  region: '11111111111111111111111111222222222222223333333333333344444444444444444444555555555555555555555555555555666666666666677777777777',
  index: { '1': 0, '2': 1, '3': 2, '4': 3, '5': 4, '6': 5, '7': 6 },
  rels: { 0: 100, 1: 124, 2: 160, 3: 196, 4: 226, 5: 244, 6: 283 },
  n: 7,
  excl: [1, 3, 5],
};

function numberWolfguyHeavy(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, WOLFGUY_H.state, WOLFGUY_H.region, WOLFGUY_H.index, WOLFGUY_H.rels, WOLFGUY_H.n, WOLFGUY_H.excl,
  );
  const numbering = [regions[0], [], regions[2], [], regions[4], [], regions[6]];
  {
    const orderedDeletions = [151];
    for (let p1 = 152, p2 = 199; p1 < 176; p1++, p2--) orderedDeletions.push(p1, p2);
    const length = regions[1].length;
    const annotations = orderedDeletions.slice(0, length).sort((a, b) => a - b);
    numbering[1] = regions[1].map(([, aa], i) => [{ n: annotations[i], ins: ' ' }, aa]);
  }
  {
    let orderedDeletions = [251];
    for (let p1 = 252, p2 = 290; p1 < 271; p1++, p2--) orderedDeletions.push(p1, p2);
    orderedDeletions.push(271);
    orderedDeletions = Array.from({ length: 9 }, (_, i) => 299 - i).concat(orderedDeletions);
    const length = regions[3].length;
    const annotations = orderedDeletions.slice(0, length).sort((a, b) => a - b);
    numbering[3] = regions[3].map(([, aa], i) => [{ n: annotations[i], ins: ' ' }, aa]);
  }
  {
    let orderedDeletions = [];
    for (let p1 = 356, p2 = 391; p1 < 374; p1++, p2--) orderedDeletions.push(p1, p2);
    orderedDeletions = [354, 394, 355, 393, 392].concat(orderedDeletions);
    orderedDeletions = [331, 332, 399, 398, 351, 352, 397, 353, 396, 395].concat(orderedDeletions);
    const length = regions[5].length;
    if (length > orderedDeletions.length) return { numbering: [], startIndex, endIndex };
    const annotations = orderedDeletions.slice(0, length).sort((a, b) => a - b);
    numbering[5] = regions[5].map(([, aa], i) => [{ n: annotations[i], ins: ' ' }, aa]);
  }
  return { numbering: numbering.flat(), startIndex, endIndex };
}

const WOLFGUY_L = {
  state: 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
  region: '1111111AAABBBBBBBBBBBBB222222222222222223333333333333334444444444444455555555555666677777777777777777777888888888888899999999999',
  index: { '1': 0, A: 1, B: 2, '2': 3, '3': 4, '4': 5, '5': 6, '6': 7, '7': 8, '8': 9, '9': 10 },
  rels: { 0: 500, 1: 500, 2: 500, 3: 527, 4: 560, 5: 595, 6: 631, 7: 630, 8: 630, 9: 646, 10: 683 },
  n: 11,
  excl: [1, 3, 5, 7, 9],
};

const WOLFGUY_L1 = {
  9: [['9', 'XXXXXXXXX', [551, 552, 554, 556, 563, 572, 597, 598, 599]]],
  10: [['10', 'XXXXXXXXXX', [551, 552, 553, 556, 561, 562, 571, 597, 598, 599]]],
  11: [
    ['11a', 'RASQDISSYLA', [551, 552, 553, 556, 561, 562, 571, 596, 597, 598, 599]],
    ['11b', 'GGNNIGSKSVH', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]],
    ['11b.2', 'SGDQLPKKYAY', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]],
  ],
  12: [
    ['12a', 'TLSSQHSTYTIE', [551, 552, 553, 554, 555, 556, 561, 563, 572, 597, 598, 599]],
    ['12b', 'TASSSVSSSYLH', [551, 552, 553, 556, 561, 562, 571, 595, 596, 597, 598, 599]],
    ['12c', 'RASQSVxNNYLA', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]],
    ['12d', 'rSShSIrSrrVh', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]],
  ],
  13: [
    ['13a', 'SGSSSNIGNNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 597, 598, 599]],
    ['13b', 'TRSSGSLANYYVQ', [551, 552, 553, 554, 556, 561, 562, 563, 571, 572, 597, 598, 599]],
  ],
  14: [
    ['14a', 'RSSTGAVTTSNYAN', [551, 552, 553, 554, 555, 561, 562, 563, 564, 571, 572, 597, 598, 599]],
    ['14b', 'TGTSSDVGGYNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 596, 597, 598, 599]],
  ],
  15: [['15', 'XXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 594, 595, 596, 597, 598, 599]]],
  16: [['16', 'XXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 594, 595, 596, 597, 598, 599]]],
  17: [['17', 'XXXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 584, 594, 595, 596, 597, 598, 599]]],
};

function blosum(a, b) {
  const key = a.toUpperCase() + b.toUpperCase();
  if (key in blosum62) return blosum62[key];
  const swapped = b.toUpperCase() + a.toUpperCase();
  if (swapped in blosum62) return blosum62[swapped];
  return 0;
}

function getWolfguyL1(seqRegion, length) {
  if (length in WOLFGUY_L1) {
    let best = null;
    let bestScore = -1e9;
    for (const canonical of WOLFGUY_L1[length]) {
      let sub = 0;
      for (let i = 0; i < length; i++) sub += blosum(seqRegion[i][1], canonical[1][i]);
      if (sub > bestScore) {
        bestScore = sub;
        best = canonical;
      }
    }
    return best[2];
  }
  const ordered = [];
  for (let p1 = 551, p2 = 599; p1 < 575; p1++, p2--) ordered.push(p2, p1);
  ordered.push(575);
  return ordered.slice(0, length).sort((a, b) => a - b);
}

function numberWolfguyLight(stateVector, sequence) {
  const { regions, startIndex, endIndex } = numberRegions(
    sequence, stateVector, WOLFGUY_L.state, WOLFGUY_L.region, WOLFGUY_L.index, WOLFGUY_L.rels, WOLFGUY_L.n, WOLFGUY_L.excl,
  );
  const numbering = [regions[0], [], regions[2], [], regions[4], [], regions[6], [], regions[8], [], regions[10]];
  {
    const length = regions[1].length;
    const head = [{ n: 510, ins: ' ' }, { n: 509, ins: ' ' }, { n: 508, ins: ' ' }].slice(0, length);
    const insArr = [];
    for (let i = 0; i < Math.max(0, length - 3); i++) insArr.push({ n: 508, ins: alphaAt(i) });
    const annotations = head.concat(insArr).sort((a, b) => a.n - b.n || (a.ins === ' ' ? -1 : b.ins === ' ' ? 1 : a.ins < b.ins ? -1 : 1));
    numbering[1] = regions[1].map(([, aa], i) => [annotations[i], aa]);
  }
  {
    const length = regions[3].length;
    const annotations = getWolfguyL1(regions[3], length);
    numbering[3] = regions[3].map(([, aa], i) => [{ n: annotations[i], ins: ' ' }, aa]);
  }
  {
    let ordered = [];
    for (let p1 = 652, p2 = 694; p1 < 673; p1++, p2--) ordered.push(p2, p1);
    ordered = [651].concat(Array.from({ length: 5 }, (_, i) => 699 - i), ordered, [673]);
    const length = regions[5].length;
    const annotations = ordered.slice(0, length).sort((a, b) => a - b);
    numbering[5] = regions[5].map(([, aa], i) => [{ n: annotations[i], ins: ' ' }, aa]);
  }
  {
    const length = regions[7].length;
    const insertions = Math.max(0, length - 4);
    const head = [{ n: 711, ins: ' ' }, { n: 712, ins: ' ' }, { n: 713, ins: ' ' }, { n: 714, ins: ' ' }].slice(0, length);
    const insArr = [];
    for (let i = 0; i < insertions; i++) insArr.push({ n: 714, ins: alphaAt(i) });
    numbering[7] = regions[7].map(([, aa], i) => [head.concat(insArr)[i], aa]);
  }
  {
    const ordered = [];
    for (let p1 = 751, p2 = 799; p1 < 775; p1++, p2--) ordered.push(p1, p2);
    ordered.push(775);
    const length = regions[9].length;
    if (length > ordered.length) return { numbering: [], startIndex, endIndex };
    const annotations = ordered.slice(0, length).sort((a, b) => a - b);
    numbering[9] = regions[9].map(([, aa], i) => [{ n: annotations[i], ins: ' ' }, aa]);
  }
  return { numbering: numbering.flat(), startIndex, endIndex };
}

// ---------------------------------------------------------------------------
// Dispatch + formatting
// ---------------------------------------------------------------------------

/**
 * Convert a 128-state IMGT state vector to a numbering scheme.
 * @returns {Array<{n:number, ins:string, aa:string}>} residues in sequence order
 */
export function convertStateVector(scheme, chainLetter, stateVector, sequence) {
  const ct = chainLetter === 'H' ? 'H' : chainLetter;
  let res;
  switch (scheme) {
    case 'kabat':
      res = ct === 'H' ? numberKabatHeavy(stateVector, sequence) : numberKabatLight(stateVector, sequence);
      break;
    case 'chothia':
      res = ct === 'H' ? numberChothiaHeavy(stateVector, sequence) : numberChothiaLight(stateVector, sequence);
      break;
    case 'martin':
      res = ct === 'H' ? numberMartinHeavy(stateVector, sequence) : numberChothiaLight(stateVector, sequence);
      break;
    case 'aho':
      res = numberAho(stateVector, sequence, ct);
      break;
    case 'wolfguy':
      res = ct === 'H' ? numberWolfguyHeavy(stateVector, sequence) : numberWolfguyLight(stateVector, sequence);
      break;
    default:
      throw new Error(`Unknown scheme: ${scheme}`);
  }
  return res.numbering
    .map((r) => (Array.isArray(r) ? { n: r[0].n, ins: r[0].ins, aa: r[1] } : r))
    .filter((r) => r.aa && r.aa !== '-');
}

// ---------------------------------------------------------------------------
// IMGT state vector reconstruction from anarcism domain output
// ---------------------------------------------------------------------------

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
  for (const r of domain.numbering) if (!r.insertionCode) matchByPos.set(r.position, r);

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
    // CDR3 insertions (IMGT 111x/112x) go between match states 111 and 112.
    // Kabat/Chothia/Martin CDR3 renumbering is length-based, so the exact attach
    // point does not affect labels, but it keeps them inside region 105-117.
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
 * Number an anarcism domain with a scheme (except IMGT, which uses anarcism directly).
 * @returns {Array<{n:number, ins:string, aa:string}>}
 */
export function schemeNumberingFromDomain(scheme, chainLetter, domain) {
  const { stateVector, sequence } = domainToStateVector(domain);
  return convertStateVector(scheme, chainLetter, stateVector, sequence);
}
