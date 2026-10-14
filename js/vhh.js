// VHH (camelid single-domain) hallmark detection — port of is_vhh from
// unplage/humanization scripts/humanize/numbering.py (MIT License).
// Hallmarks at Kabat positions 37/44/45/47 in FR2.

const VHH_HALLMARK = {
  37: ['F', 'Y'],
  44: ['E', 'Q'],
  45: ['R'],
  47: ['G', 'F', 'S'],
};

/**
 * @param {{ pos: string, aa: string }[]} residues Kabat-numbered heavy residues
 * @returns {{ isVhh: boolean, matched: number, positions: string[] }}
 */
export function isVhh(residues) {
  const byPos = new Map(residues.map((r) => [r.pos, r.aa.toUpperCase()]));
  const matched = [];
  for (const [posStr, allowed] of Object.entries(VHH_HALLMARK)) {
    const pos = Number(posStr);
    const aa = byPos.get(`H${pos}`);
    if (aa && allowed.includes(aa)) matched.push(`H${pos}${aa}`);
  }
  return { isVhh: matched.length >= 3, matched: matched.length, positions: matched };
}
