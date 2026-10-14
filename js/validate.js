// Input validation and preprocessing.

const STANDARD_AA = new Set('ACDEFGHIKLMNPQRSTVWY'.split(''));
const AMBIGUOUS_AA = new Set(['X', 'B', 'Z', 'U', 'O', 'J']);
const DNA_RE = /^[ACGTUN]+$/i;

const KNOWN_TAGS = [
  { name: 'poly-His (His6)', re: /H{4,}/g },
  { name: 'FLAG', re: /DYKDDDDK/g },
  { name: 'HA', re: /YPYDVPDYA/g },
  { name: 'myc', re: /EQKLISEEDL/g },
  { name: 'V5', re: /GKPIPNPLLGLDST/g },
  { name: 'Strep-tag II', re: /WSHPQFEK/g },
];

export const ERROR_CODES = {
  EMPTY: 'EMPTY',
  DNA: 'DNA_INPUT',
  INVALID_CHARS: 'INVALID_CHARS',
};

/**
 * Validate and normalize a raw sequence.
 * @returns {{ ok:boolean, sequence:string, errors:object[], warnings:string[] }}
 */
export function validateSequence(raw) {
  const errors = [];
  const warnings = [];
  const compact = String(raw || '').replace(/\s+/g, '').toUpperCase();

  if (compact.length === 0) {
    errors.push({ code: ERROR_CODES.EMPTY, message: 'Empty sequence.' });
    return { ok: false, sequence: '', errors, warnings };
  }

  if (DNA_RE.test(compact) && compact.length >= 30 && /[ACGT]/.test(compact)) {
    errors.push({
      code: ERROR_CODES.DNA,
      message: 'Input looks like a nucleotide (DNA/RNA) sequence. Provide an amino-acid sequence.',
    });
    return { ok: false, sequence: compact, errors, warnings };
  }

  const bad = new Set();
  for (const ch of compact) if (!/[A-Z]/.test(ch)) bad.add(ch);
  if (bad.size > 0) {
    errors.push({
      code: ERROR_CODES.INVALID_CHARS,
      message: `Non-letter characters in sequence: ${[...bad].join(' ')}`,
    });
  }

  const ambiguous = new Set();
  for (const ch of compact) {
    if (AMBIGUOUS_AA.has(ch)) ambiguous.add(ch);
    else if (!STANDARD_AA.has(ch) && /[A-Z]/.test(ch)) ambiguous.add(ch);
  }
  if (ambiguous.size > 0) {
    warnings.push(`Non-standard residue(s) present: ${[...ambiguous].join(' ')} (numbering may be unreliable).`);
  }

  // N/C-terminal tag detection.
  for (const tag of KNOWN_TAGS) {
    tag.re.lastIndex = 0;
    if (tag.re.test(compact.slice(0, 25)) || (tag.re.lastIndex = 0, tag.re.test(compact.slice(-25)))) {
      warnings.push(`Possible ${tag.name} tag detected near a terminus; confirm it should be numbered as part of the domain.`);
    }
  }
  // Signal-peptide-like hydrophobic leader after the initiator Met.
  if (compact.length > 25 && compact[0] === 'M' && /[AILVFM]{8,}/.test(compact.slice(1, 21))) {
    warnings.push('Input begins with a hydrophobic stretch resembling a signal peptide; numbering starts at the first residue given.');
  }

  return { ok: errors.length === 0, sequence: compact, errors, warnings };
}
