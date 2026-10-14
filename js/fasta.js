// Minimal FASTA parser (multi-record, wrapped lines, .fa/.fasta/.txt).

/**
 * Parse FASTA text into records.
 * @returns {{ id: string, description: string, sequence: string }[]}
 */
export function parseFasta(text) {
  const records = [];
  let current = null;
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      if (current) records.push(current);
      const header = line.slice(1).trim();
      const spaceIdx = header.search(/\s/);
      const id = spaceIdx === -1 ? header : header.slice(0, spaceIdx);
      const description = spaceIdx === -1 ? '' : header.slice(spaceIdx + 1).trim();
      current = { id: id || `seq${records.length + 1}`, description, sequence: '' };
    } else if (current) {
      current.sequence += line.replace(/\s+/g, '');
    } else {
      // Sequence data before any header — treat as a single anonymous record.
      current = { id: 'sequence', description: '', sequence: line.replace(/\s+/g, '') };
    }
  }
  if (current) records.push(current);
  return records.filter((r) => r.sequence.length > 0);
}
