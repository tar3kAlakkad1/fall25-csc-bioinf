from typing import Dict, List, Optional


def parse_fasta(filename: str) -> Dict[str, str]:
    """
    Parse a FASTA file into a mapping from sequence identifier to sequence.

    - Supports multi-FASTA files with sequences spanning multiple lines.
    - Preserves sequence character casing without modification.
    - Uses the first whitespace-delimited token on the header line (after '>')
      as the sequence identifier. This is conventional for FASTA where the
      first token is the identifier and the remainder is a description.
    """
    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    current_chunks: List[str] = []

    with open(filename, "r", encoding="utf-8") as fasta_file:
        for raw_line in fasta_file:
            line = raw_line.rstrip("\n\r")
            if not line:
                # Ignore empty lines between records or within sequences
                continue

            if line.startswith(">"):
                # Commit the previous record before starting a new one
                if current_id is not None:
                    sequences[current_id] = "".join(current_chunks)
                # Extract identifier (first token after '>')
                header_body = line[1:].strip()
                identifier = header_body.split()[0] if header_body else ""
                current_id = identifier
                current_chunks = []
            else:
                # Sequence line; collect without whitespace
                current_chunks.append(line.strip())

    # Commit the final record if present
    if current_id is not None:
        sequences[current_id] = "".join(current_chunks)

    return sequences
