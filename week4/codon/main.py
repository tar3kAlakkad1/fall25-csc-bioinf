import sys
import time
from alignment import (
    global_align,
    local_align,
    semi_global_align,
    affine_global_align,
)


def parse_fasta(filename: str):
    # Simple FASTA parser returning a dict[str, str]
    sequences = {"": ""}  # placeholder to help Codon infer dict types
    with open(filename, "r") as fh:
        current_id = ""
        chunks = []
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id != "":
                    sequences[current_id] = "".join(chunks)
                header = line[1:].strip()
                parts = header.split()
                current_id = parts[0] if len(parts) > 0 else ""
                chunks = []
            else:
                chunks.append(line.strip())
        if current_id != "":
            sequences[current_id] = "".join(chunks)
    if "" in sequences:
        del sequences[""]
    return sequences


def run(args):
    if len(args) < 6:
        sys.exit(1)
    mode = args[0]
    if mode not in ("run", "time"):
        sys.exit(1)
    _, algo, fasta1, id1, fasta2, id2 = args

    seqs1 = parse_fasta(fasta1)
    seqs2 = parse_fasta(fasta2)
    if id1 not in seqs1 or id2 not in seqs2:
        sys.exit(2)

    seq1 = seqs1[id1]
    seq2 = seqs2[id2]

    name = algo.lower()
    if mode == "run":
        if name in ("global", "needleman", "needleman-wunsch", "nw"):
            global_align(seq1, seq2)
        elif name in ("local", "smith-waterman", "sw"):
            local_align(seq1, seq2)
        elif name in ("semi-global", "semiglobal", "fitting", "fit"):
            semi_global_align(seq1, seq2)
        elif name in ("affine", "affine-global", "affine_global"):
            affine_global_align(seq1, seq2)
        else:
            sys.exit(3)
    else:
        start = time.perf_counter()
        if name in ("global", "needleman", "needleman-wunsch", "nw"):
            global_align(seq1, seq2)
        elif name in ("local", "smith-waterman", "sw"):
            local_align(seq1, seq2)
        elif name in ("semi-global", "semiglobal", "fitting", "fit"):
            semi_global_align(seq1, seq2)
        elif name in ("affine", "affine-global", "affine_global"):
            affine_global_align(seq1, seq2)
        else:
            sys.exit(3)
        end = time.perf_counter()
        # print milliseconds only
        print(int((end - start) * 1000))


if __name__ == "__main__":
    run(sys.argv[1:])
