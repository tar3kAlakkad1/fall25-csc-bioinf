from utils import parse_fasta
from alignment import (
    global_align,
    local_align,
    semi_global_align,
    affine_global_align,
)


def _demo(paths: list[str]) -> None:
    for p in paths:
        parsed = parse_fasta(p)
        print(f"{p} -> {len(parsed)} records")
        for k, v in parsed.items():
            print(f"  {k}: {len(v)} bp")


def _run(algo: str, fasta1: str, id1: str, fasta2: str, id2: str) -> int:
    # Load sequences
    seqs1 = parse_fasta(fasta1)
    seqs2 = parse_fasta(fasta2)
    if id1 not in seqs1:
        raise SystemExit(f"ID '{id1}' not found in {fasta1}")
    if id2 not in seqs2:
        raise SystemExit(f"ID '{id2}' not found in {fasta2}")
    seq1 = seqs1[id1]
    seq2 = seqs2[id2]

    # Dispatch
    name = algo.lower()
    if name in ("global", "needleman", "needleman-wunsch", "nw"):
        return global_align(seq1, seq2)
    if name in ("local", "smith-waterman", "sw"):
        return local_align(seq1, seq2)
    if name in ("semi-global", "semiglobal", "fitting", "fit"):
        return semi_global_align(seq1, seq2)
    if name in ("affine", "affine-global", "affine_global"):
        return affine_global_align(seq1, seq2)

    raise SystemExit(f"Unknown algorithm: {algo}")


if __name__ == "__main__":
    # Modes:
    #  - demo: python main.py <fasta1> [<fasta2> ...]
    #  - run:  python main.py run <algo> <fasta1> <id1> <fasta2> <id2>
    import sys

    args = sys.argv[1:]
    if not args:
        print("Usage: python main.py run <algo> <fasta1> <id1> <fasta2> <id2>")
        sys.exit(0)
    if args[0] == "run":
        if len(args) != 6:
            raise SystemExit(
                "Usage: python main.py run <algo> <fasta1> <id1> <fasta2> <id2>"
            )
        # Execute without printing to keep benchmarking output clean
        _run(args[1], args[2], args[3], args[4], args[5])
    else:
        _demo(args)
