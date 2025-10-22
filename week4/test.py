try:
    from python import os as _os
    from python import sys as _sys
except ImportError:
    import os as _os
    import sys as _sys
import time

os = _os
sys = _sys


def run_python(algo, f1, id1, f2, id2):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    py_dir = os.path.join(base_dir, "python")
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)
    import alignment as alignment_py
    import utils as utils_py

    seqs1 = utils_py.parse_fasta(f1)
    seqs2 = utils_py.parse_fasta(f2)
    if id1 not in seqs1:
        raise SystemExit(f"ID '{id1}' not in {f1}")
    if id2 not in seqs2:
        raise SystemExit(f"ID '{id2}' not in {f2}")
    s1 = seqs1[id1]
    s2 = seqs2[id2]

    name = algo.lower()
    start = time.perf_counter()
    if name in ("global", "needleman", "needleman-wunsch", "nw"):
        _ = alignment_py.global_align(s1, s2)
    elif name in ("local", "smith-waterman", "sw"):
        _ = alignment_py.local_align(s1, s2)
    elif name in ("semi-global", "semiglobal", "fitting", "fit"):
        _ = alignment_py.semi_global_align(s1, s2)
    elif name in ("affine", "affine-global", "affine_global"):
        _ = alignment_py.affine_global_align(s1, s2)
    else:
        raise SystemExit(f"Unknown algo: {algo}")
    end = time.perf_counter()
    return int((end - start) * 1000)


def run_codon(algo, f1, id1, f2, id2):
    try:
        from python import subprocess as subprocess
    except ImportError:
        import subprocess
    base_dir = os.path.dirname(os.path.abspath(__file__))
    codon_exe = os.path.join(base_dir, "codon", "align")
    if not os.path.exists(codon_exe):
        raise SystemExit("codon executable not found; run evaluate.sh to build it")
    proc = subprocess.run(
        [codon_exe, "time", algo, f1, id1, f2, id2],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(proc.returncode)
    out = proc.stdout.strip()
    try:
        return int(out)
    except Exception:
        raise SystemExit(f"unexpected codon output: {out!r}")


def main(argv):
    if len(argv) != 7 or argv[0] != "run":
        print("Usage: python test.py run <python|codon> <algo> <f1> <id1> <f2> <id2>")
        sys.exit(2)
    _, lang, algo, f1, id1, f2, id2 = argv
    if lang == "python":
        ms = run_python(algo, f1, id1, f2, id2)
    elif lang == "codon":
        ms = run_codon(algo, f1, id1, f2, id2)
    else:
        raise SystemExit(f"Unknown language: {lang}")
    print(ms)


if __name__ == "__main__":
    main(sys.argv[1:])


