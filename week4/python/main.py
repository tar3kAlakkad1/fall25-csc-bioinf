from utils import parse_fasta


def _demo(paths: list[str]) -> None:
	for p in paths:
		parsed = parse_fasta(p)
		print(f"{p} -> {len(parsed)} records")
		for k, v in parsed.items():
			print(f"  {k}: {len(v)} bp")


if __name__ == "__main__":
	# Quick manual test: run `python main.py data/q1.fa data/t1.fa ...`
	import sys
	if len(sys.argv) > 1:
		_demo(sys.argv[1:])
	else:
		print("Usage: python main.py <fasta1> [<fasta2> ...]")
