import sys
from py_nucflag import run_nucflag

def main():
    if len(sys.argv) != 5:
        raise ValueError(f"Usage: {sys.argv[0]} <bam> <bed> <config> <threads>")

    bam = sys.argv[1]
    bed = sys.argv[2]
    config = sys.argv[3]
    threads = int(sys.argv[4])

    results = run_nucflag(bam, bed, threads, config)

    for res in results:
        print(res.regions)


if __name__ == "__main__":
    raise SystemExit(main())