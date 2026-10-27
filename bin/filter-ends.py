
#!/usr/bin/env python
"""Filter non-overlapping RNA-end loci from Paraclu."""

import argparse
import csv
from collections import defaultdict
import sys


def parse_candidates(filename, max_width):
    candidates = []

    with open(filename) as handle:
        for line in handle:
            if not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 8:
                    raise ValueError(f"Invalid Paraclu row: {line}")

                chrom, strand = fields[:2]
                start, end, sites, support = map(int, fields[2:6])
                min_density, max_density = map(float, fields[6:8])

                # Paraclu coordinates are inclusive, 1-based.
                width = end - start + 1
                if width > max_width:
                    continue

                candidates.append({
                    "contig": chrom,
                    "strand": strand,
                    "start1": start,
                    "end1": end,
                    "width": width,
                    "n_positions": sites,
                    "support": support,
                    "min_density": min_density,
                    "max_density": max_density,
                })

    return candidates


def select_nonoverlapping(candidates):
    """Greedy selection favouring high local molecule density.

    Nested and partially overlapping candidates cannot both
    enter the final set. Single-position peaks are permitted.
    """
    ranked = sorted(
        candidates,
        key=lambda x: (
            -(x["support"] / x["width"]),
            -x["support"],
            x["width"],
            x["contig"],
            x["strand"],
            x["start1"],
        ),
    )

    chosen = defaultdict(list)
    selected = []

    for row in ranked:
        key = (row["contig"], row["strand"])
        overlaps = any(
            row["start1"] <= old["end1"]
            and old["start1"] <= row["end1"]
            for old in chosen[key]
        )
        if overlaps:
            continue

        chosen[key].append(row)
        selected.append(row)

    return sorted(
        selected,
        key=lambda x: (
            x["contig"], x["strand"], x["start1"]
        ),
    )


def write_loci(rows, filename):
    columns = [
        "end5_locus_id", "contig", "strand",
        "start1", "end1", "width",
        "n_positions", "support",
        "min_density", "max_density",
    ]

    with open(filename, "w", newline="") as out:
        writer = csv.DictWriter(
            out, delimiter="\t", fieldnames=columns
        )
        writer.writeheader()

        for row in rows:
            row = dict(row)
            row["end5_locus_id"] = (
                f'{row["contig"]}:'
                f'{row["start1"]}-{row["end1"]}:'
                f'{row["strand"]}'
            )
            writer.writerow(row)


def write_bed(rows, filename):
    with open(filename, "w") as out:
        for row in rows:
            out.write(
                f'{row["contig"]}\t'
                f'{row["start1"] - 1}\t'
                f'{row["end1"]}\t'
                f'{row["contig"]}:'
                f'{row["start1"]}-{row["end1"]}:'
                f'{row["strand"]}\t'
                f'{row["support"]}\t'
                f'{row["strand"]}\n'
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", default=sys.stdin, nargs="?")
    parser.add_argument("--max-width", type=int, default=200)
    args = parser.parse_args()

    if args.max_width < 1:
        parser.error("Width must be positive")

    candidates = parse_candidates(
        args.input, 
        args.max_width,
    )
    selected = select_nonoverlapping(candidates)

    write_loci(
        selected, "loci.tsv"
    )
    write_bed(
        selected, "loci.bed"
    )

    print(
        f"{len(candidates)} candidates; "
        f"{len(selected)} non-overlapping loci"
    )


if __name__ == "__main__":
    main()
