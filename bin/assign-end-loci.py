
#!/usr/bin/env python
"""Assign molecule 5' ends to established genomic loci."""

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path
import sys


def read_loci(handle):
    loci = defaultdict(list)
    reader = csv.DictReader(handle, delimiter="\t")
    required = {
        "end5_locus_id", "contig", "strand",
        "start1", "end1",
    }
    if not required.issubset(reader.fieldnames or []):
        raise ValueError("Missing required locus columns")

    for row in reader:
        start = int(row["start1"])
        end = int(row["end1"])
        if start < 1 or end < start:
            raise ValueError(f"Invalid locus: {row}")

        key = (row["contig"], row["strand"])
        loci[key].append(
            (start, end, row["end5_locus_id"])
        )

    lookup = {}
    for key, intervals in loci.items():
        intervals.sort()

        previous_end = 0
        for start, end, locus_id in intervals:
            if start <= previous_end:
                raise ValueError(
                    f"Overlapping loci on {key}: "
                    f"{locus_id}"
                )
            previous_end = end

        lookup[key] = (
            [x[0] for x in intervals],
            intervals,
        )

    return lookup


def locate(lookup, contig, strand, pos1):
    group = lookup.get((contig, strand))
    if group is None:
        return None

    starts, intervals = group
    index = bisect_right(starts, pos1) - 1

    if index < 0:
        return None

    start, end, locus_id = intervals[index]
    return locus_id if start <= pos1 <= end else None


def assign(paths, lookup, counts_output, qc_output):
    counts = Counter()
    qc = defaultdict(Counter)
    seen = set()

    for path in paths:
        with gzip.open(path, "rt", newline="") as handle:
            reader = csv.DictReader(
                handle, 
                delimiter="\t",
            )

            for row in reader:
                sample = row["sample_id"]
                cell = row["cell_barcode"]
                molecule = row["molecule_id"]

                key = (sample, cell)
                identity = key + (molecule,)
                if identity in seen:
                    raise ValueError(
                        f"Duplicate molecule end: {identity}"
                    )
                seen.add(identity)

                qc[key]["total"] += 1

                locus = locate(
                    lookup,
                    row["contig"],
                    row["strand"],
                    int(row["pos1"]),
                )

                if locus is None:
                    qc[key]["unassigned"] += 1
                else:
                    qc[key]["assigned"] += 1
                    counts[key + (locus,)] += 1

    writer = csv.writer(counts_output, delimiter="\t")
    writer.writerow([
        "sample_id", "cell_barcode",
        "end5_locus_id", "molecule_count",
    ])

    for key, n in sorted(counts.items()):
        writer.writerow([*key, n])

    
    writer = csv.writer(qc_output, delimiter="\t")
    writer.writerow([
        "sample_id", "cell_barcode",
        "total_ends_count", "assigned_ends_count",
        "unassigned_ends_count", "assigned_ends_fraction",
    ])

    for (sample, cell), stats in sorted(qc.items()):
        total = stats["total"]
        assigned = stats["assigned"]
        unassigned = stats["unassigned"]

        assert assigned + unassigned == total

        writer.writerow([
            sample, cell, total, assigned, unassigned,
            assigned / total if total else "",
        ])

    print(
        f"Assigned {sum(counts.values())} molecule ends "
        f"to {len({k[2] for k in counts})} loci",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--loci", 
        required=True,
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "--ends", 
        nargs="+", 
        required=True,
    )
    parser.add_argument(
        "--counts", 
        default="counts.tsv",
        type=argparse.FileType("w"),
    )
    parser.add_argument(
        "--qc", 
        default="assignment-qc.tsv",
        type=argparse.FileType("w"),
    )
    args = parser.parse_args()

    lookup = read_loci(args.loci)
    assign(
        args.ends,
        lookup,
        args.counts,
        args.qc,
    )


if __name__ == "__main__":
    main()
