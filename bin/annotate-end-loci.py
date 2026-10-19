
#!/usr/bin/env python
"""Annotate RNA-end loci using BED genomic intervals."""

import argparse
import csv
import sys


def annotate(result, output):

    columns = [
        "locus_id", "contig", "start1", "end1",
        "strand", "support",
        "feature_type", "feature_strand",
        "feature_attributes", "overlap_bp",
        "relative_orientation",
        "distance_to_feature_5prime",
    ]
    writer = csv.DictWriter(
        output, 
        fieldnames=columns, 
        delimiter="\t",
    )
    writer.writeheader()

    for line in result:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 13:
            raise ValueError(
                f"Unexpected number ({len(parts)}) of bedtools columns: {parts}"
            )

        (
            chrom, start0, end0,
            locus_id, support, strand,
            fchrom, fstart0, fend0,
            attrs, _score, fstrand, overlap,
        ) = parts

        start0, end0 = int(start0), int(end0)
        overlap = int(overlap)
        orientation = (
            "sense" if strand == fstrand
            else "antisense"
        )

        writer.writerow({
            "locus_id": locus_id,
            "contig": chrom,
            "start1": start0 + 1,
            "end1": end0,
            "strand": strand,
            "support": support,
            "feature_strand": fstrand,
            "feature_attributes": attrs,
            "relative_orientation": orientation,
            "distance_to_feature_5prime": overlap,
        })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", 
        default=sys.stdin, 
        type=argparse.FileType("r"), 
        nargs="?",
    )
    parser.add_argument(
        "--output", 
        default=sys.stdout, 
        type=argparse.FileType("w"),
    )
    args = parser.parse_args()

    annotate(args.input, args.output)


if __name__ == "__main__":
    main()
