#!/usr/bin/env python
"""Annotate RNA-end loci using BEDTools closest output."""

import argparse
import csv
import sys


COLUMNS = [
    "locus_id",
    "contig",
    "start1",
    "end1",
    "strand",
    "support",
    "feature_strand",
    "feature_attributes",
    "overlap_bp",
    "relative_orientation",
    "distance_to_feature",
    "distance_to_feature_5prime",
]


def annotate(source, output):
    writer = csv.DictWriter(
        output,
        fieldnames=COLUMNS,
        delimiter="\t",
    )
    writer.writeheader()

    for line in source:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 13:
            raise ValueError(
                f"Expected 13 BEDTools columns, got "
                f"{len(fields)}: {fields}"
            )

        (
            chrom, start0, end0,
            locus_id, support, strand,
            feature_chrom, feature_start,
            feature_end, attrs, _score,
            feature_strand, distance,
        ) = fields

        start0 = int(start0)
        end0 = int(end0)

        # No matching same-strand feature.
        missing = (
            feature_chrom == "."
            or feature_start == "-1"
        )

        if missing:
            overlap_bp = 0
            feature_strand = ""
            attrs = ""
            orientation = ""
            distance = ""
            distance_5prime = ""
        else:
            feature_start = int(feature_start)
            feature_end = int(feature_end)

            # BED intervals are half-open.
            overlap_bp = max(
                0,
                min(end0, feature_end)
                - max(start0, feature_start),
            )

            orientation = (
                "sense" if strand == feature_strand
                else "antisense"
            )

            # Representative end position:
            # midpoint of the selected locus.
            midpoint0 = (start0 + end0 - 1) // 2

            # 0-based coordinate of the annotated
            # gene's first base in transcript orientation.
            feature_5prime0 = (
                feature_start
                if feature_strand == "+"
                else feature_end - 1
            )

            # Positive means downstream in the
            # annotated transcript's orientation.
            distance_5prime = (
                midpoint0 - feature_5prime0
            )
            if feature_strand == "-":
                distance_5prime *= -1

            distance = int(distance)

        writer.writerow({
            "end5_locus_id": locus_id,
            "contig": chrom,
            "start1": start0 + 1,
            "end1": end0,
            "strand": strand,
            "support": support,
            "feature_strand": feature_strand,
            "feature_attributes": attrs,
            "overlap_bp": overlap_bp,
            "relative_orientation": orientation,
            "distance_to_feature": distance,
            "distance_to_feature_5prime": distance_5prime,
        })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    parser.add_argument(
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
    )
    args = parser.parse_args()

    annotate(args.input, args.output)


if __name__ == "__main__":
    main()
