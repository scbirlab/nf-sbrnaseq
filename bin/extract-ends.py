
#!/usr/bin/env python3
"""Extract one protocol-selected genomic end per molecule."""

import argparse
import csv
import gzip
from pathlib import Path

import pysam


def transcript_strand(read, strandedness):
    """FeatureCounts-style strandedness: 1 forward, 2 reverse."""
    forward = not read.is_reverse

    # In a paired-end library, R2 is opposite R1.
    if read.is_paired and read.is_read2:
        forward = not forward

    if strandedness == 2:
        forward = not forward

    return "+" if forward else "-"


def selected_boundary(read, boundary):
    """Return 1-based position and terminal soft-clip length.

    SAM CIGAR order is reference-oriented.
    For reverse alignments the sequenced 5' end is at
    the right-hand reference boundary.
    """
    use_left = (boundary == "5") != read.is_reverse
    pos1 = (
        read.reference_start + 1
        if use_left else read.reference_end
    )

    cigar = read.cigartuples or []
    terminal = cigar if use_left else cigar[::-1]

    # Skip terminal hard clips to inspect adjacent soft clips.
    clips = 0
    for operation, length in terminal:
        if operation == 5:  # hard clip
            continue
        if operation == 4:  # soft clip
            clips += length
        break

    return pos1, clips


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--read", type=int, choices=(0, 1, 2),
                        default=1)
    parser.add_argument("--boundary", choices=("5", "3"),
                        default="5")
    parser.add_argument("--strand", type=int, choices=(1, 2),
                        default=1)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--exclude-softclip",
                        action="store_true")
    args = parser.parse_args()

    columns = [
        "sample_id", "cell_barcode", "molecule_id",
        "contig", "pos1", "strand", "mapq",
        "terminal_softclip",
    ]

    seen = set()
    emitted = 0
    skipped_softclip = 0

    Path(args.output).parent.mkdir(
        parents=True, exist_ok=True
    )

    with pysam.AlignmentFile(args.bam, "rb") as bam, \
         gzip.open(args.output, "wt", newline="") as output:

        writer = csv.DictWriter(
            output, fieldnames=columns, delimiter="\t"
        )
        writer.writeheader()

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary:
                continue
            if read.is_supplementary or read.is_duplicate:
                continue
            if read.is_qcfail:
                continue
            if read.mapping_quality < args.min_mapq:
                continue

            if args.read == 1 and (
                not read.is_paired or not read.is_read1
            ):
                continue
            if args.read == 2 and (
                not read.is_paired or not read.is_read2
            ):
                continue
            if args.read == 0 and read.is_paired:
                continue

            if not read.has_tag("CB"):
                raise ValueError(
                    f"Missing CB tag: {read.query_name}"
                )

            cell = read.get_tag("CB")
            if not cell:
                raise ValueError(
                    f"Empty CB tag: {read.query_name}"
                )

            position, softclip = selected_boundary(
                read, args.boundary
            )

            if args.exclude_softclip and softclip:
                skipped_softclip += 1
                continue

            if position is None or position < 1:
                continue

            molecule = (cell, read.query_name)
            if molecule in seen:
                raise ValueError(
                    "Repeated molecule identity among eligible "
                    f"reads: {molecule}"
                )
            seen.add(molecule)

            writer.writerow({
                "sample_id": args.sample_id,
                "cell_barcode": cell,
                "molecule_id": read.query_name,
                "contig": bam.get_reference_name(
                    read.reference_id
                ),
                "pos1": position,
                "strand": transcript_strand(
                    read, args.strand
                ),
                "mapq": read.mapping_quality,
                "terminal_softclip": softclip,
            })
            emitted += 1

    print(
        f"Emitted {emitted} molecule ends; "
        f"excluded {skipped_softclip} "
        "soft-clipped terminal alignments"
    )


if __name__ == "__main__":
    main()
