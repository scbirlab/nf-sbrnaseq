
#!/usr/bin/env python3
"""Count transcript-oriented bases at resolved sites per cell."""

import argparse
import csv
import gzip
from collections import Counter, defaultdict

import pysam
import streq

BASES = "ACGT"

COUNT_COLUMNS = [
    "sample_id", "cell_barcode",
    "count_A", "count_C", "count_G", "count_T",
    "ref_count", "alt_count", "other_count",
    "callable_depth", "alt_fraction",
    "mapping_primary_reads",
    "mapping_confident_reads",
    "mapping_confident_fraction",
    "mapping_multihit_reads",
    "mapping_multihit_fraction",
    "mapping_min_mapq",
]



def load_sites(filename):
    with open(filename, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        required = {
            "site_id", "contig", "pos1", "strand",
            "ref_base", "alt_base",
            "genomic_ref", "genomic_alt",
        }
        missing = required - set(columns)
        if missing:
            raise ValueError(f"Missing site columns: {missing}")

        sites = list(reader)

    for site in sites:
        site["pos1"] = int(float(site["pos1"]))
        if site["pos1"] < 1:
            raise ValueError("Site positions must be 1-based")
        if site["strand"] not in ("+", "-"):
            raise ValueError("Site strand must be + or -")
        for field in (
            "ref_base", "alt_base",
            "genomic_ref", "genomic_alt",
        ):
            site[field] = site[field].upper()
            if site[field] not in BASES:
                raise ValueError(f"Invalid {field}: {site}")

        if site["ref_base"] == site["alt_base"]:
            raise ValueError(
                f"REF and ALT must differ: {site['site_id']}"
            )

        convert = (
            lambda b: b if site["strand"] == "+"
            else streq.reverse_complement(b)
        )
        if convert(site["genomic_ref"]) != site["ref_base"]:
            raise ValueError(f"REF orientation mismatch: {site}")
        if convert(site["genomic_alt"]) != site["alt_base"]:
            raise ValueError(f"ALT orientation mismatch: {site}")

    duplicates = len(sites) - len({
        site["site_id"] for site in sites
    })
    if duplicates:
        raise ValueError(f"{duplicates} duplicate site IDs")

    collisions = set(columns) & set(COUNT_COLUMNS)
    if collisions:
        raise ValueError(
            f"Site metadata conflicts with output columns: {collisions}"
        )

    return columns, sites


def count_site(bam, site, min_mapq, min_baseq):
    """Return per-cell counters; overlapping mates form one call."""
    pos0 = site["pos1"] - 1
    calls = defaultdict(list)
    inspected = 0
    missing_cb = 0

    for read in bam.fetch(site["contig"], pos0, pos0 + 1):
        if read.flag & 3844:
            continue
        if read.mapping_quality < min_mapq:
            continue

        inspected += 1
        if not read.has_tag("CB"):
            missing_cb += 1
            continue

        cb = read.get_tag("CB")
        if not cb:
            missing_cb += 1
            continue

        seq = read.query_sequence
        qual = read.query_qualities
        if seq is None or qual is None:
            continue

        for qpos, rpos in read.get_aligned_pairs(
            matches_only=True
        ):
            if rpos != pos0:
                continue

            base = seq[qpos].upper()
            quality = int(qual[qpos])
            if base not in BASES or quality < min_baseq:
                break

            if site["strand"] == "-":
                base = streq.reverse_complement(base)

            # Same QNAME in one cell represents one fragment.
            calls[(cb, read.query_name)].append(
                (base, quality)
            )
            break

    if missing_cb:
        raise ValueError(
            f"{missing_cb}/{inspected} passing alignments "
            f"at {site['site_id']} lack CB tags"
        )

    counts = defaultdict(Counter)
    for (cb, _name), observations in calls.items():
        highest_quality = max(q for _, q in observations)
        best_bases = {
            base for base, quality in observations
            if quality == highest_quality
        }

        # Equal-quality disagreement between mates:
        # don't guess the molecule's base.
        if len(best_bases) != 1:
            continue

        counts[cb][next(iter(best_bases))] += 1

    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sites", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=25)

    parser.add_argument("--raw-bam", required=True)
    parser.add_argument(
        "--mapper",
        required=True,
        choices=("star", "bowtie2", "minimap2"),
    )
    parser.add_argument(
        "--mapping-min-mapq", type=int, required=True
    )
    parser.add_argument(
        "--allow-multimap",
        type=lambda x: x.lower() == "true",
        required=True,
    )

    args = parser.parse_args()

    if args.min_mapq < 0 or args.min_baseq < 0:
        parser.error("Quality thresholds must be nonnegative")

    metadata_columns, sites = load_sites(args.sites)
    output_columns = COUNT_COLUMNS + metadata_columns

    nrows = 0
    
    with pysam.AlignmentFile(args.bam, "rb") as bam, \
         pysam.AlignmentFile(args.raw_bam, "rb") as raw_bam, \
         gzip.open(args.output, "wt", newline="") as handle:

        writer = csv.DictWriter(
            handle, fieldnames=output_columns, delimiter="\t"
        )
        writer.writeheader()

        for site in sites:
            if site["contig"] not in bam.references:
                raise ValueError(
                    f"Missing BAM contig: {site['contig']}"
                )

            counts = count_site(
                bam, site, args.min_mapq, args.min_baseq
            )
            qc = mapping_qc(
                raw_bam,
                site,
                args.mapper,
                args.mapping_min_mapq,
                args.allow_multimap,
            )

            for cell, bases in sorted(counts.items()):
                ref = bases[site["ref_base"]]
                alt = bases[site["alt_base"]]
                depth = sum(bases.values())

                row = {
                    "sample_id": args.sample_id,
                    "cell_barcode": cell,
                    "count_A": bases["A"],
                    "count_C": bases["C"],
                    "count_G": bases["G"],
                    "count_T": bases["T"],
                    "ref_count": ref,
                    "alt_count": alt,
                    "other_count": depth - ref - alt,
                    "callable_depth": depth,
                    "alt_fraction": (
                        alt / (ref + alt)
                        if ref + alt > 0 else ""
                    ),
                    **qc,
                    **site,
                }
                writer.writerow(row)
                nrows += 1

    print(
        f"Counted {len(sites)} sites; "
        f"wrote {nrows} covered cell-site rows"
    )



def mapping_qc(
    raw_bam, site, mapper, min_mapq, allow_multimap
):
    """Observed read-alignment confidence at a genomic site.

    Counts alignment records, not deduplicated molecules.
    Excludes unmapped, secondary, QC-failed and
    supplementary alignments.
    """
    pos0 = site["pos1"] - 1

    primary = 0
    confident = 0
    multihit = 0
    nh_known = 0

    for read in raw_bam.fetch(
        site["contig"], pos0, pos0 + 1
    ):
        if read.flag & 2820:
            # 0x4 unmapped, 0x100 secondary,
            # 0x200 QC fail, 0x800 supplementary
            continue

        # Only alignments with a base at the site.
        if not any(
            refpos == pos0
            for _, refpos in read.get_aligned_pairs(
                matches_only=True
            )
        ):
            continue

        primary += 1

        nh = read.get_tag("NH") if read.has_tag("NH") else None
        if nh is not None:
            nh_known += 1
            if nh > 1:
                multihit += 1

        passes = read.mapping_quality >= min_mapq

        # Match the current downstream filter policy.
        if mapper == "star" and not allow_multimap:
            passes = passes and nh == 1

        if passes:
            confident += 1

    return {
        "mapping_primary_reads": primary,
        "mapping_confident_reads": confident,
        "mapping_confident_fraction": (
            confident / primary if primary else ""
        ),
        "mapping_multihit_reads": multihit,
        "mapping_multihit_fraction": (
            multihit / nh_known if nh_known else ""
        ),
        "mapping_min_mapq": min_mapq,
    }


if __name__ == "__main__":
    main()
