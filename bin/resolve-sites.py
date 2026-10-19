
#!/usr/bin/env python
"""Expand annotation-relative site rules into genomic sites."""

import argparse
import sys

import pandas as pd
import pysam
import streq
from bioino import GffFile


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--sites", required=True)
    p.add_argument("--join-on", default="feature")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    # bioino handles parsing and flattening GFF attributes.
    annotations = pd.DataFrame(
        GffFile.from_file(args.gff).as_dict()
    )
    rules = pd.read_csv(args.sites, sep="\t", dtype=str)

    required = {
        args.join_on, "position", "ref_base", "alt_base"
    }
    missing = required - set(rules.columns)
    if missing:
        raise ValueError(f"Missing site columns: {missing}")

    if args.join_on not in annotations.columns:
        raise ValueError(
            f"GFF lacks join column: {args.join_on}"
        )

    # Only one positional rule column; all other columns
    # are arbitrary user metadata.
    rules["position"] = pd.to_numeric(
        rules["position"], errors="raise"
    ).astype(int)

    if (rules["position"] < 1).any():
        raise ValueError("position must be >= 1")

    for col in ("ref_base", "alt_base"):
        rules[col] = rules[col].str.upper()
        if not rules[col].isin(list("ACGT")).all():
            raise ValueError(
                f"{col} must contain single A/C/G/T bases"
            )

    # GFF is the left table; unmatched annotations
    # are not emitted as nominated sites.
    merged = annotations.merge(
        rules,
        on=args.join_on,
        how="left",
        suffixes=("", "_rule"),
        indicator=True,
    )
    merged = merged.loc[
        merged["_merge"] == "both"
    ].drop(columns="_merge").copy()

    if merged.empty:
        raise ValueError("No GFF features matched site rules")

    if not merged["strand"].isin(["+", "-"]).all():
        raise ValueError("Matched GFF features lack +/- strand")

    merged["start"] = merged["start"].astype(int)
    merged["end"] = merged["end"].astype(int)

    length = merged["end"] - merged["start"] + 1
    inside = merged["position"] <= length

    if not inside.all():
        invalid = merged.loc[
            ~inside, [args.join_on, "position", "start", "end"]
        ]
        raise ValueError(
            "Site position outside annotated feature:\n"
            + invalid.to_string(index=False)
        )

    merged["contig"] = merged["seqid"]
    merged["pos1"] = merged["start"] + merged["position"] - 1

    reverse = merged["strand"] == "-"
    merged.loc[reverse, "pos1"] = (
        merged.loc[reverse, "end"]
        - merged.loc[reverse, "position"] + 1
    )

    # streq owns the sequence orientation conversion.
    def genomic_base(base, strand):
        if strand == "+":
            return base
        return streq.reverse_complement(base)

    merged["genomic_ref"] = [
        genomic_base(base, strand)
        for base, strand in zip(
            merged["ref_base"], merged["strand"]
        )
    ]
    merged["genomic_alt"] = [
        genomic_base(base, strand)
        for base, strand in zip(
            merged["alt_base"], merged["strand"]
        )
    ]

    # Reference validation against the exact assembly.
    with pysam.FastaFile(args.fasta) as fasta:
        actual = [
            fasta.fetch(contig, int(pos) - 1, int(pos)).upper()
            for contig, pos in zip(
                merged["contig"], merged["pos1"]
            )
        ]

    
    # Select sites whose actual genomic base matches the
    # requested reference base in the correct orientation.
    matched = [
        observed == expected
        for observed, expected in zip(
            actual, merged["genomic_ref"]
        )
    ]

    n_skipped = len(matched) - sum(matched)
    if n_skipped:
        print(
            f"Skipped {n_skipped} sites with a different "
            "reference nucleotide",
            file=sys.stderr,
        )

    merged = merged.loc[matched].copy()

    if merged.empty:
        raise ValueError(
            "No candidate sites match the requested reference base"
        )


    # Stable and unique within this reference annotation.
    merged["site_id"] = [
        f"{contig}:{pos}:{strand}:{feature}:{i}"
        for i, (contig, pos, strand, feature) in enumerate(
            zip(
                merged["contig"],
                merged["pos1"],
                merged["strand"],
                merged[args.join_on],
            ),
            start=1,
        )
    ]

    merged.to_csv(args.output, sep="\t", index=False)

    print(
        f"Resolved {len(merged)} genomic sites",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
