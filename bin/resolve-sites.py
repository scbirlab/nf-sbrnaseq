
#!/usr/bin/env python
"""Expand annotation-relative site rules into genomic sites."""

import argparse
import sys

from carabiner import print_err
import numpy as np
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
        GffFile.from_file(args.gff)
        .as_dict()
    )
    rules = pd.read_csv(args.sites, sep="\t")

    required = {
        args.join_on, 
        "position", 
        "ref_base", 
        "alt_base", 
        "modification_name",
    }
    missing = required - set(rules.columns)
    if missing:
        raise ValueError(f"Missing site columns: {missing}")

    if args.join_on not in annotations.columns:
        raise ValueError(
            f"GFF lacks join column: {args.join_on}"
        )

    rules["position"] = rules["position"].astype(int)

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
        how="inner",
        suffixes=("", "_rule"),
    ).dropna(axis=1, how="all")

    if merged.empty:
        raise ValueError("No GFF features matched site rules")

    if not merged["strand"].isin(["+", "-"]).all():
        raise ValueError("Matched GFF features lack +/- strand")
    
    feature_id = [
        "Name",
        "gene",
        "locus_tag",
        "old_locus_tag",
        "gene_synonym",
        "ID",
    ]
    feature_id = next(_id for _id in feature_id if _id in merged)
    merged = (
        merged
        .assign(
            contig=lambda x: x["seqid"],
            mod_feature_name=lambda x: x[feature_id],
            start=lambda x: x["start"].astype(int),
            end=lambda x: x["end"].astype(int),
            length=lambda x: x["end"] - x["start"] + 1,
            inside=lambda x: x["position"] <= x["length"],
        )
    )

    if not merged["inside"].all():
        invalid = merged.query("not inside")
        raise ValueError(
            "Site position outside annotated feature:\n"
            + invalid.to_string(index=False)
        )

    merged = (
        merged
        .assign(
            pos1=lambda x: np.where(
                x["strand"] == "+",
                x["start"] + x["position"] - 1,
                x["end"] - x["position"] + 1,
            ),
            genomic_ref=lambda x: [
                streq.reverse_complement(base) if strand == "-"
                else base
                for base, strand in zip(
                    x["ref_base"],
                    x["strand"]
                )
            ],
            genomic_alt=lambda x: [
                streq.reverse_complement(base) if strand == "-"
                else base
                for base, strand in zip(
                    x["alt_base"],
                    x["strand"]
                )
            ],
        )  
    )

    # Reference validation against the exact assembly.
    with pysam.FastaFile(args.fasta) as fasta:
        actual = [
            fasta.fetch(contig, int(pos) - 1, int(pos)).upper()
            for contig, pos in zip(
                merged["contig"], merged["pos1"]
            )
        ]

    merged = (
        merged
        .assign(
            actual_base=actual,
            actual_matches_ref=lambda x: [
                obs == exp
                for obs, exp in zip(
                    x["actual_base"],
                    x["genomic_ref"],
                )
            ]
        )
    )

    n_skipped = len(merged["actual_matches_ref"]) - sum(merged["actual_matches_ref"])
    if n_skipped:
        print_err(
            f"Skipped {n_skipped} sites with a different "
            "reference nucleotide",
        )

    merged = merged.query("actual_matches_ref")

    if merged.empty:
        raise ValueError(
            "No candidate sites match the requested reference base"
        )

    # Stable and unique within this reference annotation.
    merged = (
        merged
        .assign(
            site_id=lambda x: [
                f"{contig}:{feat_name}:{int(pos)}:{mod_name}"
                for (contig, feat_name, pos, mod_name) in zip(
                        x["contig"],
                        x["mod_feature_name"],
                        x["pos1"],
                        x["modification_name"],
                    )
            ],
        )
    )
    merged.to_csv(args.output, sep="\t", index=False)

    print_err(
        f"Resolved {len(merged)} genomic sites",
    )


if __name__ == "__main__":
    main()
