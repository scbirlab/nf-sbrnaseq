
#!/usr/bin/env python3
"""Map eggNOG annotations onto featureCounts/AnnData gene IDs."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from bioino import GffFile
from carabiner import print_err


# Must match anndata-utils.py and featurecounts.nf.
GENE_ID = ("genome_accession", "assembly", "chr", "gene_id")
ANNOTATIONS = (
    "COG_category",
    "EC",
    "KEGG_Pathway",
    "KEGG_Module",
    "BRITE",
)
FIELDS = (
    *GENE_ID, "contig", "gff_gene_id", "locus_tag",
    *ANNOTATIONS,
)


def ec_subclasses(value):
    return tuple(sorted({
        ".".join(ec.split(".")[:2])
        for ec in value.split(",")
        if len(ec.split(".")) >= 2
        and all(x.isdigit() for x in ec.split(".")[:2])
    }))

def reference_metadata(gff):
    """Associate each GFF contig with its NCBI assembly metadata.

    NCBI GFF3 uses:
      #!genome-build ASM584v2
      #!genome-build-accession NCBI_Assembly:GCF_000005845.2

    Metadata may change within concatenated GFF files.
    """
    metadata = {}
    assembly = None
    accession = None

    with open(gff) as handle:
        for line in handle:
            if line.startswith("#!genome-build "):
                assembly = line.split(None, 1)[1].strip()

            elif line.startswith("#!genome-build-accession "):
                value = line.split(None, 1)[1].strip()
                accession = value.split(":", 1)[-1]

            elif line and not line.startswith("#"):
                seqid = line.split("\t", 1)[0]

                if not assembly or not accession:
                    raise ValueError(
                        f"Missing NCBI assembly metadata for {seqid}"
                    )

                key = (accession, assembly)

                if seqid in metadata and metadata[seqid] != key:
                    raise ValueError(
                        f"Contig {seqid} belongs to multiple assemblies"
                    )

                metadata[seqid] = key

    return metadata


def read_emapper(path):
    """Read the commented eggNOG-mapper annotation table."""
    with open(path) as handle:
        header = None

        for line in handle:
            if line.startswith("#query\t"):
                header = line[1:].rstrip("\n").split("\t")
                break

        if header is None:
            raise ValueError("Missing emapper #query header")

        required = {"query", *ANNOTATIONS}
        if not required.issubset(header):
            raise ValueError(
                f"Missing columns: {required - set(header)}"
            )

        reader = csv.DictReader(
            (
                line for line in handle
                if line.strip() and not line.startswith("#")
            ),
            fieldnames=header,
            delimiter="\t",
        )

        result = {}
        for row in reader:
            query = row["query"]

            if query in result:
                raise ValueError(
                    f"Duplicate emapper protein query: {query}"
                )

            result[query] = row

    return result


def parse_terms(value, field):
    """Parse multi-valued eggNOG functional annotations."""
    if not value or value in ("-", "."):
        return set()

    if field == "COG_category":
        return set(value.strip())

    terms = {
        term.strip() for term in value.split(",")
        if term.strip() and term.strip() != "-"
    }

    if field == "KEGG_Pathway":
        # Normalise equivalent map00260 and ko00260 entries.
        terms = {
            "ko" + term[3:] if term.startswith("map") else term
            for term in terms
        }

    return terms


def gene_table(gff, annotations):
    """Join protein-level annotations to locus-tag-labelled genes."""
    contig_metadata = reference_metadata(gff)

    genes = {}
    parents = defaultdict(set)
    proteins = defaultdict(set)

    for row in GffFile.from_file(str(gff)).as_dict():
        seqid = row["seqid"]
        feature = row["feature"]
        identifier = row.get("ID", "")
        parent_ids = [
            p for p in row.get("Parent", "").split(",") if p
        ]

        if identifier:
            key = (seqid, identifier)
            parents[key].update(
                (seqid, p) for p in parent_ids
            )

        if feature == "gene":
            if not identifier:
                raise ValueError(
                    f"Gene without GFF ID on {seqid}"
                )

            key = (seqid, identifier)
            if key in genes:
                raise ValueError(
                    f"Duplicate gene feature: {key}"
                )
            genes[key] = row

        elif feature == "CDS":
            protein = row.get("protein_id", "")

            if protein:
                # Begin traversal at the CDS's parent.
                # CDS identifiers themselves are unnecessary.
                proteins[protein].update(
                    (seqid, p) for p in parent_ids
                )

    def ancestor_genes(start):
        found = set()
        visited = set()
        pending = list(start)

        while pending:
            key = pending.pop()

            if key in visited:
                continue
            visited.add(key)

            if key in genes:
                found.add(key)

            pending.extend(parents.get(key, ()))

        return found

    memberships = defaultdict(
        lambda: {
            field: set() for field in ANNOTATIONS
        }
    )

    matched_proteins = set()
    unmatched_proteins = set()

    for protein, start_nodes in proteins.items():
        if protein not in annotations:
            unmatched_proteins.add(protein)
            continue

        targets = ancestor_genes(start_nodes)

        if not targets:
            raise ValueError(
                f"Cannot resolve gene for protein {protein}"
            )

        matched_proteins.add(protein)
        record = annotations[protein]

        # A protein accession may annotate multiple identical
        # proteins; preserve all unambiguous parent relationships.
        for gene_key in targets:
            for field in ANNOTATIONS:
                memberships[gene_key][field].update(
                    parse_terms(record.get(field, ""), field)
                )

    output = []
    skipped = 0

    for (seqid, gff_id), row in sorted(genes.items()):
        locus_tag = row.get("locus_tag", "")

        # featureCounts uses -g locus_tag. Genes without
        # locus_tag cannot be represented by that convention.
        if not locus_tag:
            skipped += 1
            continue

        accession, assembly = contig_metadata[seqid]
        values = memberships[(seqid, gff_id)]

        output.append({
            "genome_accession": accession,
            "assembly": assembly,
            "chr": seqid,
            "gene_id": locus_tag,
            "contig": seqid,
            "gff_gene_id": gff_id,
            "locus_tag": locus_tag,
            "COG_category": "".join(
                sorted(values["COG_category"])
            ),
            "EC": ",".join(sorted(values["EC"])),
            "KEGG_Pathway": ",".join(
                sorted(values["KEGG_Pathway"])
            ),
        })

    keys = [
        tuple(row[field] for field in GENE_ID)
        for row in output
    ]

    if len(keys) != len(set(keys)):
        raise ValueError(
            "Duplicate composite gene IDs after locus_tag mapping"
        )

    print_err(
        f"GFF genes: {len(genes)}; "
        f"with locus_tag: {len(output)}; "
        f"skipped: {skipped}"
    )
    print_err(
        f"Matched eggNOG proteins: {len(matched_proteins)}; "
        f"GFF proteins without eggNOG hit: "
        f"{len(unmatched_proteins)}"
    )

    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--emapper", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    annotations = read_emapper(args.emapper)
    rows = gene_table(args.gff, annotations)

    if not rows:
        raise ValueError("No locus-tag-labelled genes found")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    for field in ANNOTATIONS:
        n = sum(bool(row[field]) for row in rows)
        print_err(
            f"{field}: {n}/{len(rows)} genes "
            f"({n / len(rows):.1%})"
        )


if __name__ == "__main__":
    main()
