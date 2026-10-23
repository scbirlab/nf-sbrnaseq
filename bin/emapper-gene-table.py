
#!/usr/bin/env python3
"""Map eggNOG-mapper annotations onto NCBI GFF gene IDs."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from bioino import GffFile
from carabiner import print_err


FIELDS = (
    "contig", "gene_id", "locus_tag",
    "COG_category", "EC", "KEGG_Pathway",
)

ANNOTATIONS = ("COG_category", "EC", "KEGG_Pathway")


def read_emapper(path):
    """Read emapper annotations with their commented header."""
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
            (line for line in handle
             if line.strip() and not line.startswith("#")),
            fieldnames=header,
            delimiter="\t",
        )

        result = {}
        for row in reader:
            query = row["query"]
            if query in result:
                print_err(f"Duplicate emapper query: {query}")
            result[query] = row

        return result


def gene_table(gff, annotations):
    """Resolve protein -> CDS -> parent gene relationships."""
    records = list(GffFile.from_file(str(gff)).as_dict())

    nodes = {}
    parents = {}
    genes = {}
    proteins = defaultdict(set)

    for row in records:
        seqid = row["seqid"]
        identifier = row.get("ID", "")

        if identifier:
            key = (seqid, identifier)
            if key in nodes:
                print_err(f"Duplicate feature ID: {key}")
            nodes[key] = row
            parents[key] = [
                (seqid, p)
                for p in row.get("Parent", "").split(",") if p
            ]

        if row["feature"] == "gene":
            if not identifier:
                raise ValueError("Gene without ID")
            genes[(seqid, identifier)] = row

        if row["feature"] == "CDS":
            protein = row.get("protein_id", "")
            if protein and identifier:
                proteins[protein].add((seqid, identifier))

    def ancestor_genes(node):
        """Find all ancestral gene nodes."""
        stack = [node]
        visited = set()
        found = set()

        while stack:
            key = stack.pop()

            if key in visited:
                continue
            visited.add(key)

            if key in genes:
                found.add(key)

            stack.extend(parents.get(key, []))

        return found

    memberships = defaultdict(
        lambda: {key: set() for key in ANNOTATIONS}
    )

    for protein, cds_nodes in proteins.items():
        if protein not in annotations:
            continue

        targets = set()
        for cds in cds_nodes:
            targets.update(ancestor_genes(cds))

        if len(targets) != 1:
            raise ValueError(
                f"Protein {protein} maps to {len(targets)} genes"
            )

        gene = next(iter(targets))
        row = annotations[protein]

        for field in ANNOTATIONS:
            value = row.get(field, "")
            if value in ("", "-", "."):
                continue

            if field == "COG_category":
                terms = set(value.strip())
            else:
                terms = set(value.split(","))

            if field == "KEGG_Pathway":
                # 'map' and 'ko' entries duplicate pathway IDs.
                terms = {
                    t for t in terms if t.startswith("ko")
                }

            memberships[gene][field].update(terms)

    output = []
    for key, row in sorted(genes.items()):
        values = memberships[key]
        output.append({
            "contig": key[0],
            "gene_id": key[1],
            "locus_tag": row.get("locus_tag", ""),
            "COG_category": "".join(
                sorted(values["COG_category"])
            ),
            "EC": ",".join(sorted(values["EC"])),
            "KEGG_Pathway": ",".join(
                sorted(values["KEGG_Pathway"])
            ),
        })

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
        raise ValueError("No genes found")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    for field in ANNOTATIONS:
        n = sum(bool(row[field]) for row in rows)
        print(f"{field}: {n}/{len(rows)} genes ({n/len(rows):.1%})")


if __name__ == "__main__":
    main()
