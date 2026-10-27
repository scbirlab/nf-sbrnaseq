
#!/usr/bin/env python3
"""Test eggNOG-to-AnnData identifier consistency."""

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "bin/emapper-gene-table.py"

GFF = """##gff-version 3
#!genome-build ASM584v2
#!genome-build-accession NCBI_Assembly:GCF_000005845.2
NC_000913.3	NCBI	gene	1	100	.	+	.	ID=gene-b0001;locus_tag=b0001
NC_000913.3	NCBI	CDS	1	100	.	+	0	ID=cds-NP_1;Parent=gene-b0001;protein_id=NP_1
NC_000913.3	NCBI	gene	200	300	.	+	.	ID=gene-b0002;locus_tag=b0002
NC_000913.3	NCBI	mRNA	200	300	.	+	.	ID=rna2;Parent=gene-b0002
NC_000913.3	NCBI	CDS	200	300	.	+	0	ID=cds-NP_2;Parent=rna2;protein_id=NP_2
#!genome-build ASM19595v2
#!genome-build-accession NCBI_Assembly:GCF_000195955.2
NC_000962.3	NCBI	gene	1	100	.	+	.	ID=gene-Rv0001;locus_tag=Rv0001
NC_000962.3	NCBI	CDS	1	100	.	+	0	ID=cds-NP_3;Parent=gene-Rv0001;protein_id=NP_3
"""

EMAP = """## emapper-2.1.13
#query	COG_category	EC	KEGG_Pathway
NP_1	JK	1.1.1.3,2.7.2.4	ko00260,map00260
NP_2	E	-	ko00300
NP_3	L	2.7.7.7	ko03030
"""

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    gff = tmp / "test.gff"
    emap = tmp / "emap.tsv"
    output = tmp / "functional.tsv"

    gff.write_text(GFF)
    emap.write_text(EMAP)

    subprocess.run([
        sys.executable, str(SCRIPT),
        "--gff", str(gff),
        "--emapper", str(emap),
        "--output", str(output),
    ], check=True)

    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(rows) == 3

    def composite(row):
        return ":".join(row[key] for key in (
            "genome_accession", "assembly", "chr", "gene_id"
        ))

    by_id = {composite(row): row for row in rows}

    eco = "GCF_000005845.2:ASM584v2:NC_000913.3:b0001"
    mtb = "GCF_000195955.2:ASM19595v2:NC_000962.3:Rv0001"

    assert eco in by_id
    assert mtb in by_id
    assert by_id[eco]["gff_gene_id"] == "gene-b0001"
    assert by_id[eco]["gene_id"] == "b0001"
    assert by_id[eco]["COG_category"] == "JK"
    assert by_id[eco]["EC"] == "1.1.1.3,2.7.2.4"
    assert by_id[eco]["KEGG_Pathway"] == "ko00260"
    assert by_id[mtb]["COG_category"] == "L"

print("PASS: assembly metadata, locus tags and eggNOG joins")
