
#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "bin/emapper-gene-table.py"

GFF = """##gff-version 3
chr1	NCBI	gene	1	100	.	+	.	ID=gene1;locus_tag=L1
chr1	NCBI	mRNA	1	100	.	+	.	ID=rna1;Parent=gene1
chr1	NCBI	CDS	1	100	.	+	0	ID=cds1;Parent=rna1;protein_id=NP_1.1
chr1	NCBI	gene	200	300	.	+	.	ID=gene2;locus_tag=L2
chr1	NCBI	CDS	200	300	.	+	0	ID=cds2;Parent=gene2;protein_id=NP_2.1
chr2	NCBI	gene	1	100	.	+	.	ID=gene1;locus_tag=L3
"""

EMAP = """## emapper-2.1.13
#query	COG_category	EC	KEGG_Pathway
NP_1.1	JK	1.1.1.3,2.7.2.4	ko00260,map00260,ko01100
NP_2.1	E	-	ko00300
"""

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    gff = tmp / "test.gff"
    emap = tmp / "emap.tsv"
    output = tmp / "out.tsv"

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
    by_key = {(r["contig"], r["gene_id"]): r for r in rows}

    assert by_key["chr1", "gene1"]["COG_category"] == "JK"
    assert by_key["chr1", "gene1"]["EC"] == "1.1.1.3,2.7.2.4"
    assert by_key["chr1", "gene1"]["KEGG_Pathway"] == (
        "ko00260,ko01100"
    )
    assert by_key["chr1", "gene2"]["COG_category"] == "E"
    assert by_key["chr1", "gene2"]["EC"] == ""
    assert by_key["chr2", "gene1"]["COG_category"] == ""

print("PASS: protein-to-gene joins and functional memberships")
