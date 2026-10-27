
"""Human-readable gene and RNA-locus labels."""

from urllib.parse import unquote


def gene_label(identifier, symbol=None, locus_tag=None):
    """Format a unique, readable label without changing IDs."""
    identifier = str(identifier)
    locus = str(locus_tag or identifier.rsplit(":", 1)[-1])

    if symbol is None:
        return locus

    symbol = str(symbol).strip()
    if symbol in ("", ".", "-", "nan", locus, identifier):
        return locus

    return f"{symbol} ({locus})"


def gff_attributes(value):
    """Decode the attribute string from a GFF3 gene feature."""
    if value is None:
        return {}

    attributes = {}
    for field in str(value).split(";"):
        if "=" in field:
            key, val = field.split("=", 1)
            attributes[key] = unquote(val)
    return attributes



def end_locus_gene(identifier, feature_attributes):
    """Return stable gene grouping ID and human-readable symbol."""
    attrs = gff_attributes(feature_attributes)

    gene_id = (
        attrs.get("ID")
        or attrs.get("locus_tag")
        or ""
    )

    symbol = (
        attrs.get("gene")
        or attrs.get("Name")
        or ""
    )

    if not symbol and gene_id.startswith("gene-"):
        symbol = gene_id.removeprefix("gene-")

    # Unannotated loci must not all form one group.
    if not gene_id:
        gene_id = f"unassigned:{identifier}"

    return gene_id, symbol


def end_locus_label(identifier, feature_attributes):
    """Return symbol plus stable end-locus identifier."""
    _, symbol = end_locus_gene(
        identifier, feature_attributes
    )
    return gene_label(identifier, symbol, identifier)

