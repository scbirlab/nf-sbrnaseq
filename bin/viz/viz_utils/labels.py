
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


def end_locus_label(identifier, feature_attributes):
    """Label an end locus with its nearest annotated gene."""
    attributes = gff_attributes(feature_attributes)
    symbol = attributes.get("gene", "")

    # NCBI gene symbols generally occur under 'gene'.
    # Name is a fallback when it differs from locus_tag.
    if not symbol:
        symbol = attributes.get("Name", "")

    return gene_label(identifier, symbol, identifier)
