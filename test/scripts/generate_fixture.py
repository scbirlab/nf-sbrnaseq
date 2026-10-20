#!/usr/bin/env python3
"""Deterministic paired-end sbRNA-seq smoke fixture (stdlib only).

Modes:
    ecoli       E. coli K-12 MG1655 with mRNA, tRNA and rRNA genes
    barnyard    Mixed E. coli MG1655 / B. subtilis 168 cells in one sample

Gene sequences and biotypes are extracted from the exact NCBI RefSeq FASTA/GFF
assemblies used by the pipeline. The first run downloads ~genome-sized cache
files unless supplied locally in --reference-dir. Cache files are deliberately
excluded from the distributed fixture.

A --manifest-only run can produce the deterministic E. coli count design
without downloading a reference; it does not write FASTQs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import random
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
SEED = 42
READ_LENGTH = 70
FRAGMENT_LENGTH = 160

# Genomes are pinned to the same RefSeq assemblies as the Nextflow sample sheet.
REFERENCES = {
    'ecoli': ('GCF_000005845.2', 'ASM584v2', 'NC_000913.3'),
    'bsub': ('GCF_000009045.1', 'ASM904v1', 'NC_000964.3'),
}

# E. coli targets deliberately include 16S and 5S rRNA and two tRNAs.
# Protein genes preserve the original fixture's feature selection.
PROTEIN_TARGETS = {
    'ecoli': ('lacZ', 'lacI', 'lacY', 'rssA', 'tehB'),
    'bsub': ('dnaA', 'gyrA', 'rpoB', 'sigA', 'spo0A'),
}
ECOLI_RNA_TARGETS = ('pheU', 'thrW', 'rrsA', 'rrfA')

BC = [
    ['ACGTGCA', 'TGCATGC', 'GATCCGA', 'CTAGTCT'],
    ['TTGACCA', 'CCATGGT', 'AGCTAGC', 'GGCATTA'],
    ['ATGCCGA', 'CGTATGC', 'TCGGAAT', 'GACTTCC'],
]
LEFT_R1, LEFT_R2 = 'ACGTTCGA', 'TGCAGTCA'
RIGHT_R1, RIGHT_R2 = 'GTCAGTCA', 'CTAGGCTA'

# These read structures match modules/trimming.nf + UMItools_extract.
UMI_PATTERNS = (
    r'^(?P<discard_1>.{8})(?P<cell_1>.{7})(?P<discard_2>.{3}).*',
    r'^(?P<discard_3>.{8})(?P<cell_2>.{7})(?P<cell_3>.{7})(?P<umi_1>.{8})(?P<discard_4>.{2}).*',
)

@dataclass(frozen=True)
class Feature:
    species: str
    chr: str
    gene: str
    locus_tag: str
    biotype: str
    start: int   # zero-based inclusive
    end: int     # zero-based exclusive
    strand: str


def rc(s: str) -> str:
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def check_barcodes() -> None:
    for choices in BC:
        assert len(choices) == len(set(choices))
        assert all(len(x) == 7 for x in choices)
        assert min(hamming(x, y) for i, x in enumerate(choices)
                   for y in choices[i + 1:]) >= 3


def create_barcodes(base: Path) -> None:
    inp = base / 'inputs'
    inp.mkdir(parents=True, exist_ok=True)
    check_barcodes()
    for idx, choices in enumerate(BC, 1):
        with (inp / f'bc{idx}.csv').open('w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['Barcode ID', 'Sequence'])
            for j, seq in enumerate(choices, 1):
                writer.writerow([f'bc{j:02d}', seq])


def reference_names(species: str) -> tuple[str, str, str]:
    accession, assembly, chromosome = REFERENCES[species]
    basename = f'{accession}_{assembly}_genomic'
    directory = '/'.join(accession.split('_')[1].split('.')[0][i:i+3]
                         for i in (0, 3, 6))
    url_root = (f'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{directory}/'
                f'{accession}_{assembly}')
    return basename, url_root, chromosome


def ensure_reference(species: str, cache: Path, local_dir: Path | None,
                     allow_download: bool) -> tuple[Path, Path]:
    basename, url_root, _ = reference_names(species)
    names = (basename + '.fna.gz', basename + '.gff.gz')
    paths = []
    for name in names:
        found = None
        for directory in (local_dir, cache):
            if directory is None:
                continue
            for possibility in (directory / name, directory / name.removesuffix('.gz')):
                if possibility.is_file() and possibility.stat().st_size:
                    found = possibility
                    break
            if found:
                break
        if found is None:
            if not allow_download:
                raise FileNotFoundError(
                    f'{name} is missing. Provide --reference-dir or allow NCBI downloads.'
                )
            cache.mkdir(parents=True, exist_ok=True)
            url = f'{url_root}/{name}'
            target = cache / name
            tmp = cache / (name + '.part')
            print(f'Downloading {url} to {tmp}')
            request = Request(url, headers={'User-Agent': 'nf-sbrnaseq-smoke/2.0'})
            try:
                with urlopen(request, timeout=120) as stream, tmp.open('wb') as fh:
                    while True:
                        data = stream.read(1 << 20)
                        if not data:
                            break
                        fh.write(data)
                if tmp.stat().st_size == 0:
                    raise ValueError(f'Empty download: {url}')
                tmp.replace(target)
            finally:
                tmp.unlink(missing_ok=True)
            found = target
        paths.append(found)
    return paths[0], paths[1]


def open_text(path: Path):
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open('r')


def parse_fasta(path: Path) -> dict[str, str]:
    records, chrom, chunks = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if chrom is not None:
                    records[chrom] = ''.join(chunks).upper()
                chrom, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if chrom is not None:
        records[chrom] = ''.join(chunks).upper()
    if not records:
        raise ValueError(f'No sequences in {path}')
    return records


def parse_attrs(value: str) -> dict[str, str]:
    fields = {}
    for part in value.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            fields[k] = v
    return fields


def parse_features(path: Path, species: str) -> list[Feature]:
    genes = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9 or parts[2] != 'gene':
                continue
            a = parse_attrs(parts[8])
            biotype = a.get('gene_biotype', '')
            if biotype not in ('protein_coding', 'tRNA', 'rRNA'):
                continue
            name = a.get('gene') or a.get('Name') or a.get('locus_tag')
            if not name:
                continue
            genes.append(Feature(
                species=species, chr=parts[0], gene=name,
                locus_tag=a.get('locus_tag', ''), biotype=biotype,
                start=int(parts[3]) - 1, end=int(parts[4]),
                strand=parts[6],
            ))
    if not genes:
        raise ValueError(f'No annotated genes with biotypes found in {path}')
    return genes


def pick_named(genes: list[Feature], names: tuple[str, ...]) -> list[Feature]:
    selected = []
    for name in names:
        candidates = [x for x in genes if x.gene == name]
        if not candidates:
            raise ValueError(f'Expected gene {name} not present in GFF')
        selected.append(sorted(candidates, key=lambda x: (x.start, x.locus_tag))[0])
    return selected


def get_window(gene: Feature, chromosomes: dict[str, str],
               molecule: int, gene_index: int) -> str:
    genome = chromosomes[gene.chr]
    if gene.biotype == 'protein_coding':
        # Fully intragenic; vary coordinates to exercise positional UMI collapse.
        shift = (gene_index * 7 + molecule * 17) % min(
            35, gene.end - gene.start - FRAGMENT_LENGTH + 1)
        if gene.strand == '+':
            start = gene.start + 25 + shift
            if start + FRAGMENT_LENGTH > gene.end:
                start = gene.end - FRAGMENT_LENGTH
            insert = genome[start:start + FRAGMENT_LENGTH]
        else:
            end = gene.end - 25 - shift
            if end - FRAGMENT_LENGTH < gene.start:
                end = gene.start + FRAGMENT_LENGTH
            insert = rc(genome[end - FRAGMENT_LENGTH:end])
    else:
        # RNA genes may be ~76 nt (tRNA) or highly repeated (rRNA).
        # Include a short upstream flank to improve unique mapping and gene
        # assignment. These simulate *precursors*, not mature RNA ends.
        flank = 42 + molecule % 5
        if gene.strand == '+':
            start = gene.start - flank
            insert = genome[start:start + FRAGMENT_LENGTH]
        else:
            end = gene.end + flank
            insert = rc(genome[end - FRAGMENT_LENGTH:end])
    if len(insert) != FRAGMENT_LENGTH or any(c not in 'ACGT' for c in insert):
        raise ValueError(f'Invalid {FRAGMENT_LENGTH}-nt fragment for {gene}')
    return insert


def pick_bsub_rna(genes: list[Feature], chromosomes: dict[str, str]) -> list[Feature]:
    """Pick annotated tRNAs and 5S/16S rRNAs with useful unique flanks.

    Prefer nonoverlapping loci whose 70-bp mate anchors have few identical
    matches in their chromosome; stable sort makes the choice reproducible.
    """
    selected = []
    allgenes = sorted(genes, key=lambda x: x.start)
    genome = chromosomes[REFERENCES['bsub'][2]]
    for biotype, count in (('tRNA', 2), ('rRNA', 2)):
        candidates = [x for x in allgenes if x.biotype == biotype and x.chr == REFERENCES['bsub'][2]]
        if not candidates:
            raise ValueError(f'B. subtilis annotation contains no {biotype} genes')
        ranked = []
        for x in candidates:
            try:
                seq = get_window(x, chromosomes, 0, 0)
            except ValueError:
                continue
            # Lower scores are preferable. Count both strand orientations.
            a, b = seq[:READ_LENGTH], seq[-READ_LENGTH:]
            duplicates = sum(genome.count(s) for s in (a, b, rc(a), rc(b)))
            near = sum(1 for other in allgenes if other is not x and other.chr == x.chr
                       and other.start < x.end + 80 and other.end > x.start - 80)
            # Prefer 5S rRNAs for the first rRNA, longer rRNA for the second.
            sizeclass = 0 if x.end - x.start < 200 else 1
            ranked.append(((near * 10 + duplicates, x.start), sizeclass, x))
        ranked.sort(key=lambda item: item[0])
        if biotype == 'rRNA':
            short = next((v[2] for v in ranked if v[1] == 0), None)
            long = next((v[2] for v in ranked if v[1] == 1), None)
            first = short or (ranked[0][2] if ranked else None)
            second = long or next((v[2] for v in ranked if v[2] != first), None)
            chosen = [x for x in (first, second) if x is not None]
        else:
            chosen = []
            for _, _, x in ranked:
                if x.gene not in {g.gene for g in chosen}:
                    chosen.append(x)
                if len(chosen) == count:
                    break
        if len(chosen) != count:
            raise ValueError(f'Expected {count} B. subtilis {biotype} genes, found {len(chosen)}')
        selected.extend(chosen)
    return selected


def make_sheet(base: Path, mode: str) -> Path:
    sample = 'synthetic_illumina' #if mode == 'ecoli' else 'synthetic_barnyard'
    sheet = base / 'sample-sheet.csv' #if mode == 'ecoli' else 'sample-sheet-barnyard.csv')
    genomes = REFERENCES['ecoli'][0] if mode == 'ecoli' else '+'.join(
        REFERENCES[x][0] for x in ('ecoli', 'bsub')
    )
    with sheet.open('w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['sample_id', 'fastq_pattern', 'genome_accession',
                    'bc1', 'bc2', 'bc3',
                    'adapter_read1_5prime', 'adapter_read2_5prime',
                    'adapter_read1_3prime', 'adapter_read2_3prime',
                    'umi_read1', 'umi_read2'])
        w.writerow([sample, sample, genomes, 'bc1.csv', 'bc2.csv', 'bc3.csv',
                    LEFT_R1, LEFT_R2, RIGHT_R1, RIGHT_R2, *UMI_PATTERNS])
    return sheet


def target_plan(mode: str, gene_data: dict[str, list[Feature]] | None = None):
    if gene_data is None:
        # Offline manifest only: the B. subtilis RNA targets are discovered
        # using the real GFF/FASTA, so this mode is intentionally E. coli only.
        if mode != 'ecoli':
            raise ValueError('--manifest-only is only supported for --mode ecoli')
        return {
            'ecoli': [Feature('ecoli', '', gene, '',
                              'protein_coding' if gene in PROTEIN_TARGETS['ecoli']
                              else 'tRNA' if gene in ECOLI_RNA_TARGETS[:2] else 'rRNA',
                              0, 0, '+')
                      for gene in (*PROTEIN_TARGETS['ecoli'], *ECOLI_RNA_TARGETS)]
        }
    return gene_data


def write_design(
    base: Path, 
    mode: str, 
    plans: dict[str, list[Feature]],
    *, 
    write_reads: bool, 
    chromosomes: dict[str, dict[str, str]] | None = None
):
    sample = 'synthetic_illumina' #if mode == 'ecoli' else 'synthetic_barnyard'
    counts_name = 'design-counts.tsv' #if mode == 'ecoli' else 'design-counts-barnyard.tsv'
    cells_name = 'design-cells.tsv' # if mode == 'ecoli' else 'design-cells-.tsv'
    base.joinpath('fastq').mkdir(exist_ok=True)
    rng = random.Random(SEED + (0 if mode == 'ecoli' else 1))
    expected = Counter()
    pair_counts = Counter()
    label_counts = Counter()
    cell_types = {}
    n_pairs = 0
    n_duplicates = 0
    r1 = r2 = None
    try:
        if write_reads:
            r1 = gzip.open(base / 'fastq' / f'{sample}_R1.fastq.gz', 'wt', compresslevel=6)
            r2 = gzip.open(base / 'fastq' / f'{sample}_R2.fastq.gz', 'wt', compresslevel=6)
        index = 0
        for i, b1 in enumerate(BC[0]):
            for j, b2 in enumerate(BC[1]):
                for k, b3 in enumerate(BC[2]):
                    cell = b1 + b2 + b3
                    # Eight mixed barcodes; the remainder are single species.
                    if mode == 'ecoli':
                        species_present = ('ecoli',)
                    else:
                        slot = index % 8
                        species_present = ('ecoli', 'bsub') if slot == 0 else (
                            ('bsub',) if slot in (1, 2, 3) else ('ecoli',))
                    cell_types[cell] = '+'.join(species_present)
                    for species in species_present:
                        for gidx, feature in enumerate(plans[species]):
                            spidx = 0 if species == 'ecoli' else 1
                            n_molecules = 1 + ((i + 2*j + 3*k + 2*gidx + spidx) % 4)
                            expected[(cell, species, feature.gene, feature.biotype,
                                      REFERENCES[species][0], feature.locus_tag)] += n_molecules
                            label_counts[species, feature.biotype] += n_molecules
                            for mol in range(n_molecules):
                                copies = 2 if (i + j + k + gidx + mol + spidx) % 3 == 0 else 1
                                pair_counts[(cell, species)] += copies
                                if write_reads:
                                    insert = get_window(feature, chromosomes[species], mol, gidx)
                                    u = ''.join(rng.choices('ACGT', k=8))
                                    read1 = LEFT_R1 + b1 + 'AGT' + insert[:READ_LENGTH] + RIGHT_R1
                                    read2 = (LEFT_R2 + b2 + b3 + u + 'GG' +
                                             rc(insert[-READ_LENGTH:]) + RIGHT_R2)
                                    for copy in range(copies):
                                        n_pairs += 1
                                        name = f'SYN{n_pairs:09d}'
                                        for handle, seq in ((r1, read1), (r2, read2)):
                                            handle.write(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n')
                                        if copy:
                                            n_duplicates += 1
                    index += 1
    finally:
        if r1:
            r1.close()
        if r2:
            r2.close()

    with (base / counts_name).open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['cell_barcode', 'gene', 'designed_molecules', 'species',
                    'gene_biotype', 'genome_accession', 'locus_tag'])
        for (cell, species, gene, biotype, acc, locus), n in sorted(expected.items()):
            w.writerow([cell, gene, n, species, biotype, acc, locus])

    with (base / cells_name).open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['cell_barcode', 'designed_species', 'designed_molecules',
                    'paired_fastq_records'])
        for cell, label in sorted(cell_types.items()):
            molecules = sum(n for (bc, *_), n in expected.items() if bc == cell)
            pairs = sum(n for (bc, _), n in pair_counts.items() if bc == cell)
            w.writerow([cell, label, molecules, pairs])

    print(f'{mode}: {len(cell_types)} cells, {sum(expected.values())} designed molecules'
          + (f', {n_pairs} read pairs ({n_duplicates} PCR duplicates)' if write_reads else
             ', manifest only (FASTQs not generated)'))
    print('RNA/protein design:', dict(sorted(label_counts.items())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('ecoli', 'barnyard', 'both'), default='ecoli')
    parser.add_argument('--reference-dir', type=Path,
                        help='Directory containing pinned NCBI FASTA/GFF(.gz) references')
    parser.add_argument('--cache-dir', type=Path, default=HERE / '.reference_cache')
    parser.add_argument('--output-dir', type=Path, default=HERE)
    parser.add_argument('--no-download', action='store_true',
                        help='Require references in cache or --reference-dir')
    parser.add_argument('--manifest-only', action='store_true',
                        help='Generate E. coli design-counts.tsv without FASTQs or references')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    create_barcodes(args.output_dir)

    modes = ('ecoli', 'barnyard') if args.mode == 'both' else (args.mode,)
    if args.manifest_only:
        for mode in modes:
            plan = target_plan(mode)
            make_sheet(args.output_dir, mode)
            write_design(args.output_dir, mode, plan, write_reads=False)
        return

    loaded: dict[str, tuple[dict[str, str], list[Feature]]] = {}
    for species in ('ecoli', 'bsub') if 'barnyard' in modes else ('ecoli',):
        fasta, gff = ensure_reference(species, args.cache_dir, args.reference_dir,
                                      allow_download=not args.no_download)
        genomes = parse_fasta(fasta)
        genes = parse_features(gff, species)
        expected_chr = REFERENCES[species][2]
        if expected_chr not in genomes:
            raise ValueError(f'{species}: {expected_chr} absent from {fasta}')
        if not any(g.chr == expected_chr for g in genes):
            raise ValueError(f'{species}: {expected_chr} absent from GFF')
        loaded[species] = (genomes, genes)

    for mode in modes:
        ec_genomes, ec_genes = loaded['ecoli']
        ecoli_plan = pick_named(ec_genes,
                                 (*PROTEIN_TARGETS['ecoli'], *ECOLI_RNA_TARGETS))
        plans = {'ecoli': ecoli_plan}
        chromosomes = {'ecoli': ec_genomes}
        if mode == 'barnyard':
            bs_genomes, bs_genes = loaded['bsub']
            plans['bsub'] = (pick_named(bs_genes, PROTEIN_TARGETS['bsub']) +
                             pick_bsub_rna(bs_genes, bs_genomes))
            chromosomes['bsub'] = bs_genomes
            print('Barnyard B. subtilis RNA genes:',
                  [(x.gene, x.locus_tag, x.biotype) for x in plans['bsub']
                   if x.biotype != 'protein_coding'])
        make_sheet(args.output_dir, mode)
        write_design(args.output_dir, mode, plans, write_reads=True,
                     chromosomes=chromosomes)


if __name__ == '__main__':
    main()
