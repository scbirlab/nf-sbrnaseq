#!/usr/bin/env python3
"""Deterministic paired-end sbRNA-seq Illumina smoke fixture, stdlib only.

The gene fragments are from E. coli K-12 MG1655 BacMap reference sequences.
This is a synthetic library with simplified adapter locations designed to
exercise the pipeline's actual two-stage cutadapt + UMI-tools regex interface.
"""
from collections import Counter
from pathlib import Path
import csv
import gzip
import random

HERE = Path(__file__).resolve().parent
OUT = HERE / 'fastq'
INP = HERE / 'inputs'
SAMPLE = 'synthetic_illumina'
RNG = random.Random(20261016)

# 7-nt barcodes in three sequential rounds, 8-nt UMI.
BC = [
    ['ACGTGCA', 'TGCATGC', 'GATCCGA', 'CTAGTCT'],
    ['TTGACCA', 'CCATGGT', 'AGCTAGC', 'GGCATTA'],
    ['ATGCCGA', 'CGTATGC', 'TCGGAAT', 'GACTTCC'],
]
LEFT_R1 = 'ACGTTCGA'
LEFT_R2 = 'TGCAGTCA'
RIGHT_R1 = 'GTCAGTCA'
RIGHT_R2 = 'CTAGGCTA'

# DNA is oriented 5'->3' on the annotated gene's transcript strand.
# First ~400 bp of BacMap coding sequences; the same genes occur in ASM584v2.
GENES = {
'lacZ': (
'ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCG'
'CCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCA'
'GCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCT'
'GAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCC'
'CATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCT'
),
'lacI': (
'GTGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTTTCCCGCGTGGTGAACCAGGC'
'CAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCGGCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCAC'
'AACAACTGGCGGGCAAACAGTCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTC'
'GCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAACGAAGCGGCGTCGAAGCCTG'
),
'lacY': (
'ATGTACTATTTAAAAAACACAAACTTTTGGATGTTCGGTTTATTCTTTTTCTTTTACTTTTTTATCATGGGAGCCTACTT'
'CCCGTTTTTCCCGATTTGGCTACATGACATCAACCATATCAGCAAAAGTGATACGGGTATTATTTTTGCCGCTATTTCTC'
'TGTTCTCGCTATTATTCCAACCGCTGTTTGGTCTGCTTTCTGACAAACTCGGGCTGCGCAAATACCTGCTGTGGATTATT'
'ACCGGCATGTTAGTGATGTTTGCGCCGTTCTTTATTTTTATCTTCGGGCCACTGTTACAATACAACATTTTAGTAGGATC'
),
'rssA': (
'ATGAGAAAGATAAAAATAGGGCTGGCGCTGGGATCTGGCGCGGCGAGAGGTTGGTCGCATATTGGCGTTATTAATGCGCT'
'AAAAAAAGTGGGTATTGAAATTGATATCGTTGCAGGATGTTCAATTGGTTCGCTGGTGGGCGCTGCCTATGCATGCGATC'
'GATTATCTGCGCTGGAAGATTGGGTGACCTCTTTCAGTTATTGGGATGTTTTACGCCTGATGGATCTCTCCTGGCAGCGC'
'GGTGGGTTACTGCGCGGCGAGCGTGTCTTCAATCAATACCGCGAAATAATGCCGGAAACAGAGATCGAAAATTGTTCCCG'
),
'tehB': (
'ATGATCATTCGTGACGAAAACTATTTTACTGATAAATATGAATTAACCCGCACACACTCTGAAGTACTGGAAGCGGTGAA'
'AGTGGTTAAACCGGGTAAAACGCTGGATCTGGGCTGTGGCAATGGTCGTAACAGTCTTTACCTGGCAGCCAATGGTTATG'
'ATGTTGACGCATGGGATAAAAATGCCATGAGTATCGCCAACGTCGAGCGCATTAAATCCATTGAAAATCTGGATAATTTA'
'CACACCCGAGTCGTTGATCTGAATAACCTCACATTTGATAGACAGTACGATTTTATTCTTTCGACTGTGGTGCTGATGTT'
),
}


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def umi():
    return ''.join(RNG.choices('ACGT', k=8))


def write_fastq(f, name, seq):
    f.write(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    INP.mkdir(parents=True, exist_ok=True)
    for round_idx, choices in enumerate(BC, 1):
        assert len(choices) == len(set(choices))
        assert all(len(s) == 7 for s in choices)
        assert min(hamming(x, y) for i, x in enumerate(choices) for y in choices[i+1:]) >= 3
        with (INP / f'bc{round_idx}.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['Barcode ID', 'Sequence'])
            for i, seq in enumerate(choices, 1):
                w.writerow([f'bc{i:02d}', seq])

    # Regex columns match retained 5' adapters and then delete them using discard_ groups.
    patterns = (
        r'^(?P<discard_1>.{8})(?P<cell_1>.{7})(?P<discard_2>.{3}).*',
        r'^(?P<discard_3>.{8})(?P<cell_2>.{7})(?P<cell_3>.{7})(?P<umi_1>.{8})(?P<discard_4>.{2}).*',
    )
    with (HERE / 'sample-sheet.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['sample_id', 'fastq_pattern', 'genome_accession',
                    'bc1', 'bc2', 'bc3',
                    'adapter_read1_5prime', 'adapter_read2_5prime',
                    'adapter_read1_3prime', 'adapter_read2_3prime',
                    'umi_read1', 'umi_read2'])
        w.writerow([SAMPLE, SAMPLE, 'GCF_000005845.2',
                    'bc1.csv', 'bc2.csv', 'bc3.csv',
                    LEFT_R1, LEFT_R2, RIGHT_R1, RIGHT_R2,
                    patterns[0], patterns[1]])

    ids = 0
    expected = Counter()
    duplicate_reads = 0
    with gzip.open(OUT / f'{SAMPLE}_R1.fastq.gz', 'wt', compresslevel=6) as r1, \
         gzip.open(OUT / f'{SAMPLE}_R2.fastq.gz', 'wt', compresslevel=6) as r2:
        for i, b1 in enumerate(BC[0]):
            for j, b2 in enumerate(BC[1]):
                for k, b3 in enumerate(BC[2]):
                    cell = b1+b2+b3
                    for gidx, (gene, seq) in enumerate(GENES.items()):
                        for molecule in range(1 + ((i + 2*j + 3*k + 2*gidx) % 4)):
                            start = 18 + (gidx * 2) + molecule * 25
                            insert = seq[start:start + 160]
                            assert len(insert) == 160
                            read1 = LEFT_R1 + b1 + 'AGT' + insert[:70] + RIGHT_R1
                            read2 = LEFT_R2 + b2 + b3 + umi() + 'GG' + rc(insert[-70:]) + RIGHT_R2
                            # Each molecule has one or two identical PCR observations.
                            copies = 2 if (i+j+k+gidx+molecule)%3 == 0 else 1
                            for c in range(copies):
                                ids += 1
                                name = f'SYN{ids:09d}'
                                write_fastq(r1, name, read1)
                                write_fastq(r2, name, read2)
                                if c: duplicate_reads += 1
                            expected[(cell, gene)] += 1
    with (HERE/'design-counts.tsv').open('w', newline='') as f:
        w=csv.writer(f, delimiter='\t')
        w.writerow(['cell_barcode','gene','designed_molecules'])
        for (cell,gene), value in sorted(expected.items()):
            w.writerow([cell,gene,value])
    print(f'Wrote {ids} paired reads, {sum(expected.values())} designed molecules, '
          f'{duplicate_reads} PCR duplicate pairs, {len(BC[0])*len(BC[1])*len(BC[2])} barcodes')

if __name__ == '__main__':
    main()
