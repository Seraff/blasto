#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
from pprint import pprint

'''
All triplet distribution inside all hits
'''

parser = ArgumentParser(description = 'Gets codon statistics from fasta file and csv BLAST results')
parser.add_argument('-f', '--fasta', required = True, help = 'path to FASTA file with predicted genes')
parser.add_argument('-o', '--output', required = True, help = 'output path (.csv)')
options = parser.parse_args()

def group_by_nucleotides(rec):
    if len(rec.seq) % 3 != 0:
        print(f"WARNING: sequence {rec.id} is not a multiple of two")
    return [''.join(e).upper() for e in list(zip(*(iter(rec.seq),) * 3))]

STOPS = ('TAA', 'TAG', 'TGA')

stats = {}

for rec in tqdm(SeqIO.parse(options.fasta, 'fasta')):
    triplets = group_by_nucleotides(rec)

    seq_len = len(triplets)

    for i, triplet in enumerate(triplets):
        if triplet in STOPS:
            if triplet not in stats:
                stats[triplet] = []

            value = (i+1)/seq_len
            stats[triplet].append(value)

with open(options.output, 'w') as out_f:
    for triplet, values in stats.items():
        for value in values:
            out_f.write(f"{triplet},{value}\n")

