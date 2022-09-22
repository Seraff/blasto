#!/usr/bin/env python3

# Script 1

from Bio import SeqIO
from tqdm import tqdm

NA_PATH = "/Users/serafim/bio/blasto/data/new_assembly/results/annotation_nucleotides.fasta"
AA_PATH = "/Users/serafim/bio/blasto/data/new_assembly/results/annotation_peptides.fasta"
OUT_PATH = "/Users/serafim/bio/blasto/data/ribo_proteins_stats/scatter_plot/stop_share.csv"
GLU_TRP_CODONS = ['TAA', 'TAG', 'TGA', 'GAA', 'GAG', 'TGG']

def in_group_of(s, num):
    return [s[i:i+3] for i in range(0, len(s), num)]

def calc_stop_share(na_seq):
    codons = in_group_of(na_seq, 3)
    stats = dict([[e, 0] for e in GLU_TRP_CODONS])

    for codon in codons:
        if codon in GLU_TRP_CODONS:
            stats[codon] += 1

    stops_cnt = stats['TAA']+stats['TAG']+stats['TGA']
    all_cnt = stops_cnt+stats['GAA']+stats['GAG']+stats['TGG']

    if all_cnt != 0:
        share = stops_cnt/all_cnt
    else:
        print(f"Check the seq `{na_seq}`, no Glu/Trp codons found.")
        share = 0

    return share

def main():
    fasta_dict = {}

    for rec in SeqIO.parse(NA_PATH, 'fasta'):
        fasta_dict[rec.id] = {'na': rec.seq}
        fasta_dict[rec.id]['aa'] = rec.seq

    with open(OUT_PATH, 'w') as f:
        for header, data in tqdm(fasta_dict.items()):
            share = calc_stop_share(data['na'])
            f.write(f"{header},{share}\n")

if __name__ == "__main__":
    main()
