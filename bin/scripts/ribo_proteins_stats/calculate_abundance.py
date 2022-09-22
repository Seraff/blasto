#!/usr/bin/env python3

# Script 2
# makeblastdb -in ../../new_assembly/results/annotation_peptides.fasta -parse_seqids -dbtype prot -out annot_aa_db
# blastp -query ../../proteomics/peptides.fasta -db annot_aa_blast_db/annot_aa_db -outfmt '10 qseqid,sseqid,evalue,pident,bitscore' -num_threads 2 -out annot_aa_vs_peptides_blast_results.csv
# awk -F ',' 'BEGIN {OFS=","} { if ($3 >= 100) print }' annot_aa_vs_peptides_blast_results.csv > annot_aa_vs_peptides_blast_results_filtered.csv

from Bio import SeqIO
from tqdm import tqdm

BLAST_PATH = "/Users/serafim/bio/blasto/data/ribo_proteins_stats/scatter_plot/annot_aa_vs_peptides_blast_results_filtered.csv"
STOP_SHARE_PATH = "/Users/serafim/bio/blasto/data/ribo_proteins_stats/scatter_plot/stop_share.csv"
OUT_PATH = "/Users/serafim/bio/blasto/data/ribo_proteins_stats/scatter_plot/scatter_plot.csv"

def main():
    hits_stats = {}

    with open(BLAST_PATH) as f:
        for line in f:
            gene_id = line.split(',')[1].strip()

            if gene_id not in hits_stats:
                hits_stats[gene_id] = 0

            hits_stats[gene_id] += 1

    max_abundance = max(hits_stats.values())

    with open(OUT_PATH, 'w') as out_f:
        with open(STOP_SHARE_PATH) as f:
            for line in f:
                gene_id = line.split(',')[0].strip()

                if gene_id not in hits_stats:
                    abundance = 0
                else:
                    abundance = hits_stats[gene_id]

                abundance = abundance/max_abundance

                line = f"{line.strip()},{abundance}\n"

                out_f.write(line)


if __name__ == "__main__":
    main()
