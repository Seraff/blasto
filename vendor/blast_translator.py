#!/usr/bin/python
import sys
import csv
from Bio import SeqIO

def parse_csv(csv_path, delimiter=','):
    with open(csv_path) as handle_file:
        handle_csv = csv.reader(handle_file, delimiter=delimiter)
        results = []
        for row in handle_csv:
            results.append(row)
        handle_file.close()
    return results

def back_translate(frame, aa_start, aa_end, add_nt=0):
    if frame in (1, 2, 3):
        nt_start = (aa_start-1) * 3 + (frame-1)
        nt_end = aa_end * 3 + (frame-1) + add_nt
        return [nt_start, nt_end]
    elif frame in (4, 5, 6):
        nt_start = (aa_start-1) * 3 + (frame-4)
        nt_end = aa_end * 3 + (frame-4) + add_nt
        return [nt_start, nt_end]
    else:
        print "Error in frame"
        return 1

def cds_from_query(blast_hits, add_nt=0):
    results = []
    for bh in blast_hits:
        query_id = bh[0]
        if query_id == "qseqid":
            continue

        print query_id
        nt_id = query_id[:-2]
        frame = int(query_id[-1])
        qstart = int(bh[10])
        qend = int(bh[11])

        nt_indexes = back_translate(frame, qstart, qend, add_nt=add_nt)
        new_bh = bh
        new_bh[10] = str(nt_indexes[0])
        new_bh[11] = str(nt_indexes[1])
        results.append(','.join(new_bh))
    return results

blast_csv_path = 'data/hits_RNA_translated_best.csv'
outpath = 'data/hits_RNA_best.csv'

blast_hits = parse_csv(blast_csv_path)

outfile = open(outpath, 'w')
outfile.write("qseqid,qlen,sseqid,slen,length,evalue,pident,bitscore,mismatch,gaps,qstart,qend,sstart,send,End-send,alen/slen\n")
for result in cds_from_query(blast_hits):
    outfile.write(result + "\n")
outfile.close()
