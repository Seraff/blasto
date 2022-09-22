#!/usr/bin/env python3

import argparse
from pprint import pprint
from Bio import SeqIO
from tqdm import tqdm

# FQ1_FILENAME = 'test_1.fq'
# FQ2_FILENAME = 'test_2.fq'
# OUT1_FILENAME = 'test_1_ANY_TA_CLEANED.fq'
# OUT2_FILENAME = 'test_2_ANY_TA_CLEANED.fq'
FQ1_FILENAME = 'p57_3-end_trimmed_1.fq'
FQ2_FILENAME = 'p57_3-end_trimmed_2.fq'
OUT1_FILENAME = 'p57_3-end_trimmed_1_ANY_TA_CLEANED.fq'
OUT2_FILENAME = 'p57_3-end_trimmed_2_ANY_TA_CLEANED.fq'
POLY_MIN_CNT = 3


def has_poly_t(rec):
    return str(rec.seq).startswith('T'*POLY_MIN_CNT)


def has_poly_a(rec):
    return str(rec.seq).endswith('A'*POLY_MIN_CNT)


def has_any_poly(rec):
    return has_poly_t(rec) or has_poly_a(rec)


def has_any_TA_at_the_ends(rec):
    seq = str(rec.seq)
    return seq.startswith('T') or seq.startswith('A') or \
        seq.endswith('T') or seq.endswith('A')


def clean_from_poly(rec):
    if has_poly_t(rec):
        stripped = rec.seq.lstrip('T')
        cnt = len(rec.seq) - len(stripped)
        new_qual_list = rec.letter_annotations['phred_quality'][cnt:]
        rec.letter_annotations = {}
        rec.seq = stripped
        rec.letter_annotations['phred_quality'] = new_qual_list

    if has_poly_a(rec):
        stripped = rec.seq.rstrip('A')
        cnt = len(rec.seq) - len(stripped)
        new_qual_list = rec.letter_annotations['phred_quality'][:-cnt]
        rec.letter_annotations = {}
        rec.seq = stripped
        rec.letter_annotations['phred_quality'] = new_qual_list

    return rec


def clean_from_any_TA(rec):
    stripped = rec.seq.lstrip('TA')
    cnt = len(rec.seq) - len(stripped)
    if cnt > 0:
        new_qual_list = rec.letter_annotations['phred_quality'][cnt:]
        rec.letter_annotations = {}
        rec.seq = stripped
        rec.letter_annotations['phred_quality'] = new_qual_list

    stripped = rec.seq.rstrip('TA')
    cnt = len(rec.seq) - len(stripped)
    if cnt > 0:
        new_qual_list = rec.letter_annotations['phred_quality'][:-cnt]
        rec.letter_annotations = {}
        rec.seq = stripped
        rec.letter_annotations['phred_quality'] = new_qual_list

    return rec


def main():
    to_save = []

    f_left = SeqIO.parse(FQ1_FILENAME, 'fastq')
    f_right = SeqIO.parse(FQ2_FILENAME, 'fastq')


    with open(OUT1_FILENAME, 'w') as out_f_1, open(OUT2_FILENAME, 'w') as out_f_2:
        for rec_left, rec_right in tqdm(zip(f_left, f_right)):
            if str(rec_left.seq.lstrip('TA')) == '' or str(rec_right.seq.lstrip('TA')) == '':
                continue

            clean_from_any_TA(rec_left)
            clean_from_any_TA(rec_right)

            SeqIO.write(rec_left, out_f_1, 'fastq')
            SeqIO.write(rec_right, out_f_2, 'fastq')


if __name__ == '__main__':
    main()
