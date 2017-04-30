#!/usr/bin/python
import sys
import os
import math
import time
import re
import ntpath
from Bio import SeqIO
import matplotlib.pyplot as plt
import progressbar

SL = 'AACGCATTTTTTGTTACAGTTTCTGTACTTTATTG'
SL_HAMM_THRESHOLD = 0.05
SL_THRESHOLD = 7

POLY_T = 'T'*10
POLY_T_PATTERN = re.compile("\AT+")
POLY_A_THRESHOLD = 5

def hamming_distance(s1, s2):
  if len(s1) != len(s2):
    raise ValueError("Undefined for sequences of unequal length")
  return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def begins_with(read, substring, threshold):
  substr_len = len(substring)
  h_dist = hamming_distance(read[0:substr_len], substring)
  err_percentage = h_dist/float(substr_len)
  # return read.find(substring) == 0 # non hamming
  return err_percentage <= threshold

def max_intersection(seq, reference, mi, hamm_threshold):
  target_len = len(seq) if len(seq) < len(reference) else len(reference)

  for i in range(SL_THRESHOLD, target_len):
    tail = reference[len(reference)-i:]
    if begins_with(seq, tail, hamm_threshold):
      return i

  return 0

def process_sl(record, stats, outfile, full_outfile):
  na_count = max_intersection(str(record.seq), SL, SL_THRESHOLD, SL_HAMM_THRESHOLD)
  rev_record = deep_reverse_complement(record)
  na_count_rev = max_intersection(str(rev_record.seq), SL, SL_THRESHOLD, SL_HAMM_THRESHOLD)

  for i, cnt in enumerate([na_count, na_count_rev]):
    if cnt >= SL_THRESHOLD:
      rec = record if i == 0 else rev_record
      add_to_stats(cnt, rec[cnt:], stats, outfile)
      SeqIO.write(record, full_outfile, "fastq")

def process_poly_a(record, stats, outfile, full_outfile):
  seq = str(record.seq)
  rev_record = deep_reverse_complement(record)
  rev_seq = str(rev_record.seq)

  for i, s in enumerate([seq, rev_seq]):
    result = POLY_T_PATTERN.search(s)
    if result:
      cnt = len(result.group(0))
      if cnt >= POLY_A_THRESHOLD:
        rec = rev_record[:-cnt] if i == 0 else record[:-cnt]
        add_to_stats(cnt, rec, stats, outfile)
        SeqIO.write(record, full_outfile, "fastq")

def add_to_stats(cnt, rec, stats, outfile):
  if not cnt in stats:
    stats[cnt] = 0
  stats[cnt] += 1

  SeqIO.write(rec, outfile, "fastq")

def deep_reverse_complement(rec):
  new_rec = rec.reverse_complement()
  new_rec.id = rec.id
  new_rec.name = rec.name
  new_rec.description = rec.description
  return new_rec

reads_paths = sys.argv[1:]

for path in reads_paths:
  f_size = int(os.popen('wc -l %s' % path).read().split()[0])
  # bar = progressbar.ProgressBar(max_value=int(math.ceil(f_size/4.0)))

  sl_stats = {}
  poly_t_stats = {}

  outfile_name = ntpath.basename(path).split('.')[0]
  sl_outfile = open("results/%s_SL.fq" % outfile_name, "w")
  sl_full_outfile = open("results/%s_SL_full.fq" % outfile_name, "w")
  pa_outfile = open("results/%s_poly_a.fq" % outfile_name, "w")
  pa_full_outfile = open("results/%s_poly_a_full.fq" % outfile_name, "w")

  count = 0
  for i, record in enumerate(SeqIO.parse(path, "fastq")):
    process_sl(record, sl_stats, sl_outfile, sl_full_outfile)
    process_poly_a(record, poly_t_stats, pa_outfile, pa_full_outfile)

    count += 1
    if count % 100000 == 0:
      print count
    # bar.update(i+1)

  print sl_stats
  # plt.plot(sl_stats.keys(), sl_stats.values(), color='blue')
  # for x, y in sl_stats.iteritems():
  #   plt.annotate(y, xy=(x, y), xytext=(x, y+1000), fontsize=8, color='blue')

  print poly_t_stats
  # plt.plot(poly_t_stats.keys(), poly_t_stats.values(), color='red')
  # for x, y in poly_t_stats.iteritems():
  #   plt.annotate(y, xy=(x, y), xytext=(x, y+5000), fontsize=8, color='red')

  # plt.axis([0, 35, 0, 100000])
  # plt.show()
