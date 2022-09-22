#!/usr/bin/env python3

from pprint import pprint
import pdb
import json
import re
import subprocess

import tqdm

DATA_PATH = 'output/data.json'
SEC_STR_REGEXP = re.compile(r'^[><\.]+$')
DOTS_REGEXP = re.compile(r'\.+')

def get_anticodon_real_seq(full_string, start, end):
    before = full_string[start-2:start]
    after = full_string[end+1:end+3]

    if before == '>.':
        start -= 2

        while full_string[start] == '>':
            start -= 1
        start += 1

    if after == '.>':
        end += 2

        while full_string[end] == '>':
            end += 1
        end -= 1

    return full_string[start:end+1]

def get_org_url(org_name):
  output = subprocess.check_output(f'grep -rl \'{org_name}\' output/html', shell=True)
  return output.split('\n')[0]

with open(DATA_PATH) as f:
  data = f.read()

data = json.loads(data)

# { '<org name>': { 1: 0, 2: 2, 3: 42 } }
stats = {}

for org_name, org_data in tqdm.tqdm(data.items()):

  for aa_name, aa_data in org_data.items():
    secs = [e for e in aa_data['lines'] if re.match(SEC_STR_REGEXP, e)]

    for i, anticodon in enumerate(aa_data['anticodons']):
      current_seq = secs[i]
      anticodon_seq = get_anticodon_real_seq(current_seq, anticodon['start'], anticodon['end'])
      anticodon_len = len(anticodon_seq)

      if org_name not in stats:
        stats[org_name] = {}
      if aa_name not in stats[org_name]:
        stats[org_name][aa_name] = {}

      if anticodon_len not in stats[org_name][aa_name]:
        stats[org_name][aa_name][anticodon_len] = 0

      stats[org_name][aa_name][anticodon_len] += 1

org_urls = {}
with open('output/org_name_urls.csv') as f:
  for line in f:
    line = line.strip()
    if line == '':
      continue

    splitted = line.split(',')
    org_urls[splitted[0]] = splitted[1]

with open('output/anticodon_lengths.csv', 'w') as f:
  for org_name, data in stats.items():
    # url = org_urls[org_name]
    for aa_name, codon_len_data in data.items():
      for anticodon_len, cnt in codon_len_data.items():
        f.write(f'{org_name},{aa_name},{anticodon_len},{cnt}\n')
