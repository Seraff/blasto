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


org_urls = {}
with open('output/org_name_urls.csv') as f:
  for line in f:
    line = line.strip()
    if line == '':
      continue

    splitted = line.split(',')
    org_urls[splitted[0]] = splitted[1]

with open('output/anticodon_stats.csv', 'w') as f:
  for org_name, aa_data in data.items():
    url = org_urls[org_name]

    for aa_name, data in aa_data.items():
      for anticodon_data in data['anticodons']:
        if 'name' not in anticodon_data:
          continue

        aa_name = aa_name.strip().split(' ')[0]
        name = anticodon_data['name']
        anticodon_len = len(anticodon_data['anticodon'])
        score = anticodon_data['score']

        f.write(f"{url},{org_name},{name},{aa_name},{anticodon_len},{score}\n")
