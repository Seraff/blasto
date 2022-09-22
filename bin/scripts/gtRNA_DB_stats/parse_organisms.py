#!/usr/bin/env python3
# HTML -> JSON

from pprint import pprint
import pdb
import json
import glob

import requests
from bs4 import BeautifulSoup
from bs4.element import Tag
from bs4.element import NavigableString
import tqdm
import re

HTML_PATH = 'output/html'
# HTML_PATH = 'output/html_test'
OUT_PATH = 'output/data.json'
SEC_STR_REGEXP = re.compile(r'^[><\.]+$')
ANTICODON_ID = 'b3'

def is_float(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def process_link(file):
  result = {}

  with open(file) as f:
    page = f.read()

  soup = BeautifulSoup(page, 'html.parser')
  org_name = soup.find('header').h5.get_text()

  result[org_name] = {}

  for aln in soup.find_all(attrs = { 'class': 'seq_alignment' }):
    aa_name = aln.find_parent().find_previous('a').get_text().strip()

    if aa_name == '' or not aa_name:
      print(f'Unable to parse amino-acid in file {link}')
      exit(1)

    if aa_name in result[org_name]:
      print(f'Unexpected error!')
      exit(1)

    content = aln.get_text().strip().split('\n')
    if aa_name not in result[org_name]:
      result[org_name][aa_name] = {}

    result[org_name][aa_name]['lines'] = content
    result[org_name][aa_name]['anticodons'] = []

    # finding anticodon arms
    for line in str(aln).split('\n'):
      soupized = BeautifulSoup(line, 'html.parser')
      line_str = soupized.get_text()

      if re.match(SEC_STR_REGEXP, line_str):
        # line which looks like >>>>>>>..>>>>........<<<<.>>>>>....

        cnt = 0
        anticodon_start = None
        anticodon_end = None

        for el in soupized:
          if type(el) == Tag:
            if el.attrs['id'] == ANTICODON_ID and el.get_text() == '>':
              if anticodon_start == None:
                anticodon_start = cnt

            else:
              if anticodon_start != None and anticodon_end == None:
                anticodon_end = cnt - 1

            cnt += 1

          elif type(el) == NavigableString:
            if anticodon_start != None and anticodon_end == None:
              anticodon_end = cnt - 1

            cnt += len(el)

          else:
            raise(Exception(f"Strange element: {el} in {aa_name}, {org_name}"))

        anticodon = ''.join([e.get_text() for e in soupized.find_all('span', id='b3')]).replace('<', '')
        data = {'anticodon': anticodon,
                'start': anticodon_start,
                'end': anticodon_end}

        if anticodon.strip() != '':
          result[org_name][aa_name]['anticodons'].append(data)

      elif (len(line_str.strip()) > 0):
        a_codons_list = result[org_name][aa_name]['anticodons']

        if len(a_codons_list) > 0:
          last_data = a_codons_list[-1]

          # check score line
          splitted = line_str.strip().split(' ')

          if (len(splitted) == 4) and ('name' not in last_data):
            name = splitted[1]
            score = float(splitted[3])

            last_data['name'] = name
            last_data['score'] = score


  return result


final_result = {}

links = glob.glob(HTML_PATH + "/*.html")

# cnt = 0
for link in tqdm.tqdm(links):
  parsed_data = process_link(link)
  final_result.update(parsed_data)

  # cnt += 1
  # if cnt >= 5:
  #   break

with open(OUT_PATH, 'w') as fp:
  json.dump(final_result, fp, indent=4)

# pprint(final_result)
print('Done')
