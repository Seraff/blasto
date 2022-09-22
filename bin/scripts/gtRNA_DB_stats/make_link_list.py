#!/usr/bin/env python3

from pprint import pprint

import requests
from bs4 import BeautifulSoup

BASE_URL = 'http://gtrnadb.ucsc.edu'
LIST_URL = f'{BASE_URL}/cgi-bin/trna_chooseorg?org=&genelist=+++Go+to+Gene+List++'
OUT_PATH = 'output/url_list.txt'

page = requests.get(LIST_URL)

soup = BeautifulSoup(page.content, 'html.parser')

links = []
for li in soup.find_all('li'):
  link = li.a.get('href').strip('.')
  link = f'{BASE_URL}{link}'
  link = link.replace('-gene-list.html', '-align.html')
  links.append(link)

with open(OUT_PATH, 'w') as f:
  for link in links:
    f.write(f'{link}\n')

# call wget -i ../url_list.txt afterwards
