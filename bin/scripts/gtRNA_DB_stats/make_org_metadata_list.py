#!/usr/bin/env python3

from pprint import pprint
from glob import glob
import subprocess
from tqdm import tqdm

def get_org_name(html_path):
    cmd = f"grep -A10 'page-header' {html_path} | grep 'h5'"
    line = subprocess.check_output(cmd, shell=True)
    line = line.decode("utf-8").strip()
    line = line[4:-5]
    return line

def process_organism():
    pass

URL_LIST_PATH = 'output/url_list.txt'

def main():
    url_list = {}
    with open(URL_LIST_PATH) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue

            html_filename = line.split('/')[-1]
            url_list[html_filename] = line

    result = []

    htmls = glob('output/html/*.html')
    for html_path in tqdm(htmls):
        org_name = get_org_name(html_path)
        html_filename = html_path.split('/')[-1]

        url = url_list[html_filename]

        result.append(f"{org_name},{url}")

    with open('output/org_name_urls.csv', 'w') as out_f:
        out_f.write('\n'.join(result) + '\n')



if __name__ == '__main__':
    main()
