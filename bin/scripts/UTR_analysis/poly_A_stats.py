#!/usr/bin/env python3

import json
import re
from pprint import pprint

OUTPUT_FOLDER = 'output_06_07'
UTR_STATS_PATH = f'{OUTPUT_FOLDER}/utr_statistics.json'
POLY_A_REGEXP = re.compile('.+?(A+)$')
OUTPUT_POLY_A_STATS_FILENAME = f'{OUTPUT_FOLDER}/poly_a_stats.json'
OUTPUT_POLY_A_FASTA_FILENAME = f'{OUTPUT_FOLDER}/trans_UTRs_with_polyass.fasta'

def get_ploy_A_cnt(seq):
    match = re.match(POLY_A_REGEXP, seq)
    if match:
        return len(match.group(1))
    else:
        return 0

def main():
    result = {}

    with open(UTR_STATS_PATH) as f:
        stats = json.load(f)

    for group in ('NOT_EQUAL', 'EQUAL'):
        for rec in stats[group]["records"]:
            trn_seq = rec['trn_utr']
            gen_seq = rec['gen_utr']

            trn_cnt = get_ploy_A_cnt(trn_seq)
            gen_cnt = get_ploy_A_cnt(gen_seq)
            cnt = trn_cnt - gen_cnt

            if cnt not in result:
                result[cnt] = []

            data = {'id': rec['transcript_seq_id'],
                    'trn_utr': rec['trn_utr'],
                    'gen_utr': rec['gen_utr']}

            result[cnt].append(data)

    result = [[k, v] for k, v in result.items()]
    result = sorted(result, key=lambda x: x[0])
    result = dict(result)

    with open(OUTPUT_POLY_A_STATS_FILENAME, 'w') as f:
        json.dump(result, f, indent=4)

    s = [[k, len(v)] for k, v in result.items()]
    s = sorted(s, key=lambda x: x[0])
    s = dict(s)
    pprint(s)

    with open(OUTPUT_POLY_A_FASTA_FILENAME, 'w') as out_f:
        for cnt, data in result.items():
            if cnt <= 0:
                continue

            for entry in data:
                header = f"{entry['id']}_{cnt}_poly_A"
                seq = entry['trn_utr']
                out_f.write(header + '\n')
                out_f.write(seq + '\n')




if __name__ == '__main__':
    main()
