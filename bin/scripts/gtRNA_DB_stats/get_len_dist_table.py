#!/usr/bin/env python3

from pprint import pprint

AM_ACIDS = ['Ala tRNAs', 'Arg tRNAs', 'Asn tRNAs', 'Asp tRNAs',
            'Cys tRNAs', 'Gln tRNAs', 'Gly tRNAs', 'His tRNAs',
            'Ile tRNAs', 'Ile2 tRNAs', 'Leu tRNAs', 'Lys tRNAs',
            'Met tRNAs', 'fMet tRNAs', 'Phe tRNAs', 'Pro tRNAs',
            'Ser tRNAs', 'Thr tRNAs', 'Trp tRNAs', 'Tyr tRNAs',
            'Val tRNAs', 'Glu tRNAs', 'iMet tRNAs', 'SeC tRNAs',
            'Sup tRNAs']

anticodon_len_stats = 'output/anticodon_lengths.csv'
output_path = 'output/len_dist_table.csv'


def main():
    stats = {}

    with open(anticodon_len_stats) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue

            splitted = line.split(',')

            org_name = splitted[0]
            aa_name = splitted[1]
            anticodon_len = int(splitted[2])

            if aa_name not in stats:
                stats[aa_name] = {}

            if anticodon_len not in stats[aa_name]:
                stats[aa_name][anticodon_len] = 0

            stats[aa_name][anticodon_len] += 1

    lengths = list(range(1, 7))

    with open(output_path, 'w') as f:
        f.write(f"AA_NAME,{','.join([str(e) for e in lengths])}\n")
        for aa_name, anticodon_len_data in stats.items():
            length_stats = []

            for i in lengths:
                cnt = 0
                if i in anticodon_len_data:
                    cnt = anticodon_len_data[i]

                length_stats.append(cnt)

            line = f"{aa_name},{','.join([str(e) for e in length_stats])}\n"
            f.write(line)


    for aa in AM_ACIDS:
        pass
    pass


if __name__ == '__main__':
    main()
