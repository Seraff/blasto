#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

parser = ArgumentParser(description = 'Draws histogram')
parser.add_argument('-i', '--stats', required = True, help = 'path to CSV file with stats')
options = parser.parse_args()

CODON = 'TAA'

data = pd.read_csv(options.stats)

data = data.loc[data['codon'] == CODON]
vals = data[data.columns[1]]

sns.displot(vals, bins=50).set(title=f'{CODON} distribution')
plt.subplots_adjust(top=0.9)
plt.show()
