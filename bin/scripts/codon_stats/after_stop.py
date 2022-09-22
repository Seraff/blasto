#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

data = pd.read_csv('stop_per_triplet_stats.csv')
data = data.sort_values(by=['cnt'], ascending=False)
# data = data.iloc[0:45]

STOP_CODONS=('TAA', 'TGA', 'TAG')
colors = ['#FF0078' if e in STOP_CODONS else 'lightgrey' for e in data['codon']]
GAP = {'from': 400, 'to': 1750}
Y_LABEL_SIZE = 8
X_LABEL_SIZE = 8
PLOT_HEIGHT_RATIOS = [1, 2]

f, (ax_top, ax_bottom) = plt.subplots(ncols=1,
                                      nrows=2,
                                      sharex=True,
                                      gridspec_kw={'hspace': 0.15,
                                                   'height_ratios': PLOT_HEIGHT_RATIOS},
                                      figsize=(10,3))


## TOP
sns.barplot(x="codon",
            y="cnt",
            data=data,
            palette=colors,
            linewidth=0.5,
            edgecolor="black",
            ax=ax_top)
ax_top.set_ylim(GAP['to'],1950)
ax_top.tick_params(axis='y', labelsize=Y_LABEL_SIZE)
ax_top.tick_params(axis='x', width=0, grid_alpha=0)
sns.despine(ax=ax_top)
ax_top.spines['bottom'].set_visible(False)
ax_top.set(xlabel='', ylabel='')

## BOTTOM
sns.despine(ax=ax_bottom)
sns.barplot(x="codon",
            y="cnt",
            data=data,
            palette=colors,
            linewidth=0.5,
            edgecolor="black",
            ax=ax_bottom)

ax_bottom.set_ylim(0, GAP['from'])
ax_bottom.set_xlim(-1, len(data)-0.5)

ax_bottom.tick_params(axis='x',
                      length=0,
                      width=0,
                      labelrotation=90,
                      labelsize=X_LABEL_SIZE)

ax_bottom.tick_params(axis='y', labelsize=Y_LABEL_SIZE)
ax_bottom.set(xlabel='Codon', ylabel='    Count')

plt.subplots_adjust(bottom=0.35) # empty space

## DRAWING GAP SYMBOLS
d = .007

diag_mult = 7.5
top_diag_mult = diag_mult * PLOT_HEIGHT_RATIOS[1]
bottom_diag_mult = diag_mult * PLOT_HEIGHT_RATIOS[0]

kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False, linewidth=0.8)
ax_top.plot([-d, +d], [-d*top_diag_mult, +d*top_diag_mult], **kwargs) # top diagonal
kwargs = dict(transform=ax_bottom.transAxes, color='k', clip_on=False, linewidth=0.8)
ax_bottom.plot([-d, +d], [1-d*bottom_diag_mult, 1+d*bottom_diag_mult], **kwargs) # bottom diagonal

shift = 0.016
kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False, linewidth=0.8)
ax_top.plot([-d+shift, +d+shift], [-d*top_diag_mult, +d*top_diag_mult], **kwargs) # top diagonal
kwargs = dict(transform=ax_bottom.transAxes, color='k', clip_on=False, linewidth=0.8)
ax_bottom.plot([-d+shift, +d+shift], [1-d*bottom_diag_mult, 1+d*bottom_diag_mult], **kwargs) # bottom diagonal

## DRAWING ... AFTER THE PLOT
# kwargs = dict(transform=ax_bottom.transAxes,
#               color='k',
#               clip_on=False,
#               linewidth=0.8,
#               linestyle='dashed')

# ax_bottom.plot([1, 1+0.05], [0, 0], **kwargs) # dashed line

f.savefig("output.png", dpi=300)
# plt.subplot_tool()
# plt.show()
