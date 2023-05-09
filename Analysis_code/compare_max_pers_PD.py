import matplotlib.pyplot as plt
import numpy as np
import persim
import os
import pickle
import islet_helper as ih
import itertools as it
from joblib import Parallel, delayed
import math
import matplotlib as mpl

info_file = 'all_islets_info.csv'

PD_info_file = 'all_PD_info.p'

PD_info = pickle.load(open(PD_info_file, 'rb'))



typ = 'ad'

fname = 'figures/max_pers_' + typ + '_cells.pdf'

ptiles_dict_adPH = ih.plot_max_pers_distri_dev(PD_info, fname, typ)


typ = 'b'

fname = 'figures/max_pers_' + typ + '_cells.pdf'

ptiles_dict_bPH = ih.plot_max_pers_distri_dev(PD_info, fname, typ)

xx = []

ptile50_adPH = []
ptile95_adPH = []

ptile50_bPH = []
ptile95_bPH = []


for stage in range(4):

    this_ptiles = ptiles_dict_adPH[stage]
    ptile50_adPH.append(this_ptiles[1])
    ptile95_adPH.append(this_ptiles[2])

    this_ptiles = ptiles_dict_bPH[stage]
    ptile50_bPH.append(this_ptiles[1])
    ptile95_bPH.append(this_ptiles[2])
    
    xx.append(stage)


plt.cla()
plt.clf()
plt.close()
mpl.rcParams.update(mpl.rcParamsDefault)

plt.plot(xx, ptile50_adPH, color='blue', marker='o', ls='--')
plt.plot(xx, ptile95_adPH, color='blue', marker='o')


plt.plot(xx, ptile50_bPH, color='red', marker='o', ls='--')
plt.plot(xx, ptile95_bPH, marker='o', color='red')

plt.xticks(ticks=[0,1,2,3], labels=['Stage 0', 'Stage 1', 'Stage 2', 'Stage 3'], fontsize=16)

plt.ylabel('Maximum dimension-1 persistence', fontsize=16)

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='red', lw=2),
                ]
legend1 = plt.legend(custom_lines, [r'$\alpha\delta$ topology',
                                r'$\beta$ topology',
                                ],
                                loc=(0.65, 0.85))

custom_lines2 = [
                Line2D([0], [0], color='black', ls='-'),
                Line2D([0], [0], color='black', ls='--'),
                ]
legend2 = plt.legend(custom_lines2, [
                                r'$95\%$-tile',
                                r'$50\%$-tile'],
                            loc=(0.65, 0.72))

plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)

plt.savefig('max_pers_across_stages.pdf', dpi=600)


