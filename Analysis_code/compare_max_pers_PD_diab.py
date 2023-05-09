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


info_file = 'all_islets_info_diabetic.csv'

PD_info_file = 'all_PD_info_diab.p'

PD_info = pickle.load(open(PD_info_file, 'rb'))


typ = 'ad'

fname = 'max_pers_' + typ + '_cells_diab.pdf'

ptiles_dict_adPH = ih.plot_max_pers_distri_diab(PD_info, fname, typ)

typ = 'b'

fname = 'max_pers_' + typ + '_cells_diab.pdf'

ptiles_dict_bPH = ih.plot_max_pers_distri_diab(PD_info, fname, typ)

cats = ['C', 'D']

xx = []

ptile50_adPH = []
ptile95_adPH = []

ptile50_bPH = []
ptile95_bPH = []


for idx, cat in enumerate(cats):

    this_ptiles = ptiles_dict_adPH[cat]
    ptile50_adPH.append(this_ptiles[1])
    ptile95_adPH.append(this_ptiles[2])

    this_ptiles = ptiles_dict_bPH[cat]
    ptile50_bPH.append(this_ptiles[1])
    ptile95_bPH.append(this_ptiles[2])
    
    xx.append(idx)


plt.cla()
plt.clf()
plt.close()
mpl.rcParams.update(mpl.rcParamsDefault)

plt.plot(xx, ptile50_adPH, color='blue', marker='o', ls='--', lw=2)
plt.plot(xx, ptile95_adPH, color='blue', marker='o', lw=2)


plt.plot(xx, ptile50_bPH, color='red', marker='o', ls='--', lw=2)
plt.plot(xx, ptile95_bPH, marker='o', color='red', lw=2)


plt.xticks(ticks=[0,1], labels=['Control', 'Diabetic'], fontsize=20)

plt.ylabel('Maximum dimension-1 persistence', fontsize=20)

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='red', lw=2),
                ]
legend1 = plt.legend(custom_lines, [r'$\alpha\delta$ topology',
                                r'$\beta$ topology',
                                ],
                                loc=(0.8, 0.6))

custom_lines2 = [
                Line2D([0], [0], color='black', ls='-'),
                Line2D([0], [0], color='black', ls='--'),
                ]
legend2 = plt.legend(custom_lines2, [
                                r'$95\%$-tile',
                                r'$50\%$-tile'],
                            loc=(0.8, 0.5))

plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)

plt.savefig('max_pers_across_stages_diab.pdf', dpi=600)










