import pickle
import os
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
import subprocess



# THIS ASSUMES THAT stochastic_far_unmatched.py has been executed

islets = pickle.load(open('far_unmatched_islets.p', 'rb'))


dirr = 'stochastic'

thresh_dirr = 'Thresholds'
data_dirr = 'DevData/all_islets'

n_pert = 10


fcount = 0

for fprefix in islets:

    fcount += 1

    islet_dirr = dirr + '/' + fprefix

    pert_dirr = islet_dirr + '/Perturbations'


    for pert in range(n_pert):
        print(fprefix, pert)

        this_pert_dirr = pert_dirr + '/Perturbation'+str(pert)

        locs_file = this_pert_dirr + '/locs.tsv'

        subprocess.run(["./SimplexCode_12519_v3/Simplex_VersionForPaper"\
                        , locs_file\
                        , "321Area"\
                        , "B"\
                        , "-grdir"\
                        , this_pert_dirr
                        , "-fout"
                        , "locs"
                        , "-runGr"
                        , "y"
                        ])



