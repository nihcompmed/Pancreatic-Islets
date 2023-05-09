import numpy as np
import pickle
import os
import subprocess
import math

islets = pickle.load(open('far_unmatched_islets.p', 'rb'))


stochastic_dirr = 'stochastic'

n_pert = 50

n_perm = 50


max_pert = 0.1

for fprefix in islets:

    print(fprefix)

    islet_dirr = stochastic_dirr + '/' + fprefix

    pert_dirr = islet_dirr + '/Perturbations_ad'

    max_pert_dirr = pert_dirr + '/Perturbation'+str(max_pert)

    cluster_to_cycles = dict()

    for pert in range(n_pert):

        this_pert_dirr = max_pert_dirr + '/Perturbation' + str(pert)

        # This is of format loop:cluster
        pert_undead_file = this_pert_dirr + '/adundead_PH.csv'

        pert_undead = open(pert_undead_file, 'r')

        for line in pert_undead:
            line = line.strip('\n')
            line = line.split(':')
            if line[1] == '':
                continue

            this_cluster = line[1].split(',')[:-1]

            this_cluster = [int(x) for x in this_cluster]

            this_cluster = frozenset(this_cluster)

            this_loop = [int(x) for x in line[0].split(',')[:-1]]

            if this_cluster not in cluster_to_cycles:
                cluster_to_cycles[this_cluster] = this_loop
            else:
                if len(this_loop) < len(cluster_to_cycles[this_cluster][0]):
                    cluster_to_cycles[this_cluster] = [this_loop]
                else if len(this_loop) == len(cluster_to_cycles[this_cluster][0]):
                    cluster_to_cycles[this_cluster].append(this_loop)


        print(cluster_to_cycles)

        exit()



