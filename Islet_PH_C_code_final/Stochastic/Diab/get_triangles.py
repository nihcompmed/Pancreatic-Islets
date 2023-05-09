import numpy as np
import pickle
import os
import subprocess
import math



islets = pickle.load(open('far_unmatched_islets.p', 'rb'))

n_pert = 10

for fprefix in islets:

    print(fprefix)

    islet_dirr = stochastic_dirr + '/' + fprefix

    ad_vert_PH_file = '../VerticesEdges/' + fprefix + '_advertices_PH.csv'
    b_vert_PH_file = '../VerticesEdges/' + fprefix + '_bvertices_PH.csv'

    ad_edge_file = '../VerticesEdges/' + fprefix + '_adedges_PH.csv'
    b_edge_file = '../VerticesEdges/' + fprefix + '_bedges_PH.csv'

    ad_verts = np.loadtxt(ad_vert_PH_file, delimiter=',')
    b_verts = np.loadtxt(b_vert_PH_file, delimiter=',')

    ad_edges = np.loadtxt(ad_edge_file, delimiter=',')
    b_edges = np.loadtxt(b_edge_file, delimiter=',')

    #####################################################################
    ## ad-mantles
    #####################################################################

    if ad_edges.shape[0] < 3:
        continue

    pert_dirr = islet_dirr + '/Perturbations_ad'

    max_pert = 0.1
    max_pert_dirr = pert_dirr + '/Perturbation'+str(max_pert)

    # Perturb ad vertices
    n_pts = ad_verts.shape[0]

    for pert in range(n_pert):

        this_pert_dirr = max_pert_dirr + '/Perturbation' + str(pert)

        new_vert_file = this_pert_dirr + '/advertices_PH.csv'

        new_edge_file = this_pert_dirr + '/adedges_PH.csv'


        # Get list of triangles
        subprocess.run(["./stochastic_process.o"\
                        , this_pert_dirr
                        ])




