import pickle
import os
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
import subprocess


def get_edge_length(edge, locs):

    x1, y1 = locs[int(edge[0])]
    x2, y2 = locs[int(edge[1])]

    return math.sqrt((x1-x2)**2 + (y1-y2)**2)


# THIS ASSUMES THAT stochastic_far_unmatched.py has been executed

islets = pickle.load(open('far_unmatched_islets.p', 'rb'))


is_sorted = lambda a: np.all(a[:-1] <= a[1:])

dirr = 'stochastic'

thresh_dirr = 'Thresholds'
data_dirr = 'DevData/all_islets'

n_pert = 10

n_perm = 10

fcount = 0

for fprefix in islets:

    fcount += 1

    islet_dirr = dirr + '/' + fprefix

    pert_dirr = islet_dirr + '/Perturbations'

    print(fprefix)

    for pert in range(n_pert):

        print(pert, end='\r')

        this_pert_dirr = pert_dirr + '/Perturbation'+str(pert)

        ad_vert_PH_file = this_pert_dirr + '/advertices_PH.csv'
        b_vert_PH_file = this_pert_dirr + '/bvertices_PH.csv'

        ad_vertices = np.loadtxt(ad_vert_PH_file, delimiter=',')
        b_vertices = np.loadtxt(b_vert_PH_file, delimiter=',')

        ad_edge_file = this_pert_dirr + '/adedges_PH.csv'
        b_edge_file = this_pert_dirr + '/bedges_PH.csv'

        ad_edges = np.loadtxt(ad_edge_file, delimiter=',')
        b_edges = np.loadtxt(b_edge_file, delimiter=',')

        ad_triangles_file = this_pert_dirr + '/adtriangles_PH.csv'
        b_triangles_file = this_pert_dirr + '/btriangles_PH.csv'

        ad_triangles = np.loadtxt(ad_triangles_file, delimiter=',')
        b_triangles = np.loadtxt(b_triangles_file, delimiter=',')


        if len(ad_edges) and ad_edges.ndim == 1:
            ad_edges = np.reshape(ad_edges, (1,2))
        if len(b_edges) and b_edges.ndim == 1:
            b_edges = np.reshape(b_edges, (1,2))

        if len(ad_triangles) and ad_triangles.ndim == 1:
            ad_triangles = np.reshape(ad_triangles, (1,3))
        if len(b_triangles) and b_triangles.ndim == 1:
            b_triangles = np.reshape(b_triangles, (1,3))



        ad_dist = []
        for triangle in ad_triangles:
            e1 = ad_edges[int(triangle[0])]
            e2 = ad_edges[int(triangle[1])]
            e3 = ad_edges[int(triangle[2])]

            d1 = round(get_edge_length(e1, ad_vertices), 1)
            d2 = round(get_edge_length(e2, ad_vertices), 1)
            d3 = round(get_edge_length(e3, ad_vertices), 1)

            ad_dist.append(max(d1, d2, d3))

        ad_dist = np.array(ad_dist)

        if not is_sorted(ad_dist):
            print('Something wrong in sorting?')
            exit(1)

        uni, counts = np.unique(ad_dist, return_counts=True)
        possible_n_perm = 1
        for val in counts:
            possible_n_perm = possible_n_perm * math.factorial(val)


        local_n_perm = min(n_perm, possible_n_perm)


        triangle_indices = np.array(list(range(len(ad_triangles))))


        all_permutations = []

        perm = 0
        while (perm < local_n_perm):

            new_permutation = []

            triangle_ptr = 0.
            for idx, val in enumerate(counts):

                permuted = np.random.permutation(triangle_indices[int(triangle_ptr):int(triangle_ptr+val)])
                new_permutation += list(permuted)
                triangle_ptr += val


            if new_permutation not in all_permutations:
                all_permutations.append(new_permutation)
                perm += 1


        perm_dirr = this_pert_dirr + '/Permutations_ad'
        if not os.path.isdir(perm_dirr):
            os.mkdir(perm_dirr)

        for perm in range(local_n_perm):

            this_perm_dirr = perm_dirr + '/Permutation'+str(perm)
            if not os.path.isdir(this_perm_dirr):
                os.mkdir(this_perm_dirr)

            this_triangles = ad_triangles[all_permutations[perm]]

            #check_len = []
            #for edge in this_edges:
            #    x1,y1 = ad_vertices[int(edge[0])]
            #    x2,y2 = ad_vertices[int(edge[1])]
            #    check_len.append(round(math.sqrt((x1-x2)**2 + (y1-y2)**2), 1))

            np.savetxt(this_perm_dirr+'/triangles.csv', this_triangles, delimiter=',', fmt='%d')


        ################################################################################
        ################################################################################

        b_dist = []
        for triangle in b_triangles:
            e1 = b_edges[int(triangle[0])]
            e2 = b_edges[int(triangle[1])]
            e3 = b_edges[int(triangle[2])]

            d1 = round(get_edge_length(e1, b_vertices), 1)
            d2 = round(get_edge_length(e2, b_vertices), 1)
            d3 = round(get_edge_length(e3, b_vertices), 1)

            b_dist.append(max(d1, d2, d3))

        b_dist = np.array(b_dist)

        if not is_sorted(b_dist):
            print('Something wrong in sorting?')
            exit(1)

        uni, counts = np.unique(b_dist, return_counts=True)
        possible_n_perm = 1
        for val in counts:
            possible_n_perm = possible_n_perm * math.factorial(val)

        local_n_perm = min(n_perm, possible_n_perm)

        triangle_indices = np.array(list(range(len(b_dist))))

        all_permutations = []

        perm = 0
        while (perm < local_n_perm):

            new_permutation = []

            triangle_ptr = 0.
            for idx, val in enumerate(counts):

                permuted = np.random.permutation(triangle_indices[int(triangle_ptr):int(triangle_ptr+val)])
                new_permutation += list(permuted)
                triangle_ptr += val


            if new_permutation not in all_permutations:
                all_permutations.append(new_permutation)
                perm += 1


        perm_dirr = this_pert_dirr + '/Permutations_b'
        if not os.path.isdir(perm_dirr):
            os.mkdir(perm_dirr)

        for perm in range(local_n_perm):

            this_perm_dirr = perm_dirr + '/Permutation'+str(perm)
            if not os.path.isdir(this_perm_dirr):
                os.mkdir(this_perm_dirr)

            this_triangles = b_triangles[all_permutations[perm]]

            #check_len = []
            #for triangle in this_triangles:
            #    x1,y1 = ad_vertices[int(triangle[0])]
            #    x2,y2 = ad_vertices[int(triangle[1])]
            #    check_len.append(round(math.sqrt((x1-x2)**2 + (y1-y2)**2), 1))

            np.savetxt(this_perm_dirr+'/triangles.csv', this_triangles, delimiter=',', fmt='%d')


            








            








