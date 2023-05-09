import pickle
import os
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
import subprocess
from joblib import Parallel, delayed



# THIS ASSUMES THAT stochastic_far_unmatched.py has been executed

#islets = pickle.load(open('far_unmatched_islets.p', 'rb'))


is_sorted = lambda a: np.all(a[:-1] <= a[1:])

dirr = 'stochastic'


n_pert = 50

n_perm = 50


def single_job(job):

    fprefix = job[0]
    typ = job[1]
    pert = job[2]

    #print(fprefix)

    print(typ, fprefix, pert)

    for max_pert in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:


        islet_dirr = dirr + '/' + fprefix

        ##########################################
        # ad loops

        if typ == 'ad':
            pert_dirr = islet_dirr + '/Perturbations_ad'
        elif typ == 'b':
            pert_dirr = islet_dirr + '/Perturbations_b'


        max_pert_dirr = pert_dirr + '/Perturbation'+str(max_pert)

        # Cluster pts are not perturbed
        if typ == 'ad':
            cluster_pts_file = '../VerticesEdges/' + fprefix + '_bvertices_PH.csv'
        elif typ == 'b':
            cluster_pts_file = '../VerticesEdges/' + fprefix + '_advertices_PH.csv'

        #print(pert, end='\r')

        this_pert_dirr = max_pert_dirr + '/Perturbation' + str(pert)

        if typ == 'ad':
            loop_verts_file = this_pert_dirr + '/advertices_PH.csv'
            loop_edges_file = this_pert_dirr + '/adedges_PH.csv'
        elif typ == 'b':
            loop_verts_file = this_pert_dirr + '/bvertices_PH.csv'
            loop_edges_file = this_pert_dirr + '/bedges_PH.csv'

        if not os.path.isfile(loop_verts_file):
            continue

        if not os.path.isfile(loop_edges_file):
            continue

        flag_processed = 1
        perm_dirr = this_pert_dirr + '/Permutations'
        if not os.path.isdir(perm_dirr):
            os.mkdir(perm_dirr)



        flag_processed = 1
        for perm in range(n_perm):
            this_perm_dirr = perm_dirr + '/Permutation'+str(perm)
            if not os.path.isdir(this_perm_dirr):
                flag_processed = 0
                os.mkdir(this_perm_dirr)
                break

            if typ == 'ad':
                new_undead_file = this_perm_dirr + '/adundead_PH.csv'
            elif typ == 'b':
                new_undead_file = this_perm_dirr + '/bundead_PH.csv'

            if os.path.isfile(new_undead_file):
                continue
            else:
                flag_processed = 0
                break

        if flag_processed == 1:
            continue

        if typ == 'ad':
            # Get triangles
            pert_triangles_file = this_pert_dirr + '/adtriangles_PH.csv'
        elif typ == 'b':
            pert_triangles_file = this_pert_dirr + '/btriangles_PH.csv'

        ## MAKE PERMUTATIONS OF EDGES
        loop_vertices = np.loadtxt(loop_verts_file, delimiter=',')
        loop_edges = np.loadtxt(loop_edges_file, delimiter=',')

        loop_triangles = np.loadtxt(pert_triangles_file, delimiter=',')
        
        if not len(loop_triangles):
            continue

        if loop_triangles.ndim == 1:
            loop_triangles = np.reshape(loop_triangles, (1, 3))



        loop_dist = []
        for edge in loop_edges:

            e1, e2 = edge

            x1,y1 = loop_vertices[int(e1)]
            x2,y2 = loop_vertices[int(e2)]

            loop_dist.append(round(math.sqrt((x1-x2)**2 + (y1-y2)**2), 0))

        loop_dist = np.array(loop_dist)

        uni, counts = np.unique(loop_dist, return_counts=True)
        possible_n_perm = 1
        for val in counts:
            possible_n_perm = possible_n_perm * math.factorial(val)

        local_n_perm = min(n_perm, possible_n_perm)
        edge_indices = np.array(list(range(len(loop_edges))))
        perm_map = np.zeros(len(loop_edges), dtype=int)

        all_permutations = []

        #print('Permuting...')

        perm_dirr = this_pert_dirr + '/Permutations'
        if not os.path.isdir(perm_dirr):
            os.mkdir(perm_dirr)

        perm = 0
        while (perm < local_n_perm):

            #print(perm, end='\r')

            new_permutation = []

            edge_ptr = 0.
            for idx, val in enumerate(counts):
                permuted = np.random.permutation(edge_indices[int(edge_ptr):int(edge_ptr+val)])
                new_permutation += list(permuted)
                edge_ptr += val

            if new_permutation in all_permutations:
                continue

            all_permutations.append(new_permutation)

            this_edges = loop_edges[all_permutations[perm]]

            for new_idx, old_idx in enumerate(all_permutations[perm]):
                perm_map[old_idx] = new_idx

            this_perm_dirr = perm_dirr + '/Permutation'+str(perm)
            if not os.path.isdir(this_perm_dirr):
                os.mkdir(this_perm_dirr)

            # Save permutation
            np.savetxt(this_perm_dirr+'/permutation.csv'\
                                , all_permutations[perm]\
                                , delimiter=','\
                                , fmt='%d')


            # Save permuted edges
            this_edge_file = this_perm_dirr + '/edges.csv'
            np.savetxt(this_edge_file, this_edges, delimiter=',', fmt='%d')

            # Modify triangle boundaries with new edge-indices
            new_triangles = []
            triangle_dia = []
            for triangle in loop_triangles:
                e1, e2, e3 = triangle

                new_e1 = perm_map[int(e1)]
                new_e2 = perm_map[int(e2)]
                new_e3 = perm_map[int(e3)]

                new_triangle = [new_e1, new_e2, new_e3]

                # BOUNDARY HAS TO BE SORTED!!!!!!
                new_triangles.append(sorted(new_triangle, key=int))

                triangle_dia.append(max(new_triangle))

            new_triangles = np.array(new_triangles)

            triangle_dia = np.array(triangle_dia, dtype=int)

            idxs = np.argsort(triangle_dia)

            new_triangles = new_triangles[idxs]

            this_triangles_file = this_perm_dirr + '/triangles.csv'

            np.savetxt(this_triangles_file, new_triangles, delimiter=',', fmt='%d')

            perm += 1

            # Save undead
            if typ == 'ad':
                new_undead_file = this_perm_dirr + '/adundead_PH.csv'
            elif typ == 'b':
                new_undead_file = this_perm_dirr + '/bundead_PH.csv'

            subprocess.run(["../get_undead.o"\
                        , this_edge_file\
                        , this_triangles_file\
                        , new_undead_file\
                        , cluster_pts_file\
                        , loop_verts_file\
                        ])


#for islet in islets:
#    single_islet(islet)
#    exit()


far_ad_mantles_file = open('far_ad_mantles.csv', 'r')
far_ad_mantles = far_ad_mantles_file.readline().split(',')[:-1]
far_ad_mantles = frozenset(far_ad_mantles)

far_b_mantles_file = open('far_b_mantles.csv', 'r')
far_b_mantles = far_b_mantles_file.readline().split(',')[:-1]
far_b_mantles = frozenset(far_b_mantles)


jobs = []

for fprefix in far_ad_mantles:

    for pert in range(n_pert):

        jobs.append([fprefix, 'ad', pert])


for fprefix in far_b_mantles:

    for pert in range(n_pert):

        jobs.append([fprefix, 'b', pert])


#for job in jobs:
#    single_job(job)
#    input('w')

num_cores = 8
Parallel(n_jobs=6, verbose = 12)\
                   (delayed(single_job)\
                       (job) for job in jobs)



