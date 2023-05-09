import pickle
import os
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
import subprocess



# THIS ASSUMES THAT stochastic_far_unmatched.py has been executed

islets = pickle.load(open('far_unmatched_islets.p', 'rb'))


is_sorted = lambda a: np.all(a[:-1] <= a[1:])

dirr = 'stochastic'

thresh_dirr = 'Thresholds'
data_dirr = 'DevData/all_islets'

n_pert = 50

n_perm = 50


fcount = 0

for fprefix in islets:

    print(fprefix)

    fcount += 1

    islet_dirr = dirr + '/' + fprefix


    ##########################################
    # ad loops

    pert_dirr = islet_dirr + '/Perturbations_ad'

    max_pert = 0.5
    max_pert_dirr = pert_dirr + '/Perturbation'+str(max_pert)

    # Cluster pts are not perturbed
    cluster_pts_file = '../VerticesEdges/' + fprefix + '_bvertices_PH.csv'


    print(fprefix)

    for pert in range(n_pert):

        print(pert, end='\r')

        this_pert_dirr = max_pert_dirr + '/Perturbation' + str(pert)

        loop_verts_file = this_pert_dirr + '/advertices_PH.csv'
        loop_edges_file = this_pert_dirr + '/adedges_PH.csv'

        if not os.path.isfile(loop_verts_file):
            continue

        if not os.path.isfile(loop_edges_file):
            continue

        ## MAKE PERMUTATIONS OF EDGES
        ad_vertices = np.loadtxt(loop_verts_file, delimiter=',')
        ad_edges = np.loadtxt(loop_edges_file, delimiter=',')

        # Save triangles
        new_triangles_file = this_pert_dirr + '/adtriangles_PH.csv'

        subprocess.run(["../get_triangles.o"\
                        , loop_verts_file\
                        , loop_edges_file\
                        , cluster_pts_file\
                        , new_triangles_file\
                        ])

        ad_triangles = np.loadtxt(new_triangles_file, delimiter=',')

        ad_dist = []
        for edge in ad_edges:

            e1, e2 = edge

            x1,y1 = ad_vertices[int(e1)]
            x2,y2 = ad_vertices[int(e2)]

            ad_dist.append(round(math.sqrt((x1-x2)**2 + (y1-y2)**2), 0))

        ad_dist = np.array(ad_dist)

        uni, counts = np.unique(ad_dist, return_counts=True)
        possible_n_perm = 1
        for val in counts:
            possible_n_perm = possible_n_perm * math.factorial(val)


        local_n_perm = min(n_perm, possible_n_perm)
        edge_indices = np.array(list(range(len(ad_edges))))
        perm_map = np.zeros(len(ad_edges),)

        all_permutations = []

        print('Permuting...')

        perm_dirr = this_pert_dirr + '/Permutations'
        if not os.path.isdir(_perm_dirr):
            os.mkdir(perm_dirr)

        perm = 0
        while (perm < local_n_perm):

            print(perm, end='\r')

            new_permutation = []

            edge_ptr = 0.
            for idx, val in enumerate(counts):
                permuted = np.random.permutation(edge_indices[int(edge_ptr):int(edge_ptr+val)])
                new_permutation += list(permuted)
                edge_ptr += val

            if new_permutation in all_permutations:
                continue

            all_permutations.append(new_permutation)

            this_edges = ad_edges[all_permutations[perm]]

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
            for triangle in ad_triangles:
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
            new_undead_file = this_perm_dirr + '/adundead_PH.csv'

            subprocess.run(["../get_undead.o"\
                        , loop_edges_file\
                        , new_triangles_file\
                        , new_undead_file\
                        , cluster_pts_file\
                        , loop_verts_file\
                        ])


        ## plot cycles
        ## Check plot
        #cluster_locs = np.loadtxt(cluster_pts_file, delimiter=',')
        #loop_locs = np.loadtxt(loop_verts_file, delimiter=',')

        #geom_loop = [111.68,166.922\
        #        ,116.897,168.882\
        #        ,129.645,148.534\
        #        ,130.063,141.697\
        #        ,138.219,131.029\
        #        ,145.775,117.367\
        #        ,132.062,102.278\
        #        ,109.201,100.373\
        #        ,113.275,86.688\
        #        ,93.6056,71.4107\
        #        ,78.55,78.1606\
        #        ,75.684,86.2747\
        #        ,76.3051,103.355\
        #        ,92.7404,122.844\
        #        ,91.3707,130.286\
        #        ,102.697,144.821\
        #        ,99.3922,158.574\
        #        ,96.7553,163.666\
        #        ,89.5886,182.222\
        #        ,110.63,174.722\
        #        ,111.68,166.922]
        #
        #nn = len(geom_loop)

        #loop_edges = np.loadtxt(loop_edges_file, delimiter=',')
        #undead_loops = open(new_undead_file, 'r')
        #for line in undead_loops:
        #    line = line.strip('\n')

        #    title = line

        #    line = line.split(':')[0]
        #    line = line.split(',')
        #    line = line[:-1]
        #    line = [int(x) for x in line]

        #    plt.scatter(cluster_locs[:,0], cluster_locs[:,1], color='blue', alpha=0.5)
        #    plt.scatter(loop_locs[:,0], loop_locs[:,1], color='red', alpha=0.5)

        #    bb = 0
        #    while (bb < nn-4):
        #        x1 = geom_loop[bb]
        #        y1 = geom_loop[bb+1]

        #        x2 = geom_loop[bb+2]
        #        y2 = geom_loop[bb+3]

        #        plt.plot([x1, x2], [y1, y2], color='black')

        #        bb += 2



        #    for edge in line:
        #        v1, v2 = loop_edges[edge]

        #        x1, y1 = loop_locs[int(v1)]
        #        x2, y2 = loop_locs[int(v2)]
        #        plt.plot([x1, x2], [y1, y2], lw=4, alpha=0.7, color='orange')

        #    plt.title(line)
        #    plt.show()
        #    plt.cla()
        #    plt.clf()

        #undead_loops.close()
        #exit()




