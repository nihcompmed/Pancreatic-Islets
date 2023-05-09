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

dirr = 'stochastic_only_perms'
if not os.path.isdir(dirr):
    os.mkdir(dirr)

thresh_dirr = 'Thresholds'
data_dirr = 'DevData/all_islets'

n_perm = 100
fcount = 0

for fprefix in islets:

    fcount += 1

    islet_dirr = dirr + '/' + fprefix

    if not os.path.isdir(islet_dirr):
        os.mkdir(islet_dirr)

    perm_dirr = islet_dirr + '/Permutations'
    if not os.path.isdir(perm_dirr):
        os.mkdir(perm_dirr)

    print(fprefix)

    ad_vert_PH_file = 'VerticesEdges/' + fprefix + '_advertices_PH.csv'
    b_vert_PH_file = 'VerticesEdges/' + fprefix + '_bvertices_PH.csv'

    ad_vertices = np.loadtxt(ad_vert_PH_file, delimiter=',')
    b_vertices = np.loadtxt(b_vert_PH_file, delimiter=',')

    ad_edge_file = 'VerticesEdges/' + fprefix + '_adedges_PH.csv'
    b_edge_file = 'VerticesEdges/' + fprefix + '_bedges_PH.csv'

    ad_edges = np.loadtxt(ad_edge_file, delimiter=',')
    b_edges = np.loadtxt(b_edge_file, delimiter=',')

    ad_triangles_file = 'VerticesEdges/' + fprefix + '_adtriangles_PH.csv'
    b_triangles_file = 'VerticesEdges/' + fprefix + '_btriangles_PH.csv'




    ad_triangles = np.loadtxt(ad_triangles_file, delimiter=',')

    correct_ad_triangles = []
    for triangle in ad_triangles:
        if list(triangle) not in correct_ad_triangles:
            correct_ad_triangles.append(list(triangle))

    ad_triangles = np.array(correct_ad_triangles, dtype=int)


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
    for edge in ad_edges:

        e1, e2 = edge

        x1,y1 = ad_vertices[int(e1)]
        x2,y2 = ad_vertices[int(e2)]

        ad_dist.append(round(math.sqrt((x1-x2)**2 + (y1-y2)**2), 0))

    ad_dist = np.array(ad_dist)


    if not is_sorted(ad_dist):
        print('Something wrong in sorting?')
        exit(1)

    uni, counts = np.unique(ad_dist, return_counts=True)
    possible_n_perm = 1
    for val in counts:
        possible_n_perm = possible_n_perm * math.factorial(val)

    local_n_perm = min(n_perm, possible_n_perm)

    edge_indices = np.array(list(range(len(ad_edges))))

    perm_map = np.zeros(len(ad_edges),)

    all_permutations = []

    print('Permuting...')

    perm = 0
    while (perm < local_n_perm):


        print(perm, end='\r')

        new_permutation = []

        edge_ptr = 0.
        for idx, val in enumerate(counts):
            permuted = np.random.permutation(edge_indices[int(edge_ptr):int(edge_ptr+val)])
            new_permutation += list(permuted)
            edge_ptr += val

        if new_permutation not in all_permutations:
            all_permutations.append(new_permutation)
            perm += 1

    ad_perm_dirr = perm_dirr + '/Permutations_ad'
    if not os.path.isdir(ad_perm_dirr):
        os.mkdir(ad_perm_dirr)

    print('Permuted. Computing PH undead...')

    for perm in range(local_n_perm):

        print(perm, end='\r')

        this_perm_dirr = ad_perm_dirr + '/Permutation'+str(perm)
        if not os.path.isdir(this_perm_dirr):
            os.mkdir(this_perm_dirr)

        this_edges = ad_edges[all_permutations[perm]]

        for new_idx, old_idx in enumerate(all_permutations[perm]):
            perm_map[old_idx] = new_idx

        #print(all_permutations[perm])
        #print(perm_map)
        #exit()

        # Save permutation
        np.savetxt(this_perm_dirr+'/permutation.csv', all_permutations[perm], delimiter=',', fmt='%d')


        #check_len = []
        #for edge in this_edges:
        #    x1,y1 = ad_vertices[int(edge[0])]
        #    x2,y2 = ad_vertices[int(edge[1])]
        #    check_len.append(round(math.sqrt((x1-x2)**2 + (y1-y2)**2), 1))

        this_edge_file = this_perm_dirr + '/edges.csv'


        new_triangles = []
        triangle_dia = []
        for triangle in ad_triangles:
            e1, e2, e3 = triangle

            #  check vertices
            check_verts = []
            check_verts +=  list(ad_edges[int(e1)])
            check_verts +=  list(ad_edges[int(e2)])
            check_verts +=  list(ad_edges[int(e3)])
            check_verts = frozenset(check_verts)

            new_e1 = perm_map[int(e1)]
            new_e2 = perm_map[int(e2)]
            new_e3 = perm_map[int(e3)]

            new_triangle = [new_e1, new_e2, new_e3]

            # BOUNDARY HAS TO BE SORTED!!!!!!
            new_triangles.append(sorted(new_triangle, key=int))

            triangle_dia.append(max(new_triangle))

            #  check vertices
            check_verts2 = []
            check_verts2 +=  list(this_edges[int(new_e1)])
            check_verts2 +=  list(this_edges[int(new_e2)])
            check_verts2 +=  list(this_edges[int(new_e3)])
            check_verts2 = frozenset(check_verts2)

            if (check_verts != check_verts2):
                print("UNMATCHED")
                exit(1)

        new_triangles = np.array(new_triangles)


        triangle_dia = np.array(triangle_dia, dtype=int)

        idxs = np.argsort(triangle_dia)

        new_triangles = new_triangles[idxs]


        this_triangles_file = this_perm_dirr + '/triangles.csv'

        ##############################################################
        np.savetxt(this_edge_file, this_edges, delimiter=',', fmt='%d')
        np.savetxt(this_triangles_file, new_triangles, delimiter=',', fmt='%d')

        #check_ddist = []
        #for edge in this_edges:
        #    x1,y1 = ad_vertices[int(edge[0])]
        #    x2,y2 = ad_vertices[int(edge[1])]

        #    ddist = round(math.sqrt((x1-x2)**2 + (y1-y2)**2),1)
        #    check_ddist.append(ddist)

        #print(check_ddist)
        #exit()

        ## DUMMY CHECK -- REMOVE LATER
        #np.savetxt(this_edge_file, ad_edges, delimiter=',', fmt='%d')
        ## DUMMY CHECK -- REMOVE LATER
        #np.savetxt(this_triangles_file, ad_triangles, delimiter=',', fmt='%d')

        this_undead_file = this_perm_dirr + '/undead.csv'


        subprocess.run(["./ad_mantles_only_perms.o"\
                        , this_edge_file\
                        , this_triangles_file\
                        , this_undead_file\
                        ])
        


    ################################################################################
    ################################################################################

