import numpy as np
import pickle
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed


islets = pickle.load(open('far_unmatched_islets.p', 'rb'))


is_sorted = lambda a: np.all(a[:-1] <= a[1:])

dirr = 'stochastic'

thresh_dirr = 'Thresholds'
data_dirr = 'DevData/all_islets'

n_pert = 50

n_perm = 50


def single_max_pert(max_pert):

    all_ad_mantle_info = dict()

    for fprefix in islets:
    
        islet_dirr = dirr + '/' + fprefix
    
        print('\nProcessing', fprefix, max_pert)
    
    
        ##########################################
        # ad loops
    
        pert_dirr = islet_dirr + '/Perturbations_ad'
    
        max_pert_dirr = pert_dirr + '/Perturbation'+str(max_pert)
    
        # Cluster pts are not perturbed
        b_pts_file = '../VerticesEdges/' + fprefix + '_bvertices_PH.csv'
        ad_pts_file = '../VerticesEdges/' + fprefix + '_advertices_PH.csv'
    
        b_locs = np.loadtxt(b_pts_file, delimiter=',')
        ad_locs = np.loadtxt(ad_pts_file, delimiter=',')
    
        mantle_info = dict()
    
        for pert in range(n_pert):
    
            this_pert_dirr = max_pert_dirr + '/Perturbation' + str(pert)
    
            perm_dirr = this_pert_dirr + '/Permutations'
    
            for perm in range(n_perm):
    
                print(pert, perm, end='\r')
    
                this_perm_dirr = perm_dirr + '/Permutation'+str(perm)
                if not os.path.isdir(this_perm_dirr):
                    continue
    
                this_edge_file = this_perm_dirr + '/edges.csv'
                this_edge = np.loadtxt(this_edge_file, delimiter=',')
                
                undead_file = this_perm_dirr + '/adundead_PH.csv'
                if not os.path.isfile(undead_file):
                    continue
    
                undeads = open(undead_file, 'r')
    
                for line in undeads:
                    line = line.strip('\n')
                    line = line.split(':')
                    cluster = line[1].split(',')[:-1]
                    if not len(cluster):
                        continue
    
                    cluster = frozenset([int(x) for x in cluster])

                    loop = line[0].split(',')[:-1]
                    loop = [int(x) for x in loop]
                    edge_vertices = frozenset(np.array([this_edge[x] for x in loop], dtype=int).flatten())

                    if cluster not in mantle_info:
                        mantle_info[cluster] = [edge_vertices]
                    else:
                        if len(edge_vertices) < len(mantle_info[cluster][0]):
                            mantle_info[cluster] = [edge_vertices]
                        elif len(edge_vertices) == len(mantle_info[cluster][0]):
                            if edge_vertices not in mantle_info[cluster]:
                                mantle_info[cluster].append(edge_vertices)


    
                    #plt.scatter(b_locs[:,0], b_locs[:,1], color='blue', alpha=0.5)
                    #plt.scatter(ad_locs[:,0], ad_locs[:,1], color='red', alpha=0.5)
                    #for edge in edge_vertices:
                    #    p1, p2 = edge
                    #    x1, y1 = ad_locs[p1]
                    #    x2, y2 = ad_locs[p2]
                    #    plt.plot([x1, x2], [y1, y2], color='black')
                    #plt.show()
    
                    #if cluster not in mantle_info:
                    #    mantle_info[cluster] = edge_vertices
                    #else:
                    #    if edge_vertices.shape[0] < mantle_info[cluster].shape[0]:
                    #        mantle_info[cluster] = edge_vertices

    
    
        all_ad_mantle_info[fprefix] = mantle_info

    return all_ad_mantle_info

max_perts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

#for max_pert in max_perts:
#    single_max_pert(max_pert)


num_cores = 5
results = Parallel(n_jobs=num_cores, verbose = 12)\
                   (delayed(single_max_pert)\
                       (max_pert) for max_pert in max_perts)


all_all_ad_mantle_info = dict()

for idx, max_pert in enumerate(max_perts):
    all_all_ad_mantle_info[max_pert] = results[idx]


pickle.dump(all_all_ad_mantle_info, open('all_ad_mantle_info.p', 'wb'))

    
