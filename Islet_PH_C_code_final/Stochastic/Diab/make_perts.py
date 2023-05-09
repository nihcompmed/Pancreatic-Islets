import numpy as np
import pickle
import os
import subprocess
import math


far_ad_mantles_file = open('far_ad_mantles.csv', 'r')
far_ad_mantles = far_ad_mantles_file.readline().split(',')[:-1]
far_ad_mantles = frozenset(far_ad_mantles)

far_b_mantles_file = open('far_b_mantles.csv', 'r')
far_b_mantles = far_b_mantles_file.readline().split(',')[:-1]
far_b_mantles = frozenset(far_b_mantles)



#Note: Perturb only cells in ad-graph from geom loops

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

stochastic_dirr = 'stochastic'
if not os.path.isdir(stochastic_dirr):
    os.mkdir(stochastic_dirr)

n_pert = 50


################################################
############# FAR ad-mantles ###################
################################################

for fprefix in far_ad_mantles:

    print(fprefix)
    
    islet_dirr = stochastic_dirr + '/' + fprefix
    if not os.path.isdir(islet_dirr):
        os.mkdir(islet_dirr)
    
    
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

    # Cluster pts are not perturbed
    cluster_pts_file = '../VerticesEdges/' + fprefix + '_bvertices_PH.csv'

    pert_dirr = islet_dirr + '/Perturbations_ad'
    if not os.path.isdir(pert_dirr):
        os.mkdir(pert_dirr)



    for max_pert in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:

        max_pert_dirr = pert_dirr + '/Perturbation'+str(max_pert)
        if not os.path.isdir(max_pert_dirr):
            os.mkdir(max_pert_dirr)


        # Perturb ad vertices
        n_pts = ad_verts.shape[0]

    
        for pert in range(n_pert):
    
            print(max_pert, pert, end='\r')
    
            this_pert_dirr = max_pert_dirr + '/Perturbation' + str(pert)
            if not os.path.isdir(this_pert_dirr):
                os.mkdir(this_pert_dirr)
    
            
            this_rand_perts = np.random.uniform(size=(n_pts, 2))
            delta_x = this_rand_perts[:,0] * max_pert * np.cos(2*math.pi*this_rand_perts[:,1])
            delta_y = this_rand_perts[:,0] * max_pert * np.sin(2*math.pi*this_rand_perts[:,1])
    
            
            new_locs = np.copy(ad_verts)
            new_locs[:,0] += delta_x
            new_locs[:,1] += delta_y
    
            new_vert_file = this_pert_dirr + '/advertices_PH.csv'
            np.savetxt(new_vert_file, new_locs, delimiter=',', fmt='%.2f')
    
            new_edge_file = this_pert_dirr + '/adedges_PH.csv'
    
            new_edges = np.copy(ad_edges)
            if new_edges.ndim == 1:
                new_edges = np.reshape(new_edges, (1, 2))
            edge_lengths = np.zeros(new_edges.shape[0],)
            for idx, edge in enumerate(new_edges):
                v1, v2 = edge
    
                x1, y1 = new_locs[int(v1)]
                x2, y2 = new_locs[int(v2)]
    
                dist = (x1-x2)**2 + (y1-y2)**2
    
                edge_lengths[idx] = dist
    
            sort_idxs = np.argsort(edge_lengths)
            new_edges = new_edges[sort_idxs]
    
            # Save sorted edges
            np.savetxt(new_edge_file, new_edges, delimiter=',', fmt='%d')
    
            # Save triangles
            new_triangles_file = this_pert_dirr + '/adtriangles_PH.csv'

            loop_verts_file = new_vert_file
            loop_edges_file = new_edge_file
            cluster_pts_file = b_vert_PH_file

            subprocess.run(["../get_triangles.o"\
                            , loop_verts_file\
                            , loop_edges_file\
                            , cluster_pts_file\
                            , new_triangles_file\
                            ])

            #subprocess.run(["../get_triangles.o"\
            #                , loop_vertices\
            #                , loop_edges\
            #                , cluster_pts\
            #                , new_triangles_file\
            #                ])
    
            # Save undead
            pert_undead_file = this_pert_dirr + '/adundead_PH.csv'

            subprocess.run(["../get_undead.o"\
                        , loop_edges_file\
                        , new_triangles_file\
                        , pert_undead_file\
                        , cluster_pts_file\
                        , loop_verts_file\
                        ])
    
################################################
############# FAR b-mantles ###################
################################################

for fprefix in far_b_mantles:

    print(fprefix)
    
    islet_dirr = stochastic_dirr + '/' + fprefix
    if not os.path.isdir(islet_dirr):
        os.mkdir(islet_dirr)
    
    
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
    
    if b_edges.shape[0] < 3:
        continue

    # Cluster pts are not perturbed
    cluster_pts_file = '../VerticesEdges/' + fprefix + '_advertices_PH.csv'

    pert_dirr = islet_dirr + '/Perturbations_b'
    if not os.path.isdir(pert_dirr):
        os.mkdir(pert_dirr)



    for max_pert in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:

        max_pert_dirr = pert_dirr + '/Perturbation'+str(max_pert)
        if not os.path.isdir(max_pert_dirr):
            os.mkdir(max_pert_dirr)


        # Perturb b vertices
        n_pts = b_verts.shape[0]

    
        for pert in range(n_pert):
    
            print(max_pert, pert, end='\r')
    
            this_pert_dirr = max_pert_dirr + '/Perturbation' + str(pert)
            if not os.path.isdir(this_pert_dirr):
                os.mkdir(this_pert_dirr)
    
            
            this_rand_perts = np.random.uniform(size=(n_pts, 2))
            delta_x = this_rand_perts[:,0] * max_pert * np.cos(2*math.pi*this_rand_perts[:,1])
            delta_y = this_rand_perts[:,0] * max_pert * np.sin(2*math.pi*this_rand_perts[:,1])
    
            
            new_locs = np.copy(b_verts)
            new_locs[:,0] += delta_x
            new_locs[:,1] += delta_y
    
            new_vert_file = this_pert_dirr + '/bvertices_PH.csv'
            np.savetxt(new_vert_file, new_locs, delimiter=',', fmt='%.2f')
    
            new_edge_file = this_pert_dirr + '/bedges_PH.csv'
    
            new_edges = np.copy(b_edges)
            if new_edges.ndim == 1:
                new_edges = np.reshape(new_edges, (1, 2))
            edge_lengths = np.zeros(new_edges.shape[0],)
            for idx, edge in enumerate(new_edges):
                v1, v2 = edge
    
                x1, y1 = new_locs[int(v1)]
                x2, y2 = new_locs[int(v2)]
    
                dist = (x1-x2)**2 + (y1-y2)**2
    
                edge_lengths[idx] = dist
    
            sort_idxs = np.argsort(edge_lengths)
            new_edges = new_edges[sort_idxs]
    
            # Save sorted edges
            np.savetxt(new_edge_file, new_edges, delimiter=',', fmt='%d')
    
            # Save triangles
            new_triangles_file = this_pert_dirr + '/btriangles_PH.csv'

            loop_verts_file = new_vert_file
            loop_edges_file = new_edge_file
            cluster_pts_file = b_vert_PH_file

            subprocess.run(["../get_triangles.o"\
                            , loop_verts_file\
                            , loop_edges_file\
                            , cluster_pts_file\
                            , new_triangles_file\
                            ])

            #subprocess.run(["../get_triangles.o"\
            #                , loop_vertices\
            #                , loop_edges\
            #                , cluster_pts\
            #                , new_triangles_file\
            #                ])
    
            # Save undead
            pert_undead_file = this_pert_dirr + '/bundead_PH.csv'

            subprocess.run(["../get_undead.o"\
                        , loop_edges_file\
                        , new_triangles_file\
                        , pert_undead_file\
                        , cluster_pts_file\
                        , loop_verts_file\
                        ])
    
