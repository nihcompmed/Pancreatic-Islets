import numpy as np
import pickle
import  matplotlib.pyplot as plt
import islet_helper as ih
import math
import os
import matplotlib.pyplot
import networkx as nx


def check_condition(dist, thresh):

    #if dist <= round(thresh, 1)+1:
    if dist <= thresh:
        return 1
    else:
        return 0


#far_ad_mantle_info = pickle.load(open('Stochastic_NEW/all_ad_mantle_info.p', 'rb'))
#max_perts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

ff = open('all_islets_info_diabetic.csv', 'r')
n_cat = 2
all_prefix = []

cats = []


ff.readline()
for line in ff:

    line = line.split(',')
    prefix = line[3]
    all_prefix.append(prefix)
    if line[0] not in cats:
        cats.append(line[0])
    
ff.close()

far_ad_mantle = open('far_ad_mantles.csv', 'w')
far_b_mantle = open('far_b_mantles.csv', 'w')


count = 0

#all_results_file = 'matching_results.csv' 
#all_results = open(all_results_file, 'w')


all_far = []
all_close = []
n_unmatched = 0

max_b_in_mantle = dict()
max_ad_in_mantle = dict()

max_b_in_mantle_PH = dict()
max_ad_in_mantle_PH = dict()

far_islets = dict()

for cat in cats:
    max_b_in_mantle[cat] = []
    max_ad_in_mantle[cat] = []
    
    max_b_in_mantle_PH[cat] = []
    max_ad_in_mantle_PH[cat] = []

for fprefix in all_prefix:

    if count%1000 == 0:
        print(count, end='\r')
    count += 1

    v_file = 'VerticesEdges/' + fprefix + '.vertices'
    e_file = 'VerticesEdges/' + fprefix + '.edges'

    if not os.path.isfile(v_file) or not os.path.isfile(e_file):
        continue

    cat = fprefix.split('_')[0][-1]

    ad_PH_locs_file = 'VerticesEdges/' + fprefix + '_advertices_PH.csv'
    b_PH_locs_file = 'VerticesEdges/' + fprefix + '_advertices_PH.csv'

    ad_PH_locs = np.loadtxt(ad_PH_locs_file, delimiter=',')
    b_PH_locs = np.loadtxt(b_PH_locs_file, delimiter=',')


    # Get all edges (Islet numbers)
    vertices = open(v_file, 'r')
    v_islnum = []
    for line in vertices:
        line = line.strip('\n')
        line = line.split('\t')
        v_islnum.append(int(line[-2]))

    n_total_cells = len(v_islnum)

    edges = np.loadtxt(e_file, delimiter=',', dtype=int)
    if edges.ndim == 1:
        edges = np.reshape(edges, (1, 2))
    G = nx.Graph()
    G.add_nodes_from(v_islnum)

    for edge in edges:
        G.add_edge(edge[0], edge[1])



    far_scores = []
    close_scores = []

    # Load G(r) thresholds
    thresh_file = 'Thresholds/'+fprefix + '.smooth2.thresh.minBetPeak2And3.dat'
    threshs_data = open(thresh_file, 'r')
    threshs = threshs_data.readline()
    threshs = threshs.strip('\n')
    threshs = threshs.split('\t')
    ad_thresh = float(threshs[1])
    b_thresh = float(threshs[2])
    threshs_data.close()


    # Load geom to PH map
    # Map is of form dict[isl_num] = {'x', 'y', 'index'}
    ad_vert_geom_to_PH_file = 'VerticesEdges/' + fprefix + '_advertices_geom_to_PH.p'
    ad_map = pickle.load( open( ad_vert_geom_to_PH_file, "rb" ) )
    
    b_vert_geom_to_PH_file = 'VerticesEdges/' + fprefix + '_bvertices_geom_to_PH.p'
    b_map = pickle.load( open( b_vert_geom_to_PH_file, "rb" ) )
    
    
    # Geometric loops
    # Each line is of the form Islcell1, Islcell2, ...:x1,y1,x2,y2,.... of the loop around the cells.
    # (n,n+1), (n+1, n+2) ... are edges of the loop
    ad_geom_loops = 'Loops/' + fprefix + '.admantles'
    b_geom_loops = 'Loops/' + fprefix + '.bmantles'
    
    # PH loops
    # list [loops, info], loops is list of loops and info is dictionary dict[loop] = {cells inside loop}
    # loops are stored as geometric coordinates of edges
    ad_PH_loops_file = 'Loops/' + fprefix + '_admantles_PH.p'
    b_PH_loops_file = 'Loops/' + fprefix + '_bmantles_PH.p'
    
    ad_PH = pickle.load(open(ad_PH_loops_file, 'rb'))


    b_PH = pickle.load(open(b_PH_loops_file, 'rb'))
    
    ad_PH_loops = ad_PH[0]
    ad_PH_info = ad_PH[1]
    
    b_PH_loops = b_PH[0]
    b_PH_info = b_PH[1]
    
    if os.path.isfile(ad_geom_loops):
    
        # Get b-sets that have geometric loop around them
        far = []
        close = []
        vv = open(ad_geom_loops, 'r')
        for line in vv:
            line = line.strip('\n')
            line = line.split(':')
        
            geom_loop = line[1]
            geom_loop = geom_loop.split(',')
            geom_loop = geom_loop[:-1]
            geom_loop = [float(x) for x in geom_loop]
            geom_loop_pts = []
            nn = 0
            while (nn < len(geom_loop)):
                geom_loop_pts.append([geom_loop[nn], geom_loop[nn+1]])
                nn += 2
        
        
            bcells = line[0].split(',')
            bcells = bcells[:-1]
        
            bcells2 = []
        
            for x in bcells:
                y = int(x)
                if y < 100000:
                    bcells2.append(y)

            #bcells2 has isletnumbers of b-cells inside a geom loop
            this_b_graph = G.subgraph([int(x) for x in bcells2])

            this_b_components = nx.connected_components(this_b_graph)

            maxx_comp_size = 0

            for comp in this_b_components:
                if len(comp) > maxx_comp_size:
                    maxx_comp_size = len(comp)

            #max_b_in_mantle[stage].append([n_total_cells, maxx_comp_size])
        
            # Transform to PH vertex indexing
            geom_cells = [b_map[cell]['index'] for cell in bcells2]
        
            geom_cells = frozenset(geom_cells)
        
            loop = line[1]
        
            flag_found = 0
        
            for ii in ad_PH_info:
        
                # Have to change later
                # Right now this enforces that each set of b-cells is associated to a unique PH loop
        
                if len(ii) == 0:
                    continue
        
                PH_cells = ad_PH_info[ii]
        
                PH_cells = frozenset(PH_cells)
                
                diff1 = geom_cells.difference(PH_cells)
                diff2 = PH_cells.difference(geom_cells)
                not_common = diff1.union(diff2)
                if not len(not_common):
        
                    flag_found = 1
                    min_dist = math.inf
                    for jj in ii:
                        this_loop = ad_PH_loops[jj]
                        this_loop_pts = []
                        for edge in this_loop:
                            this_loop_pts.append(edge[0])
                            this_loop_pts.append(edge[1])
                        
                        ddist = ih.max_min_dist_loop_loop(geom_loop_pts, this_loop_pts)
                        if ddist < min_dist:
                            min_dist = ddist
        
                    #print('matched! with loop distance', min_dist)
        
                    if flag_found and check_condition(min_dist, ad_thresh):
                        close_scores.append(min_dist)
                        # This measure has to be the same because exactly same set of b-cells
                        # is inside the PH loops
                        #max_b_in_mantle_PH[stage].append([n_total_cells, maxx_comp_size])
                        break
        
            if (not flag_found) or ((flag_found) and (not check_condition(min_dist, ad_thresh))):

                far_ad_mantle.write(fprefix + ',')

                ## Look in far_ad_mantles dictionary
                #flag_found_far = 0
                #for max_pert in max_perts:
                #    this_dict = far_ad_mantle_info[max_pert][fprefix]
                #    if frozenset(geom_cells) not in this_dict:
                #        continue
                #    this_PH_loop_list = this_dict[frozenset(geom_cells)]
                #    for this_PH_loop in this_PH_loop_list:
                #        loop_pts = []
                #        for pt in this_PH_loop:
                #            loop_pts.append(ad_PH_locs[pt])
                #        loop_pts = np.array(loop_pts)
                #        ddist = ih.max_min_dist_loop_loop(geom_loop_pts, loop_pts)
                #        if ddist <= ad_thresh:
                #            max_b_in_mantle_PH[stage].append([n_total_cells, maxx_comp_size])
                #            close_scores.append(ddist)
                #            flag_found_far = 1
                #            break

                #    if flag_found_far == 1:
                #        break

                #if not flag_found_far:
                #    if not flag_found:
                #        n_unmatched += 1
                #    else:
                #        far_islets[fprefix] = 1
                #        far_scores.append(min_dist)
        
    #all_close += close_scores
    #all_far += far_scores

    #all_results.write(fprefix + ", ad-mantles:")
    #for entry in match:
    #    all_results.write(str(entry)+',')
    #all_results.write('\n')



    if os.path.isfile(b_geom_loops):
    
        # Get ad-sets that have geometric loop around them
        match = []
        vv = open(b_geom_loops, 'r')
        for line in vv:
            line = line.strip('\n')
            line = line.split(':')
        
            geom_loop = line[1]
            geom_loop = geom_loop.split(',')
            geom_loop = geom_loop[:-1]
            geom_loop = [float(x) for x in geom_loop]
            geom_loop_pts = []
            nn = 0
            while (nn < len(geom_loop)):
                geom_loop_pts.append([geom_loop[nn], geom_loop[nn+1]])
                nn += 2
        
        
            adcells = line[0].split(',')
            adcells = adcells[:-1]
        
        
            adcells2 = []
        
            for x in adcells:
                y = int(x)
                if y < 100000:
                    adcells2.append(y)

            #adcells2 has isletnumbers of b-cells inside a geom loop
            this_ad_graph = G.subgraph([int(x) for x in adcells2])

            this_ad_components = nx.connected_components(this_ad_graph)

            maxx_comp_size = 0

            for comp in this_ad_components:
                if len(comp) > maxx_comp_size:
                    maxx_comp_size = len(comp)

            #max_ad_in_mantle[stage].append([n_total_cells, maxx_comp_size])
        
            geom_cells = [ad_map[cell]['index'] for cell in adcells2]
        
            geom_cells = frozenset(geom_cells)
        
            loop = line[1]
        
            flag_found = 0
        
            for ii in b_PH_info:
        
                # Have to change later
                # Right now this enforces that each set of b-cells is associated to a unique PH loop
        
                if len(ii) == 0:
                    continue
        
                PH_cells = b_PH_info[ii]
        
                PH_cells = frozenset(PH_cells)
                
                diff1 = geom_cells.difference(PH_cells)
                diff2 = PH_cells.difference(geom_cells)
                not_common = diff1.union(diff2)
                if not len(not_common):
        
                    flag_found = 1
                    min_dist = math.inf
                    for jj in ii:
                        this_loop = b_PH_loops[jj]
                        this_loop_pts = []
                        for edge in this_loop:
                            this_loop_pts.append(edge[0])
                            this_loop_pts.append(edge[1])
                        
                        ddist = ih.max_min_dist_loop_loop(geom_loop_pts, this_loop_pts)
                        if ddist < min_dist:
                            min_dist = ddist
        
                    #print('matched! with loop distance', min_dist)
        
                    if flag_found:
                        if not check_condition(min_dist, b_thresh):
                            far_islets[fprefix] = 1
                            far_scores.append(min_dist)
                            far_b_mantle.write(fprefix + ',')
                        else:
                            close_scores.append(min_dist)
                            #max_ad_in_mantle_PH[stage].append([n_total_cells, maxx_comp_size])
                        break
        
            if not flag_found:
                far_b_mantle.write(fprefix + ',')
                far_islets[fprefix] = 0
                n_unmatched += 1
    
    
    all_close += close_scores
    all_far += far_scores
    

    #all_results.write(fprefix + ", b-mantles:")
    #for entry in match:
    #    all_results.write(str(entry)+',')
    #all_results.write('\n')

far_ad_mantle.close()
far_b_mantle.close()
exit()

#############################################################################
#############################################################################
pickle.dump(far_islets, open('far_unmatched_islets_after_stochastic.p', 'wb'))

for stage in range(4):
    max_b_in_mantle[stage]     = np.array(max_b_in_mantle[stage])
    max_ad_in_mantle[stage]    = np.array(max_ad_in_mantle[stage])
                                                                    
    max_b_in_mantle_PH[stage]  = np.array(max_b_in_mantle_PH[stage])
    max_ad_in_mantle_PH[stage] = np.array(max_ad_in_mantle_PH[stage])
    
max_components_in_mantle = [\
                             max_b_in_mantle\
                            ,max_ad_in_mantle\
                            ,max_b_in_mantle_PH\
                            ,max_ad_in_mantle_PH\
                            ]

pickle.dump(max_components_in_mantle, open('max_comp_size_in_mantle_after_stochastic.p', 'wb'))
#############################################################################
#############################################################################
    
#all_results.close() 

#all_scores = np.array(all_scores)

#unmatch_idxs = np.argwhere(all_scores == -1)
#match_idxs = np.argwhere(all_scores != -1)
#print('Unmatched = ', len(unmatch_idxs))
#print('Matched = ', len(match_idxs))

plt.hist(all_close, bins=100, color='blue', alpha=0.5, label='Close')

plt.hist(all_far, bins=100, color='red', alpha=0.5, label='Far')

plt.legend()

plt.xlabel('distance between geom and PH loop')
plt.ylabel('count')
plt.yscale('log', base = 2)
plt.title('Close '+ str(len(all_close)) + ', Far ' + str(len(all_far)) + ', Unmatched ' + str(n_unmatched))
#plt.show()
plt.savefig('figures/geom_vs_PH_loops_after_stochastic_dev.pdf')



