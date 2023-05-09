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


far_ad_mantle_info = pickle.load(open('Stochastic_NEW/all_ad_mantle_info.p', 'rb'))
max_perts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

ff = open('all_islets_info.csv', 'r')
all_prefix = []


ff.readline()
for line in ff:

    line = line.split(',')
    prefix = line[3]
    all_prefix.append(prefix)
    
ff.close()


count = 0

#all_results_file = 'matching_results.csv' 
#all_results = open(all_results_file, 'w')


all_far = []
all_close = []
n_unmatched = 0

#########################
max_b_in_mantle = dict()
max_ad_in_mantle = dict()

max_b_in_mantle_PH = dict()
max_ad_in_mantle_PH = dict()

#########################
num_b_in_mantle = dict()
num_ad_in_mantle = dict()

num_b_in_mantle_PH = dict()
num_ad_in_mantle_PH = dict()

#########################
comp_lens_b_in_mantle = dict()
comp_lens_ad_in_mantle = dict()

comp_lens_b_in_mantle_PH = dict()
comp_lens_ad_in_mantle_PH = dict()

#########################

smallest_thresh = 100000


PH_islet_comps_sizes_dict = dict()


far_islets = dict()

flag_plot_example = 1

for stage in range(4):
    max_b_in_mantle[stage] = []
    max_ad_in_mantle[stage] = []
    
    max_b_in_mantle_PH[stage] = []
    max_ad_in_mantle_PH[stage] = []

    #########################
    num_b_in_mantle[stage]           = []
    num_ad_in_mantle[stage]          = []

    num_b_in_mantle_PH[stage]        = []
    num_ad_in_mantle_PH[stage]       = []
                                       
    #########################[stage]   
    comp_lens_b_in_mantle[stage]     = []
    comp_lens_ad_in_mantle[stage]    = []
                                       
    comp_lens_b_in_mantle_PH[stage]  = []
    comp_lens_ad_in_mantle_PH[stage] = []
    
    #########################

    PH_islet_comps_sizes_dict[stage] = []
    

for fprefix in all_prefix:

    if count%1000 == 0:
        print(count, end='\r')
    count += 1

    v_file = 'VerticesEdges/' + fprefix + '.vertices'
    e_file = 'VerticesEdges/' + fprefix + '.edges'

    if not os.path.isfile(v_file) or not os.path.isfile(e_file):
        continue

    stage = int(fprefix.split('_')[0][-1])

    

    ad_PH_locs_file = 'VerticesEdges/' + fprefix + '_advertices_PH.csv'
    b_PH_locs_file = 'VerticesEdges/' + fprefix + '_bvertices_PH.csv'

    ad_PH_locs = np.loadtxt(ad_PH_locs_file, delimiter=',')
    b_PH_locs = np.loadtxt(b_PH_locs_file, delimiter=',')


    all_ad_cells = []
    all_b_cells = []

    # Get all edges (Islet numbers)
    vertices = open(v_file, 'r')
    v_islnum = []
    for line in vertices:
        line = line.strip('\n')
        line = line.split('\t')
        v_islnum.append(int(line[-2]))
        ctyp = line[-1]

        if ctyp == 'b':
            all_b_cells.append(int(line[-2]))
        else:
            all_ad_cells.append(int(line[-2]))
        

    n_total_cells = len(v_islnum)

    this_b_fraction = len(all_b_cells)/n_total_cells

    this_islet_characterization = (this_b_fraction, n_total_cells)

    PH_islet_comps_sizes_dict[stage].append({'feature':this_islet_characterization\
                                        , 'b_comp_lens':[]\
                                        , 'b_comp_lens_inside_mantle':[]\
                                        , 'b_comp_lens_not_inside_mantle':[]\
                                        , 'b_comp_lens_inside_PHmantle':[]\
                                        , 'ad_comp_lens':[]\
                                        , 'ad_comp_lens_inside_mantle':[]\
                                        , 'ad_comp_lens_not_inside_mantle':[]\
                                        , 'ad_comp_lens_inside_PHmantle':[]\
                                        })

    edges = np.loadtxt(e_file, delimiter=',', dtype=int)
    if edges.ndim == 1:
        edges = np.reshape(edges, (1, 2))
    G = nx.Graph()
    G.add_nodes_from(v_islnum)
    for edge in edges:
        G.add_edge(edge[0], edge[1])
    
    b_graph = G.subgraph(all_b_cells)
    ad_graph = G.subgraph(all_ad_cells)
    
    b_comps = nx.connected_components(b_graph)
    ad_comps = nx.connected_components(ad_graph)
    
    for b_comp in b_comps:
        PH_islet_comps_sizes_dict[stage][-1]['b_comp_lens'].append(len(b_comp))

    for ad_comp in ad_comps:
        PH_islet_comps_sizes_dict[stage][-1]['ad_comp_lens'].append(len(ad_comp))
    


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

    smallest_thresh = min([smallest_thresh, ad_thresh, b_thresh])


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

            num_b_in_mantle[stage].append(len(bcells2))

            this_comp_lens = []
            for comp in this_b_components:

                PH_islet_comps_sizes_dict[stage][-1]['b_comp_lens_inside_mantle'].append(len(comp))
                
                if len(comp) > maxx_comp_size:
                    maxx_comp_size = len(comp)

                this_comp_lens.append(len(comp))

            comp_lens_b_in_mantle[stage].append(this_comp_lens)

            # Number of bcells in this mantle

            max_b_in_mantle[stage].append([n_total_cells, maxx_comp_size])
        
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

                PH_loop_found = []

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
                            PH_loop_found = this_loop
        
                    #print('matched! with loop distance', min_dist)
        
                    if flag_found and check_condition(min_dist, ad_thresh):
                        close_scores.append(min_dist)
                        # This measure has to be the same because exactly same set of b-cells
                        # is inside the PH loops
                        num_b_in_mantle_PH[stage].append(len(geom_cells))
                        comp_lens_b_in_mantle_PH[stage].append(this_comp_lens)
                        max_b_in_mantle_PH[stage].append([n_total_cells, maxx_comp_size])

                        for comp in nx.connected_components(this_b_graph):
                            PH_islet_comps_sizes_dict[stage][-1]['b_comp_lens_inside_PHmantle'].append(len(comp))

                            # Just to plot an example
                            if flag_plot_example:
                                if len(comp) > 2:
                                    nn = 0
                                    xx = []
                                    yy = []
                                    while nn < len(geom_loop_pts):
                                        xx.append(geom_loop_pts[nn][0])
                                        yy.append(geom_loop_pts[nn][1])
                                        nn += 1
                                    
                                    for edge in PH_loop_found:
                                        x0,y0 = edge[0]
                                        x1,y1 = edge[1]
                                        plt.plot([x0, x1], [y0, y1], color='black', lw=4, alpha=0.4)
                                    plt.scatter(ad_PH_locs[:,0], ad_PH_locs[:,1], color='red')
                                    plt.scatter(b_PH_locs[:,0], b_PH_locs[:,1], color='blue')
                                    plt.plot(xx, yy, color='black', lw=2, ls='--')
                                    plt.axis('off')

                                    from matplotlib.lines import Line2D

                                    custom_lines1 = [
                                                    Line2D([0], [0], marker='o', color='None', markerfacecolor='b', markersize=15),
                                                    Line2D([0], [0], marker='o', color='None', markerfacecolor='r', markersize=15),
                                                    ]
                                    legend1 = plt.legend(custom_lines1, [
                                                                    r'$\beta$-cells',
                                                                    r'$\alpha\delta$-cells'],
                                                         fontsize=14, loc='lower left',
                                                         framealpha=0.4)


                                    custom_lines2 = [
                                                    Line2D([0], [0], color='black', ls='--', lw=2),
                                                    Line2D([0], [0], color='black', ls='-', lw=4, alpha=0.4),
                                                    ]
                                    legend2 = plt.legend(custom_lines2, [
                                                                    'Geom. mantle',
                                                                    'PH-mantle'], fontsize=14\
                                                                            , loc='upper left'\
                                                                            , framealpha=0.4)


                                    plt.gca().add_artist(legend1)
                                    plt.gca().add_artist(legend2)

                                    plt.savefig('example_geom_PH_ad_mantle.pdf', dpi=600)

                                    plt.cla()
                                    plt.clf()

                                                        

                        break
        
            if (not flag_found) or ((flag_found) and (not check_condition(min_dist, ad_thresh))):
                # Look in far_ad_mantles dictionary
                flag_found_far = 0
                for max_pert in max_perts:
                    this_dict = far_ad_mantle_info[max_pert][fprefix]
                    if frozenset(geom_cells) not in this_dict:
                        continue
                    this_PH_loop_list = this_dict[frozenset(geom_cells)]
                    for this_PH_loop in this_PH_loop_list:
                        loop_pts = []
                        for pt in this_PH_loop:
                            loop_pts.append(ad_PH_locs[pt])
                        loop_pts = np.array(loop_pts)
                        ddist = ih.max_min_dist_loop_loop(geom_loop_pts, loop_pts)
                        if ddist <= ad_thresh:

                            num_b_in_mantle_PH[stage].append(len(geom_cells))
                            comp_lens_b_in_mantle_PH[stage].append(this_comp_lens)
                            max_b_in_mantle_PH[stage].append([n_total_cells, maxx_comp_size])

                            close_scores.append(ddist)
                            flag_found_far = 1

                            for comp in nx.connected_components(this_b_graph):
                                PH_islet_comps_sizes_dict[stage][-1]['b_comp_lens_inside_PHmantle'].append(len(comp))

                            break

                    if flag_found_far == 1:
                        break

                if not flag_found_far:
                    if not flag_found:
                        n_unmatched += 1
                    else:
                        far_islets[fprefix] = 1
                        far_scores.append(min_dist)
        
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

            num_ad_in_mantle[stage].append(len(adcells2))
            this_comp_lens = []
            for comp in this_ad_components:

                PH_islet_comps_sizes_dict[stage][-1]['ad_comp_lens_inside_mantle'].append(len(comp))

                if len(comp) > maxx_comp_size:
                    maxx_comp_size = len(comp)

                this_comp_lens.append(len(comp))
            comp_lens_ad_in_mantle[stage].append(this_comp_lens)

            max_ad_in_mantle[stage].append([n_total_cells, maxx_comp_size])
        
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

                PH_loop_found = []
                
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
                            PH_loop_found = this_loop
        
                    #print('matched! with loop distance', min_dist)
        
                    if flag_found:
                        if not check_condition(min_dist, b_thresh):
                            far_islets[fprefix] = 1
                            far_scores.append(min_dist)
                        else:

                            for comp in nx.connected_components(this_ad_graph):
                                PH_islet_comps_sizes_dict[stage][-1]['ad_comp_lens_inside_PHmantle'].append(len(comp))

                            # Just to plot an example
                            if flag_plot_example:
                                if len(comp) > 2:
                                    nn = 0
                                    xx = []
                                    yy = []
                                    while nn < len(geom_loop_pts):
                                        xx.append(geom_loop_pts[nn][0])
                                        yy.append(geom_loop_pts[nn][1])
                                        nn += 1
                                    
                                    for edge in PH_loop_found:
                                        x0,y0 = edge[0]
                                        x1,y1 = edge[1]
                                        plt.plot([x0, x1], [y0, y1], color='black', lw=4, alpha=0.4)
                                    plt.scatter(ad_PH_locs[:,0], ad_PH_locs[:,1], color='red')
                                    plt.scatter(b_PH_locs[:,0], b_PH_locs[:,1], color='blue')
                                    plt.plot(xx, yy, color='black', lw=2, ls='--')
                                    plt.axis('off')

                                    from matplotlib.lines import Line2D

                                    custom_lines1 = [
                                                    Line2D([0], [0], marker='o', color='None', markerfacecolor='b', markersize=15),
                                                    Line2D([0], [0], marker='o', color='None', markerfacecolor='r', markersize=15),
                                                    ]
                                    legend1 = plt.legend(custom_lines1, [
                                                                    r'$\beta$-cells',
                                                                    r'$\alpha\delta$-cells'],
                                                         fontsize=14, loc='lower left'\
                                                                 , framealpha=0.4)


                                    custom_lines2 = [
                                                    Line2D([0], [0], color='black', ls='--', lw=2),
                                                    Line2D([0], [0], color='black', ls='-', lw=4, alpha=0.4),
                                                    ]
                                    legend2 = plt.legend(custom_lines2, [
                                                                    'Geom. mantle',
                                                                    'PH-mantle'], fontsize=14\
                                                                            , loc='upper left'\
                                                                            , framealpha=0.4)


                                    plt.gca().add_artist(legend1)
                                    plt.gca().add_artist(legend2)

                                    plt.savefig('example_geom_PH_b_mantle.pdf', dpi=600)

                                    plt.cla()
                                    plt.clf()

                                    exit()


                            close_scores.append(min_dist)
                            num_ad_in_mantle_PH[stage].append(len(geom_cells))
                            comp_lens_ad_in_mantle_PH[stage].append(this_comp_lens)
                            max_ad_in_mantle_PH[stage].append([n_total_cells, maxx_comp_size])
                        break
        
            if not flag_found:
                far_islets[fprefix] = 0
                n_unmatched += 1
    
    
    all_close += close_scores
    all_far += far_scores
    

    #all_results.write(fprefix + ", b-mantles:")
    #for entry in match:
    #    all_results.write(str(entry)+',')
    #all_results.write('\n')

pickle.dump(PH_islet_comps_sizes_dict, open('PH_mantles_comp_sizes_dev.p', 'wb'))
exit()

#############################################################################
#############################################################################
pickle.dump(far_islets, open('far_unmatched_islets_after_stochastic.p', 'wb'))

for stage in range(4):
    max_b_in_mantle[stage]     = np.array(max_b_in_mantle[stage])
    max_ad_in_mantle[stage]    = np.array(max_ad_in_mantle[stage])
                                                                    
    max_b_in_mantle_PH[stage]  = np.array(max_b_in_mantle_PH[stage])
    max_ad_in_mantle_PH[stage] = np.array(max_ad_in_mantle_PH[stage])


num_cells_in_mantle = [\
                         num_b_in_mantle\
                        ,num_ad_in_mantle\
                        ,num_b_in_mantle_PH\
                        ,num_ad_in_mantle_PH\
                        ]

comp_lens_in_mantle = [\
                         comp_lens_b_in_mantle\
                        ,comp_lens_ad_in_mantle\
                        ,comp_lens_b_in_mantle_PH\
                        ,comp_lens_ad_in_mantle_PH\
                        ]

    
max_components_in_mantle = [\
                             max_b_in_mantle\
                            ,max_ad_in_mantle\
                            ,max_b_in_mantle_PH\
                            ,max_ad_in_mantle_PH\
                            ]

pickle.dump(num_cells_in_mantle, open('num_cells_in_mantle_after_stochastic.p', 'wb'))
pickle.dump(comp_lens_in_mantle, open('comp_lens_in_mantle_after_stochastic.p', 'wb'))
pickle.dump(max_components_in_mantle, open('max_comp_size_in_mantle_after_stochastic.p', 'wb'))
#############################################################################
#############################################################################
    
#all_results.close() 

#all_scores = np.array(all_scores)

#unmatch_idxs = np.argwhere(all_scores == -1)
#match_idxs = np.argwhere(all_scores != -1)
#print('Unmatched = ', len(unmatch_idxs))
#print('Matched = ', len(match_idxs))

plt.hist(all_close, bins=50, color='blue', alpha=0.5, label='Nearby')

plt.hist(all_far, bins=50, color='red', alpha=0.5, label='Far')

plt.legend(fontsize=14)

plt.xlabel('distance of geometric mantle from closest PH mantle', fontsize=14)
plt.ylabel('count', fontsize=14)
plt.yscale('log', base = 2)
plt.title('Nearby '+ str(len(all_close)) + ', Far ' + str(len(all_far)) + ', Unmatched ' + str(n_unmatched))
#plt.show()
#exit()
plt.savefig('figures/geom_vs_PH_loops_after_stochastic_dev.pdf')

plt.cla()
plt.clf()


#mylabels = ['Nearby', 'Far', 'Unmatched']
#myexplode = [0, 0, 0.2]
#
#plt.pie([len(all_close), len(all_far), n_unmatched]\
#          , labels = mylabels\
#          , explode = myexplode\
#          ,autopct='%1.2f%%'\
#        )
##  ,pctdistance=1.25, labeldistance=1.4)
#
#plt.show()


