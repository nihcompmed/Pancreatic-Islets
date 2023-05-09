import matplotlib.pyplot as plt
import pickle
import numpy as np
import matplotlib as mpl
import scipy
import itertools as it
import seaborn as sns
import pandas as pd
import islet_helper as ih
import networkx as nx
import os


colors = [\
        'tab:blue'\
        ,'tab:orange'\
        ,'tab:olive'\
        ,'tab:green'\
          ]

dev_stages = [\
            'Gestation/'\
            ,'1-35weeks/'\
            ,'12-24months/'
            ,'28MandAfter/'\
                    ]

ff = open('all_islets_info.csv', 'r')
all_prefix = []

# For now taking islet sizes that are given
# Otherwise, estimate area from the tightest rectangle around the point-cloud
all_areas = []

all_b_fraction = []

ff.readline()
for line in ff:

    line = line.split(',')
    prefix = line[3]
    all_prefix.append(prefix)

    area = float(line[-3])

    all_areas.append(area)

    
ff.close()

count = 0


b_S_in_mantle_count       = dict() 
b_NS_in_mantle_count      = dict() 

b_S_not_in_mantle_count   = dict() 
b_NS_not_in_mantle_count  = dict() 

ad_S_in_mantle_count      = dict() 
ad_NS_in_mantle_count     = dict() 

ad_S_not_in_mantle_count  = dict() 
ad_NS_not_in_mantle_count = dict() 

n_islets = dict()

# area, total_cells, b_cell_fraction
ALL_DATA_B_ALL_IN_MANTLE = dict()
ALL_DATA_B_NS_IN_MANTLE = dict()

ALL_DATA_B_COMPS = dict()
ALL_DATA_ISLETS = dict()

islet_comps_sizes_dict = dict()
#islet_comps_sizes_not_in_mantle_dict = dict()


for stage in range(4):

    b_S_in_mantle_count[stage]       = 0 
    b_NS_in_mantle_count[stage]      = 0 
    
    b_S_not_in_mantle_count[stage]   = 0 
    b_NS_not_in_mantle_count[stage]  = 0 
    
    ad_S_in_mantle_count[stage]      = 0 
    ad_NS_in_mantle_count[stage]     = 0 
    
    ad_S_not_in_mantle_count[stage]  = 0 
    ad_NS_not_in_mantle_count[stage] = 0 

    n_islets[stage] = 0

    ALL_DATA_B_ALL_IN_MANTLE[stage] = []
    ALL_DATA_B_NS_IN_MANTLE[stage] = []

    ALL_DATA_ISLETS[stage] = []
    ALL_DATA_B_COMPS[stage] = []

    islet_comps_sizes_dict[stage] = []
    #islet_comps_sizes_not_in_mantle_dict[stage] = []

max_total_cells = 0

for count, fprefix in enumerate(all_prefix):

    if count % 1000 == 0:
        print(count)


    v_file = 'VerticesEdges/' + fprefix + '.vertices'
    e_file = 'VerticesEdges/' + fprefix + '.edges'

    if not os.path.isfile(v_file) or not os.path.isfile(e_file):
        continue

    THIS_AREA = all_areas[count]

    stage = int(fprefix.split('_')[0][-1])

    n_islets[stage] += 1

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

    all_b_fraction.append(this_b_fraction)

    this_islet_characterization = (this_b_fraction, n_total_cells)

    max_total_cells = max(max_total_cells, n_total_cells)

    islet_comps_sizes_dict[stage].append({'feature':this_islet_characterization\
                                        , 'b_comp_lens':[]\
                                        , 'b_comp_lens_inside_mantle':[]\
                                        , 'b_comp_lens_not_inside_mantle':[]\
                                        , 'ad_comp_lens':[]\
                                        , 'ad_comp_lens_inside_mantle':[]\
                                        , 'ad_comp_lens_not_inside_mantle':[]\
                                        })

    #islet_comps_sizes_not_in_mantle_dict[stage].append({'feature':this_islet_characterization\
    #                                    , 'comp_lens_not_inside_mantle':[]\
    #                                    })

    edges = np.loadtxt(e_file, delimiter=',', dtype=int)
    if edges.ndim == 1:
        edges = np.reshape(edges, (1, 2))
    G = nx.Graph()
    G.add_nodes_from(v_islnum)
    for edge in edges:
        G.add_edge(edge[0], edge[1])

    b_graph = G.subgraph(all_b_cells)
    ad_graph = G.subgraph(all_ad_cells)

    ALL_DATA_ISLETS[stage].append([this_b_fraction, n_total_cells, THIS_AREA])

    n_b_comps = len(list(nx.connected_components(b_graph)))
    for ii in range(n_b_comps):
        ALL_DATA_B_COMPS[stage].append([this_b_fraction, n_total_cells, THIS_AREA])


    b_comps = nx.connected_components(b_graph)
    ad_comps = nx.connected_components(ad_graph)

    for b_comp in b_comps:
        islet_comps_sizes_dict[stage][-1]['b_comp_lens'].append(len(b_comp))

    for ad_comp in ad_comps:
        islet_comps_sizes_dict[stage][-1]['ad_comp_lens'].append(len(ad_comp))


    b_comps_in_mantle = dict()
    ad_comps_in_mantle = dict()

    for comp in b_comps:
        b_comps_in_mantle[frozenset(comp)] = 0

    for comp in ad_comps:
        ad_comps_in_mantle[frozenset(comp)] = 0

    # Geometric loops
    # Each line is of the form Islcell1, Islcell2, ...:x1,y1,x2,y2,.... of the loop around the cells.
    # (n,n+1), (n+1, n+2) ... are edges of the loop
    ad_geom_loops = 'Loops/' + fprefix + '.admantles'
    b_geom_loops = 'Loops/' + fprefix + '.bmantles'

    bcells_in_mantle = []


    if os.path.isfile(ad_geom_loops):
        vv = open(ad_geom_loops, 'r')
        for line in vv:
            line = line.strip('\n')
            line = line.split(':')

            bcells = line[0].split(',')
            bcells = bcells[:-1]

            bcells2 = []

            for x in bcells:
                y = int(x)
                if y < 100000:
                    bcells2.append(y)


            bcells_in_mantle += [int(x) for x in bcells2]

        
            #bcells2 has isletnumbers of b-cells inside a geom loop
            this_b_graph = G.subgraph([int(x) for x in bcells2])

            this_b_components = nx.connected_components(this_b_graph)

            for this_comp in this_b_components:

                islet_comps_sizes_dict[stage][-1]['b_comp_lens_inside_mantle'].append(len(this_comp))

                b_comps_in_mantle[frozenset(this_comp)] = 1
                ALL_DATA_B_ALL_IN_MANTLE[stage].append([this_b_fraction, n_total_cells, THIS_AREA])
                if len(this_comp) > 1:
                    ALL_DATA_B_NS_IN_MANTLE[stage].append([this_b_fraction, n_total_cells, THIS_AREA])

    n_b_comps_in_mantle = len(b_comps_in_mantle)

    bcells_in_mantle_set = frozenset(bcells_in_mantle)
    all_b_cells_set = frozenset(all_b_cells)

    bcells_not_in_mantle_set = all_b_cells_set - bcells_in_mantle_set

    bgraph_not_in_mantle = G.subgraph(list(bcells_not_in_mantle_set))
    bcomps_not_in_mantle = nx.connected_components(bgraph_not_in_mantle)

    for this_comp in bcomps_not_in_mantle:
        islet_comps_sizes_dict[stage][-1]['b_comp_lens_not_inside_mantle'].append(len(this_comp))


    adcells_in_mantle = []

    if os.path.isfile(b_geom_loops):
        vv = open(b_geom_loops, 'r')
        for line in vv:
            line = line.strip('\n')
            line = line.split(':')

            adcells = line[0].split(',')
            adcells = adcells[:-1]

            adcells2 = []
        
            for x in adcells:
                y = int(x)
                if y < 100000:
                    adcells2.append(y)

            adcells_in_mantle += [int(x) for x in adcells2]

            #adcells2 has isletnumbers of b-cells inside a geom loop
            this_ad_graph = G.subgraph([int(x) for x in adcells2])

            this_ad_components = nx.connected_components(this_ad_graph)

            for this_comp in this_ad_components:
                ad_comps_in_mantle[frozenset(this_comp)] = 1
                islet_comps_sizes_dict[stage][-1]['ad_comp_lens_inside_mantle'].append(len(this_comp))


    adcells_in_mantle_set = frozenset(adcells_in_mantle)
    all_ad_cells_set = frozenset(all_ad_cells)

    adcells_not_in_mantle_set = all_ad_cells_set - adcells_in_mantle_set

    adgraph_not_in_mantle = G.subgraph(list(adcells_not_in_mantle_set))
    adcomps_not_in_mantle = nx.connected_components(adgraph_not_in_mantle)

    for this_comp in adcomps_not_in_mantle:
        islet_comps_sizes_dict[stage][-1]['ad_comp_lens_not_inside_mantle'].append(len(this_comp))





    for comp in b_comps_in_mantle:

        if b_comps_in_mantle[comp] == 0:
            if len(comp) == 1:
                b_S_not_in_mantle_count[stage] += 1
            else:
                b_NS_not_in_mantle_count[stage] += 1
        else:
            if len(comp) == 1:
                b_S_in_mantle_count[stage] += 1
            else:
                b_NS_in_mantle_count[stage] += 1

    for comp in ad_comps_in_mantle:

        if ad_comps_in_mantle[comp] == 0:
            if len(comp) == 1:
                ad_S_not_in_mantle_count[stage] += 1
            else:
                ad_NS_not_in_mantle_count[stage] += 1
        else:
            if len(comp) == 1:
                ad_S_in_mantle_count[stage] += 1
            else:
                ad_NS_in_mantle_count[stage] += 1

print(max_total_cells)

pickle.dump(ALL_DATA_B_ALL_IN_MANTLE, open('B_ALL_COMPS_IN_MANTLE_DICT.p', 'wb'))
pickle.dump(ALL_DATA_B_NS_IN_MANTLE, open('B_NS_COMPS_IN_MANTLE_DICT.p', 'wb'))
pickle.dump(ALL_DATA_B_COMPS, open('B_COMPS_DICT.p', 'wb'))
pickle.dump(ALL_DATA_ISLETS, open('ISLETS_DICT.p', 'wb'))

pickle.dump(islet_comps_sizes_dict, open('B_COMPS_SIZES_IN_MANTLE_DICT.p', 'wb'))
#pickle.dump(islet_comps_sizes_not_in_mantle_dict, open('B_COMPS_SIZES_NOT_IN_MANTLE_DICT.p', 'wb'))

exit()


mylabels = ["S in mantle", "S not in mantle", "NS in mantle", "NS not in mantle"]

for stage in range(4):

    plt.pie([b_S_in_mantle_count[stage]\
          , b_S_not_in_mantle_count[stage]\
          , b_NS_in_mantle_count[stage]\
          , b_NS_not_in_mantle_count[stage]]\
          , labels = mylabels)
    
    plt.show()

    plt.cla()
    plt.clf()

    plt.pie([ad_S_in_mantle_count[stage]\
          , ad_S_not_in_mantle_count[stage]\
          , ad_NS_in_mantle_count[stage]\
          , ad_NS_not_in_mantle_count[stage]]\
          , labels = mylabels)
    
    plt.show()

    plt.cla()
    plt.clf()


exit()

