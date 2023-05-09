import numpy as np
import os
import pickle

ff = open('all_islets_info.csv', 'r')

all_prefix = []

ff.readline()
for line in ff:

    line = line.split(',')
    prefix = line[3]
    all_prefix.append(prefix)
    
ff.close()

skip_ff = open('skip_islets.csv', 'w')

count = 0
for fprefix in all_prefix:


    print(count, end='\r')
    count += 1

    v_file = 'VerticesEdges/' + fprefix + '.vertices'
    e_file = 'VerticesEdges/' + fprefix + '.edges'

    if not os.path.isfile(v_file) or not os.path.isfile(e_file):
        skip_ff.write(fprefix+'\n')
        continue

    
    cell_info = dict()
    
    vv = open(v_file, 'r')
    
    ad_cells_info = dict()
    b_cells_info = dict()
    
    ad_idx = 0
    b_idx = 0
    
    b_cells_list = []
    ad_cells_list = []
    
    for line in vv:
    
        line = line.strip('\n')
        line = line.split('\t')
    
        xx = float(line[0])
        yy = float(line[1])
    
        isl_num = int(line[2])
        typ = line[3]
    
        if typ == 'b':
            b_cells_list.append([xx, yy])
            b_cells_info[isl_num] = {'x':xx, 'y':yy, 'index':b_idx}
            b_idx += 1
        else:
            ad_cells_list.append([xx, yy])
            ad_cells_info[isl_num] = {'x':xx, 'y':yy, 'index':ad_idx}
            ad_idx += 1
    
    vv.close()

    # Save islet num to PH index dictionary
    ad_vert_geom_to_PH_file = 'VerticesEdges/' + fprefix + '_advertices_geom_to_PH.p'
    pickle.dump( ad_cells_info, open( ad_vert_geom_to_PH_file, "wb" ) )

    b_vert_geom_to_PH_file = 'VerticesEdges/' + fprefix + '_bvertices_geom_to_PH.p'
    pickle.dump( b_cells_info, open( b_vert_geom_to_PH_file, "wb" ) )
    
    ad_vert_PH_file = 'VerticesEdges/' + fprefix + '_advertices_PH.csv'

    
    vv = open(ad_vert_PH_file, 'w')
    for xx, yy in ad_cells_list:
        vv.write(str(xx) + ',' + str(yy) + '\n')
    vv.close()

    
    
    
    b_vert_PH_file = 'VerticesEdges/' + fprefix + '_bvertices_PH.csv'
    vv = open(b_vert_PH_file, 'w')
    for xx, yy in b_cells_list:
        vv.write(str(xx) + ',' + str(yy) + '\n')
    vv.close()
    
    
    e_file = 'VerticesEdges/' + fprefix + '.edges'
    ee = open(e_file, 'r')
    
    ad_edges_list = []
    ad_edges_dist_list = []
    
    b_edges_list = []
    b_edges_dist_list = []
    
    for line in ee:
        line = line.strip('\n')
        line = line.split(',')
        line = [int(x) for x in line]
    
        c1 = line[0]
        c2 = line[1]
    
        if c1 in ad_cells_info:
    
            c1_idx = ad_cells_info[c1]['index']
            c2_idx = ad_cells_info[c2]['index']
            ad_edges_list.append([c1_idx, c2_idx])
    
            xx1 = ad_cells_info[c1]['x']
            yy1 = ad_cells_info[c1]['y']
    
            xx2 = ad_cells_info[c2]['x']
            yy2 = ad_cells_info[c2]['y']
    
            ddist = (xx1-xx2)**2 + (yy1-yy2)**2
            ad_edges_dist_list.append(ddist)
    
            #ad_edges.write(str(c1_idx) + ',' + str(c2_idx) + '\n')
    
        else:
    
            c1_idx = b_cells_info[c1]['index']
            c2_idx = b_cells_info[c2]['index']
            b_edges_list.append([c1_idx, c2_idx])
    
            xx1 = b_cells_info[c1]['x']
            yy1 = b_cells_info[c1]['y']
                  
            xx2 = b_cells_info[c2]['x']
            yy2 = b_cells_info[c2]['y']
    
            ddist = (xx1-xx2)**2 + (yy1-yy2)**2
            b_edges_dist_list.append(ddist)
    
            #_edges.write(str(c1_idx) + ',' + str(c2_idx) + '\n')
    
    
    ad_edges_list = np.array(ad_edges_list)
    b_edges_list = np.array(b_edges_list)
    
    
    ad_edges_dist_list = np.array(ad_edges_dist_list)
    b_edges_dist_list = np.array(b_edges_dist_list)
    
    ad_sort_args = np.argsort(ad_edges_dist_list)
    b_sort_args = np.argsort(b_edges_dist_list)
    
    ad_edges_dist_list = ad_edges_dist_list[ad_sort_args]
    b_edges_dist_list = b_edges_dist_list[b_sort_args]
    
    ad_edges_list = ad_edges_list[ad_sort_args]
    b_edges_list = b_edges_list[b_sort_args]
    
    
    ad_edge_file = 'VerticesEdges/' + fprefix + '_adedges_PH.csv'
    b_edge_file = 'VerticesEdges/' + fprefix + '_bedges_PH.csv'
    
    #print(ad_edges_list)
    
    np.savetxt(ad_edge_file, ad_edges_list, delimiter=',', fmt='%d')
    np.savetxt(b_edge_file, b_edges_list, delimiter=',', fmt='%d')
    
    
skip_ff.close() 
