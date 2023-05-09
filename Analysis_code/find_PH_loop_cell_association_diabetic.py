import numpy as np
import islet_helper as ih
import  math
import pickle
import os

def get_loops(ffile, cell_info, edges_arr):

    undead = open(ffile, 'r')
    all_loops = []
    for line in undead:
        line = line.strip('\n')
        line = line.split(',')
        line = line[:-1]
        line = [int(x) for x in line]
    
        loop_coords = []
    
        nn = 0
        while (nn < len(line)):
            
            this_edge = line[nn]
            c1 = edges_arr[this_edge][0]
            c2 = edges_arr[this_edge][1]
    
            x1 = cell_info[c1]['x']
            y1 = cell_info[c1]['y']
    
            x2 = cell_info[c2]['x']
            y2 = cell_info[c2]['y']
    
            loop_coords.append([[x1, y1], [x2, y2]])
    
            nn += 1
    
        all_loops.append(loop_coords)
    
    undead.close()

    return all_loops

def get_association(loops, cells):

    info = dict()

    for cell in cells:
    
        this_cell = [cells[cell]['x'], cells[cell]['y']]
    
        min_dist = math.inf
        loop_id = []
    
        for idx, loop_coords in enumerate(loops):
    
            flag = ih.geom_loop_has_cell(loop_coords, this_cell)
    
            if not flag:
                continue
    
            ddist = ih.max_dist_loop_cell(loop_coords, this_cell)
    
            if ddist < min_dist:
                min_dist = ddist
                loop_id = [idx]
            elif ddist == min_dist:
                loop_id.append(idx)
    
    
        loop_set = frozenset(loop_id)
    
        if loop_set not in info:
            info[loop_set] = [cell]
        else:
            info[loop_set].append(cell)

    return info

ff = open('all_islets_info_diabetic.csv', 'r')

all_prefix = []

ff.readline()
for line in ff:

    line = line.split(',')
    prefix = line[3]
    all_prefix.append(prefix)
    
ff.close()

count = 0
for fprefix in all_prefix:

    print(count, end='\r')
    count += 1

    v_file = 'VerticesEdges/' + fprefix + '.vertices'
    e_file = 'VerticesEdges/' + fprefix + '.edges'

    if not os.path.isfile(v_file) or not os.path.isfile(e_file):
        continue

    b_v_file = 'VerticesEdges/' + fprefix + '_bvertices_PH.csv'
    ad_v_file = 'VerticesEdges/' + fprefix + '_advertices_PH.csv'
    
    b_e_file = 'VerticesEdges/' + fprefix + '_bedges_PH.csv'
    ad_e_file = 'VerticesEdges/' + fprefix + '_adedges_PH.csv'
    
    ad_loop_file = 'Loops/' + fprefix + '_adundead_PH.csv'
    b_loop_file = 'Loops/' + fprefix + '_bundead_PH.csv'
    
    ad_edges_arr = np.loadtxt(ad_e_file, delimiter=',', dtype=int)
    b_edges_arr = np.loadtxt(b_e_file, delimiter=',', dtype=int)
    
    b_cell_info = dict()
    vv = open(b_v_file, 'r')
    idx = 0
    for line in vv:
        line = line.strip('\n')
        line = line.split(',')
    
        xx = float(line[0])
        yy = float(line[1])
    
        b_cell_info[idx] = {'x':xx, 'y':yy}
        idx += 1
    
        #all_xx.append(xx)
        #all_yy.append(yy)
        #colors.append('blue')
    
    vv.close()
    
    ad_cell_info = dict()
    vv = open(ad_v_file, 'r')
    idx = 0
    for line in vv:
        line = line.strip('\n')
        line = line.split(',')
    
        xx = float(line[0])
        yy = float(line[1])
    
        ad_cell_info[idx] = {'x':xx, 'y':yy}
        idx += 1
    
        #all_xx.append(xx)
        #all_yy.append(yy)
        #colors.append('red')
    
    vv.close()
    
    
    # Find b-cell to ad-mantle association for PH loops
    
    ad_PH_loops = get_loops(ad_loop_file, ad_cell_info, ad_edges_arr)
    ad_info = get_association(ad_PH_loops, b_cell_info)
    
    b_PH_loops = get_loops(b_loop_file, b_cell_info, b_edges_arr)
    b_info = get_association(b_PH_loops, ad_cell_info)
    
    
    ad_PH_file = 'Loops/' + fprefix + '_admantles_PH.p'
    pickle.dump( [ad_PH_loops, ad_info], open( ad_PH_file, "wb" ) )
    
    b_PH_file = 'Loops/' + fprefix + '_bmantles_PH.p'
    pickle.dump( [b_PH_loops, b_info], open( b_PH_file, "wb" ) )
    
    
    
    
