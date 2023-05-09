import numpy as np
import os
import pickle
import islet_helper as ih



ff = open('all_islets_info.csv', 'r')

all_info = dict()

all_prefix = []

ff.readline()
for line in ff:

    info = line.split(',')

    this_stage   = int(info[0])
    this_subject = info[1]
    this_islet   = int(info[2])

    fprefix = info[3]
                               
    this_area    = float(info[4])
    this_n_ad    = int(info[5])
    this_n_b     = int(info[6])

    print('Processing', fprefix)

    if this_stage not in all_info:
        all_info[this_stage] = dict()

    if this_subject not in all_info[this_stage]:
        all_info[this_stage][this_subject] = dict()

    all_info[this_stage][this_subject][this_islet] = dict()


    all_info[this_stage][this_subject][this_islet]['n_ad'] = this_n_ad
    all_info[this_stage][this_subject][this_islet]['n_b'] = this_n_b

    all_info[this_stage][this_subject][this_islet]['area'] = this_area

    e_file = 'VerticesEdges/' + fprefix + '.edges'

    if not os.path.isfile(e_file):
        all_info[this_stage][this_subject][this_islet]['edges'] = 0
        all_info[this_stage][this_subject][this_islet]['gor'] = 0

        all_info[this_stage][this_subject][this_islet]['b_in_ad_geom'] = 0
        all_info[this_stage][this_subject][this_islet]['ad_in_b_geom'] = 0

        all_info[this_stage][this_subject][this_islet]['b_in_ad_PH'] = 0
        all_info[this_stage][this_subject][this_islet]['ad_in_b_PH'] = 0
        continue

    all_info[this_stage][this_subject][this_islet]['edges'] = 1

    # Was G(r) of this islet computed?
    gr_file = 'Thresholds/' + fprefix + '.smooth2.thresh.minBetPeak2and3.dat'

    if not os.path.isfile(gr_file):
        all_info[this_stage][this_subject][this_islet]['gor'] = 0

        all_info[this_stage][this_subject][this_islet]['b_in_ad_geom'] = 0
        all_info[this_stage][this_subject][this_islet]['ad_in_b_geom'] = 0

        all_info[this_stage][this_subject][this_islet]['b_in_ad_PH'] = 0
        all_info[this_stage][this_subject][this_islet]['ad_in_b_PH'] = 0
        # G(r) is not computed, continue
        continue
    
    all_info[this_stage][this_subject][this_islet]['gor'] = 1

    # Load geom to PH map
    # Map is of form dict[isl_num] = {'x', 'y', 'index'}
    ad_vert_geom_to_PH_file = 'VerticesEdges/' + fprefix + '_advertices_geom_to_PH.p'
    ad_map = pickle.load( open( ad_vert_geom_to_PH_file, "rb" ) )
    ad_inv_map = dict()
    for ad_cell in ad_map:
        isl_num = ad_cell
        PH_index = ad_map[ad_cell]['index']
        ad_inv_map[PH_index] = isl_num

    b_vert_geom_to_PH_file = 'VerticesEdges/' + fprefix + '_bvertices_geom_to_PH.p'
    b_map = pickle.load( open( b_vert_geom_to_PH_file, "rb" ) )
    b_inv_map = dict()
    for b_cell in b_map:
        isl_num = b_cell
        PH_index = b_map[b_cell]['index']
        b_inv_map[PH_index] = isl_num
         

    

    # Get geometric loops
    ad_geom_loops = 'Loops/' + fprefix + '.admantles'
    if not os.path.isfile(ad_geom_loops):
        all_info[this_stage][this_subject][this_islet]['b_in_ad_geom'] = 0
    else:
        all_info[this_stage][this_subject][this_islet]['b_in_ad_geom'] = ih.get_geom_loops(ad_geom_loops)
    
    b_geom_loops = 'Loops/' + fprefix + '.bmantles'
    if not os.path.isfile(b_geom_loops):
        all_info[this_stage][this_subject][this_islet]['ad_in_b_geom'] = 0
    else:

        all_info[this_stage][this_subject][this_islet]['ad_in_b_geom'] = ih.get_geom_loops(b_geom_loops)


    # Get PH loops
    ad_PH_loops_file = 'Loops/' + fprefix + '_admantles_PH.p'

    if not os.path.isfile(ad_PH_loops_file):
        all_info[this_stage][this_subject][this_islet]['b_in_ad_PH'] = 0
    else:
        all_info[this_stage][this_subject][this_islet]['b_in_ad_PH'] =\
                                                    ih.get_PH_loops(ad_PH_loops_file, b_inv_map)

    b_PH_loops_file = 'Loops/' + fprefix + '_bmantles_PH.p'
    if not os.path.isfile(b_PH_loops_file):
        all_info[this_stage][this_subject][this_islet]['ad_in_b_PH'] = 0
    else:
        all_info[this_stage][this_subject][this_islet]['ad_in_b_PH'] =\
                                                    ih.get_PH_loops(b_PH_loops_file, ad_inv_map)


    #print(all_info[this_stage][this_subject][this_islet])

all_info_file = 'development_all_info.p'

pickle.dump( all_info, open( all_info_file, "wb" ) )





