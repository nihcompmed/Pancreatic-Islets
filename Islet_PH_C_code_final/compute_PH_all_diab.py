import numpy as np
import os
import gc
import pydory as dory


# Diabetic data

ff = open('all_islets_info_diabetic.csv', 'r')

# Dory pars
filetype = 1
low_thresh = 0
up_thresh = -1
threads = 4
dim = 1

compute_cycles = 0
reduce_cyc_lengths = 0
birth_thresh = 0
suppress_output = 1
hom_ws = 1000
cohom_ws = 100



ff.readline()
for line in ff:

    info = line.split(',')

    fprefix = info[3]
                               
    print('Processing', fprefix)

    v_file = 'VerticesEdges/' + fprefix + '.vertices'

    if not os.path.isfile(v_file):
        continue

    vv = open(v_file, 'r')

    PH_ad_file = 'PD_results/' + fprefix + '_ad.csv'
    PH_b_file = 'PD_results/' + fprefix + '_b.csv'

    PH_ad = open(PH_ad_file, 'w')
    PH_b = open(PH_b_file, 'w')

    for line in vv:
        line = line.strip('\n')
        line = line.split('\t')
        xx = line[0]
        yy = line[1]
        typ = line[-1]

        if typ == 'b':
            PH_b.write(xx + ',' + yy + '\n')
        else:
            PH_ad.write(xx + ',' + yy + '\n')

    

    PH_b.close()
    PH_ad.close()

    source = PH_ad_file
    target = 'PD_results/' + fprefix + '_ad_'

    #if not os.path.isfile(target + 'H1_pers_data.txt'):


    dory.compute_PH(source\
                , low_thresh\
                , up_thresh\
                , filetype\
                , threads\
                , target\
                , dim\
                , compute_cycles\
                , reduce_cyc_lengths\
                , birth_thresh\
                , suppress_output\
                , hom_ws\
                , cohom_ws)
    source = PH_b_file
    target = 'PD_results/' + fprefix + '_b_'

    #if not os.path.isfile(target + 'H1_pers_data.txt'):
    dory.compute_PH(source\
                    , low_thresh\
                    , up_thresh\
                    , filetype\
                    , threads\
                    , target\
                    , dim\
                    , compute_cycles\
                    , reduce_cyc_lengths\
                    , birth_thresh\
                    , suppress_output\
                    , hom_ws\
                    , cohom_ws)



ff.close()



