import pickle
import os
import glob
import numpy as np
import math
import matplotlib.pyplot as plt
import subprocess
from sklearn.neighbors import KDTree

islets = pickle.load(open('far_unmatched_islets.p', 'rb'))


dirr = 'stochastic'

thresh_dirr = 'Thresholds'
data_dirr = 'DevData/all_islets'

if not os.path.isdir(dirr):
    os.mkdir(dirr)

n_pert = 10

fcount = 0

for fprefix in islets:

    fcount += 1

    islet_dirr = dirr + '/' + fprefix
    if not os.path.isdir(islet_dirr):
        os.mkdir(islet_dirr)

    thresh_data_file = thresh_dirr + '/' + fprefix + '.smooth2.thresh.minBetPeak2And3.dat'
    
    thresh_file = 'Thresholds/'+fprefix + '.smooth2.thresh.minBetPeak2And3.dat'
    threshs_data = open(thresh_file, 'r')
    threshs = threshs_data.readline()
    threshs = threshs.strip('\n')
    threshs = threshs.split('\t')
    ad_thresh = float(threshs[1])
    b_thresh = float(threshs[2])
    adm_thresh = float(threshs[3])
    threshs_data.close()


    islet_file = data_dirr + '/' + fprefix + '.tsv'
    islet_data = open(islet_file, 'r')
    locs = []
    for line in islet_data:
        line = line.split('\t')
        locs.append([float(line[2]), float(line[3])])
    islet_data.close()

    locs = np.array(locs)
    n_pts = locs.shape[0]

    tree = KDTree(locs) 
    ddist,i = tree.query(locs, k=2, return_distance=True)

    ddist = ddist[:,1]

    print(np.amin(ddist), np.amax(ddist))

    exit()




    max_r = 0.1


    pert_dirr = islet_dirr + '/Perturbations'
    if not os.path.isdir(pert_dirr):
        os.mkdir(pert_dirr)

    for pert in range(n_pert):
        print(fcount, pert)

        this_pert_dirr = pert_dirr + '/Perturbation'+str(pert)
        if not os.path.isdir(this_pert_dirr):
            os.mkdir(this_pert_dirr)
        
        flag = 1
        while flag:
            this_rand_perts = np.random.uniform(size=(n_pts, 2))
            delta_x = this_rand_perts[:,0]*pert_r*np.cos(2*math.pi*this_rand_perts[:,1])
            delta_y = this_rand_perts[:,0]*pert_r*np.sin(2*math.pi*this_rand_perts[:,1])

            new_locs = np.copy(locs)

            new_locs[:,0] += delta_x
            new_locs[:,1] += delta_y

            new_pert_file = this_pert_dirr + '/locs.tsv'
            new_pert_data = open(new_pert_file, 'w')

            islet_file = data_dirr + '/' + fprefix + '.tsv'
            islet_data = open(islet_file, 'r')
            locs = []
            idx = 0
            for line in islet_data:
                line = line.split('\t')
                locs.append([float(line[2]), float(line[3])])
                new_line = line[0] + '\t' + line[1] + '\t'\
                        + str(new_locs[idx, 0]) + '\t' + str(new_locs[idx, 1]) + '\t'\
                        + line[4] + '\t' + line[5]
                new_pert_data.write(new_line)

                idx += 1

            islet_data.close()
            new_pert_data.close()


            subprocess.run(["./SimplexCode_12519_v3/Simplex_VersionForPaper"\
                            , new_pert_file\
                            , "321Area"\
                            , "B"\
                            , "-grdir"\
                            , this_pert_dirr
                            , "-fout"
                            , "locs"
                            , "-runGr"
                            , "y"
                            ])


            new_thresh_file = this_pert_dirr + '/locs.smooth2.thresh.minBetPeak2And3.dat'

            new_threshs = open(new_thresh_file, 'r')

            line = new_threshs.readline()
            line = line.split('\t')
            new_ad_thresh = float(line[1])
            new_b_thresh = float(line[2])
            new_adm_thresh = float(line[3])
            new_threshs.close()

            if (round(ad_thresh,1) == round(new_ad_thresh,1))\
                and (round(b_thresh,1) == round(new_b_thresh,1))\
                and (round(adm_thresh,1) == round(new_adm_thresh,1)):
                    flag = 0
                    #print(ad_thresh, b_thresh, adm_thresh)
                    #print(new_ad_thresh, new_b_thresh, new_adm_thresh)
                    #print('Perturbation found!')
                    #input('w')



        #exit()





        
        
        






    








