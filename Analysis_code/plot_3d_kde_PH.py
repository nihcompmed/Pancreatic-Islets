#import plotly.graph_objects as go
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
import math
from scipy import stats, special
from scipy.stats import entropy
from sklearn.cluster import MeanShift, OPTICS

def transform(arr):
    arr = np.log1p(arr)
    return arr

def get_2d_kde_plot(data, rranges, grid_reso, aspect, save_fname_prefix):


    xmin, xmax = rranges[0]
    ymin, ymax = rranges[1]


    # heat map

    plt.hist2d(data.T[:,0], data.T[:,1]\
            , range=rranges)
    plt.colorbar(label='Count')
    plt.xlabel('$\\beta$ cell fraction', fontsize=16)
    plt.ylabel(r'Total cells in islet (ln(1+val))', fontsize=16)
    plt.savefig(save_fname_prefix+'_heatmap.pdf', dpi=600)
    plt.cla()
    plt.clf()

    kernel = stats.gaussian_kde(data)

    X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
    positions = np.vstack([X.ravel(), Y.ravel()])

    kernel_estimates = kernel(positions)

    maxx = np.argmax(kernel_estimates)
    maxx_at = positions.T[maxx]


    Z = np.reshape(kernel_estimates.T, X.shape)


    maxx = np.unravel_index(np.argmax(Z, axis=None), Z.shape)

    plt.imshow(np.rot90(Z), cmap=plt.cm.gnuplot,
              extent=[xmin, xmax, ymin, ymax]\
                      , aspect=aspect)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$\beta$-fraction', fontsize=16)
    plt.ylabel(r'total cells (ln(1+val))', fontsize=16)

    # Mark max kde 
    plt.scatter(maxx_at[0], maxx_at[1], color='red', marker='*')

    plt.savefig(save_fname_prefix+'_kde.pdf', dpi=600)

    plt.cla()
    plt.clf()

    return [kernel, maxx_at]

def get_stat_ranges(all_data, stat_kind, global_stat_ranges):

    for stage in range(4):

        info = all_data[stage]
    
        for idx, islet in enumerate(info):

            data = islet['comp_lens']

            if data == []:
                islet_stat = 0
                # If you want to ignore islets without components in mantles
                continue
            else:
                data = np.array(data)
                #if np.amax(data) < 2:
                #    continue
                if stat_kind == 'mean':
                    islet_stat = np.mean(data)
                elif stat_kind == 'max':
                    islet_stat = np.amax(data)
                elif stat_kind == 'count_NS':
                    islet_stat = len(np.argwhere(data > 1).flatten())
                elif stat_kind == 'count_S':
                    islet_stat = len(np.argwhere(data == 1).flatten())

            global_stat_ranges[stat_kind]['min'] = min(global_stat_ranges[stat_kind]['min']\
                                                        , islet_stat)

            global_stat_ranges[stat_kind]['max'] = max(global_stat_ranges[stat_kind]['max']\
                                                        , islet_stat)



def get_kldiv(kernel_info_file, grid_reso):

    info = pickle.load(open(kernel_info_file, 'rb'))

    all_kernels, rranges = info

    print(all_kernels)

    xmin, xmax = rranges[0]
    ymin, ymax = rranges[1]

    X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
    positions = np.vstack([X.ravel(), Y.ravel()])

    kernels = []
    maxx = []

    for stage in all_kernels:
        kernels.append(all_kernels[stage][0])
        maxx.append(all_kernels[stage][1])

    # pairwise
    for s1_idx, s2_idx in it.combinations(list(range(4)), 2):

        k1 = kernels[s1_idx]
        k2 = kernels[s2_idx]

        data1 = k1(positions)
        data2 = k2(positions)

        entr = entropy(data1, data2)

        print(f'kl div between stage {s1_idx} and stage {s2_idx} is {entr}')







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

stage_labels = [\
                'Stage0'\
                ,'Stage1'\
                ,'Stage2'\
                ,'Stage3'\
                ]


b_comps_in_mantle_info = pickle.load(open('PH_mantles_comp_sizes_dev.p', 'rb'))


###################
### pipeline #######
###################

pipeline_get_global_rranges = 1

pipeline_get_kde_all_islets = 0

pipeline_get_kde_islets_bcomps_inmantle = 0

pipeline_get_kde_islets_NSbcomps_inmantle = 1
pipeline_get_kde_islets_NSadcomps_inmantle = 1

pipeline_get_kde_islets_with_NS_bcomps = 0
pipeline_get_kde_islets_with_NO_NS_bcomps = 0

pipeline_get_kde_islets_with_NO_NS_adcomps = 0

# 3d kde
pipeline_compare_ns_comps_inside_and_out = 0




###################

grid_reso = 100j

tot_isl_counts = dict()
NSb_in_ad_isl_counts = dict()
NSad_in_b_isl_counts = dict()

# 1. Count the number of islets ( at least 5 b and 5 ad cells) for each stage
for stage in range(4):

    n_islets = len(b_comps_in_mantle_info[stage])
    print(f'Stage {stage} has {n_islets} islets')
    tot_isl_counts[stage] = n_islets


# 1a. Get global ranges for islet characterization
if pipeline_get_global_rranges:
    global_x_min = math.inf
    global_x_max = 0
    
    global_y_min = math.inf
    global_y_max = 0

    for stage in range(4):
    
        for islet in b_comps_in_mantle_info[stage]:
    
            this_islet_characterization = islet['feature']
    
            global_x_min = min(global_x_min, this_islet_characterization[0])
            global_x_max = max(global_x_max, this_islet_characterization[0])
    
    
            global_y_min = min(global_y_min, this_islet_characterization[1])
            global_y_max = max(global_y_max, this_islet_characterization[1])


rranges = np.array([[global_x_min, global_x_max], [global_y_min, global_y_max]])
# transform
t_rranges = np.copy(rranges)
t_rranges[1,:] = transform(t_rranges[1,:])
print(t_rranges)
exit()

# 2. Heat map of islets that have at least one b-cell NS component in a mantle
if pipeline_get_kde_islets_NSbcomps_inmantle:

    character_info = dict()

    for stage in range(4):
    
        character_info[stage] = []
    
        count = 0
    
        for islet in b_comps_in_mantle_info[stage]:


            # at least one bcomp in mantle
            if not len(islet['b_comp_lens_inside_PHmantle']):
                continue

            comps_inside = np.array(islet['b_comp_lens_inside_PHmantle'])

            NS_comps_inside = np.argwhere(comps_inside > 1).flatten()

            if len(NS_comps_inside) == 0:
                continue

            count += 1
    
            this_islet_characterization = islet['feature']
    
            character_info[stage].append(this_islet_characterization)

        print(f'Stage {stage} has {count} islets with at least NS b comp in a mantle')

        NSb_in_ad_isl_counts[stage] = count

    # 2b. Plot Kde maps and save kernels
    
    all_islet_kernels = dict()

    for stage in range(4):
    
        this_info = np.array(character_info[stage]).T
    
        # transform
        this_info[1,:] = transform(this_info[1,:])
    
        save_fname_prefix = 'stage'+str(stage)+'_atleast_one_NSbcomp_in_PHmantle'
    
        all_islet_kernels[stage] = get_2d_kde_plot(this_info, t_rranges, grid_reso, 0.1, save_fname_prefix)

    
    
    save_kernel = 'islets_atleast_one_NSbcomp_in_PHmantle.p'
    
    pickle.dump([all_islet_kernels, t_rranges], open(save_kernel, 'wb'))

    print('Computing kl div of islets with at least one NS bcomp in PH mantle')
    
    get_kldiv(save_kernel, grid_reso)


# 5. Heat map of islets that have at least one ad-cell NS component in a mantle
if pipeline_get_kde_islets_NSadcomps_inmantle:

    character_info = dict()

    for stage in range(4):
    
        character_info[stage] = []
    
        count = 0
    
        for islet in b_comps_in_mantle_info[stage]:

            # at least one adcomp in mantle
            if not len(islet['ad_comp_lens_inside_PHmantle']):
                continue

            comps_inside = np.array(islet['ad_comp_lens_inside_PHmantle'])

            NS_comps_inside = np.argwhere(comps_inside > 1).flatten()

            if len(NS_comps_inside) == 0:
                continue

            count += 1
    
            this_islet_characterization = islet['feature']
    
            character_info[stage].append(this_islet_characterization)

        print(f'Stage {stage} has {count} islets with at least NS ad comp in a mantle')
    
        NSad_in_b_isl_counts[stage] = count

    # 2b. Plot Kde maps and save kernels
    
    all_islet_kernels = dict()

    for stage in range(4):
    
        this_info = np.array(character_info[stage]).T
    
        # transform
        this_info[1,:] = transform(this_info[1,:])
    
        save_fname_prefix = 'stage'+str(stage)+'_atleast_one_NSadcomp_in_PHmantle'
    
        all_islet_kernels[stage] = get_2d_kde_plot(this_info, t_rranges, grid_reso, 0.1, save_fname_prefix)

    
    
    save_kernel = 'islets_atleast_one_NSadcomp_in_PHmantle.p'
    
    pickle.dump([all_islet_kernels, t_rranges], open(save_kernel, 'wb'))

    print('Computing kl div of islets with at least one NS adcomp in mantle')
    
    get_kldiv(save_kernel, grid_reso)
    







