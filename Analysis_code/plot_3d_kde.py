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
from scipy.spatial import distance
import ndtest

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

        entr = (entropy(data1, data2) + entropy(data2, data1))/2

        #print(f'kl div between stage {s1_idx} and stage {s2_idx} is {entr}')

        mmean = (data1 + data2)/2

        test_measure = math.sqrt((entropy(data1, mmean) + entropy(data2, mmean))/2)


        jensen_shannon = distance.jensenshannon(data1, data2)

        #print(f'Jensen-Shannon div between stage {s1_idx} and stage {s2_idx} is {jensen_shannon}')

        print(f'Stages {s1_idx} and {s2_idx}: (D(p||q) + D(q||p))/2 = {round(entr,3)}, sqrt is {round(math.sqrt(entr),3)}, jensen_shannon is {round(jensen_shannon,3)}, test measure is {test_measure}')






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


b_comps_in_mantle_info = pickle.load(open('B_COMPS_SIZES_IN_MANTLE_DICT.p', 'rb'))


###################
### pipeline #######
###################

pipeline_get_global_rranges = 1

pipeline_get_kde_all_islets = 1

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

# 2. Heat map of islet characterization
if pipeline_get_kde_all_islets:
    character_info = dict()

    for stage in range(4):
    
        character_info[stage] = []
    
    
        for islet in b_comps_in_mantle_info[stage]:
    
            this_islet_characterization = islet['feature']
    
            character_info[stage].append(this_islet_characterization)
    
    ## 2b. Plot Kde maps and save kernels
    #
    #all_islet_kernels = dict()

    #for stage in range(4):
    #
    #    this_info = np.array(character_info[stage]).T
    #
    #    # transform
    #    this_info[1,:] = transform(this_info[1,:])
    #
    #    save_fname_prefix = 'stage'+str(stage)+'_all_islets'
    #
    #    all_islet_kernels[stage] = get_2d_kde_plot(this_info, t_rranges, grid_reso, 0.1, save_fname_prefix)
    
    
    #save_kernel = 'all_islets_kernels.p'
    #
    #pickle.dump([all_islet_kernels, t_rranges], open(save_kernel, 'wb'))

    #print('Computing kl div of islets')
    #
    #get_kldiv(save_kernel, grid_reso)

    # Get 2D KS-test

    for s1_idx, s2_idx in it.combinations(list(range(4)), 2):

        #s1_idx = 0
        #s2_idx = 1

        data1 = np.array(character_info[s1_idx])
        data2 = np.array(character_info[s2_idx])


        data1[:,1] = transform(data1[:,1])
        data2[:,1] = transform(data2[:,1])

        P, D = ndtest.ks2d2s(data1[:,0], data1[:,1], data2[:,0], data2[:,1], nboot=1000, extra=True)

        print(f'stages {s1_idx} {s2_idx} p-val of 2D KS test is {P}\n')

        #p, en, enboot = ndtest.estat2d(data1[:,0], data1[:,1], data2[:,0], data2[:,1])

        #print(f'stages {s1_idx} {s2_idx} p-val of 2D energy dist is {p}')
        ##print(P, D)
    

# 3. Heat map of islets that have no NS b-cell components
if pipeline_get_kde_islets_with_NO_NS_bcomps:

    info = dict()

    for stage in range(4):
    
        info[stage] = []
    
        for islet in b_comps_in_mantle_info[stage]:

            if len(islet['b_comp_lens']):

                this_islet_info = np.array(islet['b_comp_lens'])
                idxs = np.argwhere(this_islet_info > 1).flatten()

                if len(idxs):
                    continue

                this_islet_characterization = islet['feature']
    
                info[stage].append(this_islet_characterization)
    
    all_islet_kernels = dict()
    
    for stage in range(4):

        print(f'There are {len(info[stage])} islets in stage {stage} without NS b components')

# 4. Heat map of islets that have no NS ad-cell components
if pipeline_get_kde_islets_with_NO_NS_adcomps:

    info = dict()

    for stage in range(4):
    
        info[stage] = []
    
        for islet in b_comps_in_mantle_info[stage]:

            if len(islet['ad_comp_lens']):

                this_islet_info = np.array(islet['ad_comp_lens'])
                idxs = np.argwhere(this_islet_info > 1).flatten()

                if len(idxs):
                    continue

                this_islet_characterization = islet['feature']
    
                info[stage].append(this_islet_characterization)
    
    all_islet_kernels = dict()
    
    for stage in range(4):

        print(f'There are {len(info[stage])} islets in stage {stage} without NS ad components')

# 5. Heat map of islets that have at least one b-cell component in a mantle
if pipeline_get_kde_islets_bcomps_inmantle:

    character_info = dict()

    for stage in range(4):
    
        character_info[stage] = []
    
    
        for islet in b_comps_in_mantle_info[stage]:

            # at least one bcomp in mantle
            if not len(islet['b_comp_lens_inside_mantle']):
                continue
    
            this_islet_characterization = islet['feature']
    
            character_info[stage].append(this_islet_characterization)
    
    # 2b. Plot Kde maps and save kernels
    
    all_islet_kernels = dict()

    for stage in range(4):
    
        this_info = np.array(character_info[stage]).T
    
        # transform
        this_info[1,:] = transform(this_info[1,:])
    
        save_fname_prefix = 'stage'+str(stage)+'_atleast_one_bcomp_inmantle'
    
        all_islet_kernels[stage] = get_2d_kde_plot(this_info, t_rranges, grid_reso, 0.1, save_fname_prefix)

    
    
    save_kernel = 'islets_atleast_one_bcomp_inmantle.p'
    
    pickle.dump([all_islet_kernels, t_rranges], open(save_kernel, 'wb'))

    print('Computing kl div of islets with at least one bcomp in mantle')
    
    get_kldiv(save_kernel, grid_reso)


# 5. Heat map of islets that have at least one b-cell NS component in a mantle
if pipeline_get_kde_islets_NSbcomps_inmantle:

    character_info = dict()

    for stage in range(4):
    
        character_info[stage] = []
    
        count = 0
    
        for islet in b_comps_in_mantle_info[stage]:

            # at least one bcomp in mantle
            if not len(islet['b_comp_lens_inside_mantle']):
                continue

            comps_inside = np.array(islet['b_comp_lens_inside_mantle'])

            NS_comps_inside = np.argwhere(comps_inside > 1).flatten()

            if len(NS_comps_inside) == 0:
                continue

            count += 1
    
            this_islet_characterization = islet['feature']
    
            character_info[stage].append(this_islet_characterization)
    

        NSb_in_ad_isl_counts[stage] = count

    # 2b. Plot Kde maps and save kernels
    
    all_islet_kernels = dict()

    for stage in range(4):
    
        this_info = np.array(character_info[stage]).T
    
        # transform
        this_info[1,:] = transform(this_info[1,:])
    
        save_fname_prefix = 'stage'+str(stage)+'_atleast_one_NSbcomp_inmantle'
    
        all_islet_kernels[stage] = get_2d_kde_plot(this_info, t_rranges, grid_reso, 0.1, save_fname_prefix)

    
    
    save_kernel = 'islets_atleast_one_NSbcomp_inmantle.p'
    
    pickle.dump([all_islet_kernels, t_rranges], open(save_kernel, 'wb'))

    print('Computing kl div of islets with at least one NS bcomp in mantle')
    
    get_kldiv(save_kernel, grid_reso)
    
# 5. Heat map of islets that have at least one ad-cell NS component in a mantle
if pipeline_get_kde_islets_NSadcomps_inmantle:

    character_info = dict()

    for stage in range(4):
    
        character_info[stage] = []
    
        count = 0
    
        for islet in b_comps_in_mantle_info[stage]:

            # at least one adcomp in mantle
            if not len(islet['ad_comp_lens_inside_mantle']):
                continue

            comps_inside = np.array(islet['ad_comp_lens_inside_mantle'])

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
    
        save_fname_prefix = 'stage'+str(stage)+'_atleast_one_NSadcomp_inmantle'
    
        all_islet_kernels[stage] = get_2d_kde_plot(this_info, t_rranges, grid_reso, 0.1, save_fname_prefix)

    
    
    save_kernel = 'islets_atleast_one_NSadcomp_inmantle.p'
    
    pickle.dump([all_islet_kernels, t_rranges], open(save_kernel, 'wb'))

    print('Computing kl div of islets with at least one NS adcomp in mantle')
    
    get_kldiv(save_kernel, grid_reso)
    

xx = []
yy1 = []
yy2 = []

for stage in range(4):

    xx.append(stage)

    yy1.append(NSb_in_ad_isl_counts[stage]/tot_isl_counts[stage]*100)
    yy2.append(NSad_in_b_isl_counts[stage]/tot_isl_counts[stage]*100)

plt.plot(xx, yy1, color='blue', lw=2, marker='o', label=r'islets with at least one NS $\beta$ in mantle')
plt.plot(xx, yy2, color='red', lw=2, marker='o', label=r'islets with at least one NS $\alpha\delta$ in mantle')

plt.legend()

#plt.xlabel('Developmental stages', fontsize=16)
plt.ylabel('Percentage of total islets', fontsize=16)
plt.xticks(ticks=[0,1,2,3], labels=['Stage 0', 'Stage 1', 'Stage 2', 'Stage 3'], fontsize=16)
plt.savefig('compare_NS_b_ad_in_mantle.pdf', dpi=600)



exit()





global_stat_ranges = dict()


stat_kinds = ['count_comps'\
            , 'max_comp'\
            , 'mean_comp'\
            , 'count_NS'\
            , 'max_NS'
            , 'mean_NS'
            ]

if pipeline_compare_ns_comps_inside_and_out:

    colors = ['blue'\
            , 'green'\
            , 'orange'\
            , 'red'\
            ]

    global_stat_x_min = 0
    global_stat_y_min = 0
    global_stat_x_max = 0
    global_stat_y_max = 0

    info_dict = dict()

    for stage in range(4):

        info_dict[stage] = []

        count_ns_inside = 0

        count_any_inside = 0


        for islet in b_comps_in_mantle_info[stage]:

            comps_inside = islet['b_comp_lens_inside_mantle']
            comps_not_inside = islet['b_comp_lens_not_inside_mantle']

            comps_inside = np.array(comps_inside)
            comps_not_inside = np.array(comps_not_inside)

            if len(comps_inside):
                count_any_inside += 1


            ns_inside = comps_inside[comps_inside > 1]
            ns_not_inside = comps_not_inside[comps_not_inside > 1]

            if len(ns_inside) == 0:
                continue

            count_ns_inside += 1

            ## for count
            #xx = len(ns_not_inside)
            #yy = len(ns_inside)

            # for max
            try:
                xx = np.amax(ns_not_inside)
            except:
                xx = 0

            try:
                yy = np.amax(ns_inside)
            except:
                yy = 0

            info_dict[stage].append([xx, yy])

        print(f'{count_any_inside} islets in stage {stage} have at least one b-component inside mantle')
        print(f'{count_ns_inside} islets in stage {stage} have at least one NS b-component inside mantle')
        info_dict[stage] = np.array(info_dict[stage])

        this_mat = info_dict[stage]

        global_stat_x_max = max(global_stat_x_max, np.amax(this_mat[:,0]))
        global_stat_y_max = max(global_stat_y_max, np.amax(this_mat[:,1]))

    xmin = global_stat_x_min
    xmax = global_stat_x_max
    ymin = global_stat_y_min
    ymax = global_stat_y_max

    rranges = np.array([[xmin, xmax], [ymin, ymax]])
    rranges = transform(rranges)
    xmin, xmax = rranges[0]
    ymin, ymax = rranges[1]

    for stage in range(4):

        data = info_dict[stage]

        # transform
        data = transform(data)

        kernel = stats.gaussian_kde(data.T)

        X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
        positions = np.vstack([X.ravel(), Y.ravel()])

        kernel_estimates = kernel(positions)

        Z = np.reshape(kernel_estimates.T, X.shape)

        plt.imshow(np.rot90(Z), cmap=plt.cm.gnuplot,
              extent=[xmin, xmax, ymin, ymax]\
                      , aspect=1)

        plt.show()
        plt.cla()
        plt.clf()

xx = []
yy = []

for stage in range(4):

    xx.append(stage)

    yy1.append(tot_isl_counts[stage]/NSb_in_ad_isl_counts[stage]*100)
    yy2.append(tot_isl_counts[stage]/NSad_in_b_isl_counts[stage]*100)

plt.plot(xx, yy1, color='blue', lw=4, marker='o')
plt.plot(xx, yy2, color='red', lw=4, marker='o')

plt.xlabel('Developmental stages')
plt.ylabel('Percentage of total islets')

plt.show()


exit()




global_x_min = math.inf
global_y_min = math.inf

global_x_max = 0
global_y_max = 0

# Get global min max range for x (b cell fraction) and y (total cells)
for stage in range(4):
    islets_info = b_comps_in_mantle_info[stage]

    for islet in islets_info:

        this_x = islet['feature'][0]
        this_y = islet['feature'][1]

        global_x_min = min(global_x_min, this_x)
        global_x_max = max(global_x_max, this_x)

        global_y_min = min(global_y_min, this_y)
        global_y_max = max(global_y_max, this_y)


rrange = [[global_x_min, global_x_max], [global_y_min, global_y_max]]



xmin = global_x_min
xmax = global_x_max
ymin = global_y_min
ymax = global_y_max

# transform y ranges
transform_y = transform(np.array([ymin, ymax]))
ymin = transform_y[0]
ymax = transform_y[1]



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




# One kde plot
def kde_3d(all_data, stage, stat_kind, global_stat_ranges):

    info = all_data[stage]

    info_array = []

    for idx, islet in enumerate(info):

        islet_b_fraction = islet['feature'][0]
        islet_total_cells = islet['feature'][1]
        data = islet['comp_lens']

        if data == []:
            islet_stat = 0
            # If you want to ignore islets without mantles
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

        info_array.append([islet_b_fraction, islet_total_cells, islet_stat])


    info_array = np.array(info_array)

    info_array[:,1] = transform(info_array[:,1])

    if stat_kind == 'mean':
        # there are nan in mean if there were 0 (S, NS) components in mantles
        info_array[:,2] = np.nan_to_num(info_array[:,2])
        info_array[:,2] = transform(info_array[:,2])
        z_label = r'mean size of comp inside a loop (ln(1+val))'
        file_label = 'mean'
    elif stat_kind == 'max':
        info_array[:,2] = transform(info_array[:,2])
        z_label = r'max size of comp inside a loop (ln(1+val))'
        file_label = 'max'
    elif stat_kind == 'count_NS':
        z_label = '\#NS comps inside loops'
        file_label = 'numNS'
    elif stat_kind == 'count_S':
        z_label = '\#S comps inside loops'
        file_label = 'numS'



    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(info_array[:,0]\
            , info_array[:,1]\
            , info_array[:,2]\
            , alpha=0.5)
    
    ax.set_xlabel(r'$\beta$-fraction')
    ax.set_ylabel(r'\#total cells (ln(1+val))')
    ax.set_zlabel(z_label)
    
    #plt.show()
    plt.cla()
    plt.clf()
    plt.close()

    #stat_min = np.amin(info_array[:,2])
    stat_ranges = transform(np.array([global_stat_ranges[stat_kind]['min']\
                                    , global_stat_ranges[stat_kind]['max']]))
    stat_min = stat_ranges[0]
    # MAKING STAT_MIN ALWAYS 0
    stat_min = 0

    stat_max = stat_ranges[1]

    xx = info_array[:,0]
    yy = info_array[:,1]
    zz = info_array[:,2]

    # kde plots in 2 dimensions
    # x vs y

    values = np.vstack([xx, yy])
    kernel_xy = stats.gaussian_kde(values)


    X, Y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]
    positions = np.vstack([X.ravel(), Y.ravel()])

    Z = np.reshape(kernel_xy(positions).T, X.shape)

    plt.imshow(np.rot90(Z), cmap=plt.cm.gnuplot,
              extent=[xmin, xmax, ymin, ymax]\
                      , aspect=0.1)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.xlabel(r'$\beta$-fraction')
    plt.ylabel(r'total cells')

    
    plt.savefig('stage'+str(stage)+'_'+file_label+'_kde_xy.pdf')
    #plt.show()
    plt.cla()
    plt.clf()

    # kde plots in 2 dimensions
    # x vs z

    values = np.vstack([xx, zz])
    kernel_xz = stats.gaussian_kde(values)


    X, Y = np.mgrid[xmin:xmax:50j, stat_min:stat_max:50j]
    positions = np.vstack([X.ravel(), Y.ravel()])

    Z = np.reshape(kernel_xz(positions).T, X.shape)

    plt.imshow(np.rot90(Z), cmap=plt.cm.gnuplot,
              extent=[xmin, xmax, stat_min, stat_max]\
                      , aspect=0.1)
    plt.xlim([xmin, xmax])
    plt.ylim([stat_min, stat_max])
    plt.xlabel(r'$\beta$-fraction')
    plt.ylabel(z_label)

    
    plt.savefig('stage'+str(stage)+'_'+file_label+'_kde_xz.pdf')
    #plt.show()
    plt.cla()
    plt.clf()


    # kde plots in 2 dimensions
    # y vs z

    values = np.vstack([yy, zz])
    kernel_yz = stats.gaussian_kde(values)


    X, Y = np.mgrid[ymin:ymax:50j, stat_min:stat_max:50j]
    positions = np.vstack([X.ravel(), Y.ravel()])

    Z = np.reshape(kernel_yz(positions).T, X.shape)

    plt.imshow(np.rot90(Z), cmap=plt.cm.gnuplot,
              extent=[ymin, ymax, stat_min, stat_max]\
                      , aspect=1)
    plt.xlim([ymin, ymax])
    plt.ylim([stat_min, stat_max])
    plt.xlabel(r'total cells')
    plt.ylabel(z_label)

    
    plt.savefig('stage'+str(stage)+'_'+file_label+'_kde_yz.pdf')
    #plt.show()
    plt.cla()
    plt.clf()

    all_kernels = [kernel_xy, kernel_xz, kernel_yz]

    pickle.dump(all_kernels, open('stage'+str(stage)+'_'+file_label+'_kernels.p', 'wb'))



    ## Skip 3d KDE

    #X, Y, Z = np.mgrid[xmin:xmax:50j, transform_y[0]:transform_y[1]:50j, stat_min:stat_max:10j]
    #positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

    ## islet kde
    #m1 = info_array[:,0]
    #m2 = info_array[:,1]
    #m3 = info_array[:,2]
    #
    #values = np.vstack([m1, m2, m3])
    #kernel = stats.gaussian_kde(values, bw_method=1)

    #kernel_estimates = kernel(positions)

    ##S = np.reshape(kernel(positions).T, X.shape)

    #fig = go.Figure(data=go.Volume(
    #x=X.flatten(),
    #y=Y.flatten(),
    #z=Z.flatten(),
    #value=kernel_estimates,
    #opacity=0.5,
    #))
  
    #fig.show()

    ##print(S)
    ##exit()



    return info_array


# Stat kinds
# mean component size inside mantle
# max component size inside mantle
# num of NS components inside mantle
# num of singular components inside mantle

stat_kinds = ['mean'\
            , 'max'\
            , 'count_NS'\
            , 'count_S']

# Get stat ranges
# global stat ranges
global_stat_ranges = dict()

for stat_kind in stat_kinds:

    global_stat_ranges[stat_kind] = {'min':math.inf, 'max':0}

    get_stat_ranges(b_comps_in_mantle_info, stat_kind, global_stat_ranges)


for stage in range(4):

    for stat_kind in stat_kinds:

        print(stage, stat_kind)

        kde_3d(b_comps_in_mantle_info, stage, stat_kind, global_stat_ranges)










