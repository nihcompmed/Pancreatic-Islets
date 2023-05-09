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




# Load kernels for geom


# admantle around ns beta

# geom
geom_kernel_file = 'islets_atleast_one_NSbcomp_inmantle.p'
geom_kernels = pickle.load(open(geom_kernel_file, 'rb'))

geom_all_kernels = geom_kernels[0]
g_rranges = geom_kernels[1]



# ph
ph_kernel_file = 'islets_atleast_one_NSbcomp_in_PHmantle.p'
ph_kernels = pickle.load(open(ph_kernel_file, 'rb'))

ph_all_kernels = ph_kernels[0]
ph_rranges = ph_kernels[1]

assert np.all(g_rranges == ph_rranges)

grid_reso = 100j

rranges = g_rranges

xmin, xmax = rranges[0]
ymin, ymax = rranges[1]

X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
positions = np.vstack([X.ravel(), Y.ravel()])

for stage in range(4):
    
    k1 = geom_all_kernels[stage][0]
    k2 = ph_all_kernels[stage][0]


    data1 = k1(positions)
    data2 = k2(positions)

    entr = entropy(data1, data2)


    print(f'kl div between geom and ph for ad mantles around NS b comp for stage {stage} is {entr}')



# bmantle around ns alpha-delta

# geom
geom_kernel_file = 'islets_atleast_one_NSadcomp_inmantle.p'
geom_kernels = pickle.load(open(geom_kernel_file, 'rb'))

geom_all_kernels = geom_kernels[0]
g_rranges = geom_kernels[1]



# ph
ph_kernel_file = 'islets_atleast_one_NSadcomp_in_PHmantle.p'
ph_kernels = pickle.load(open(ph_kernel_file, 'rb'))

ph_all_kernels = ph_kernels[0]
ph_rranges = ph_kernels[1]

assert np.all(g_rranges == ph_rranges)

grid_reso = 100j

rranges = g_rranges

xmin, xmax = rranges[0]
ymin, ymax = rranges[1]

X, Y = np.mgrid[xmin:xmax:grid_reso, ymin:ymax:grid_reso]
positions = np.vstack([X.ravel(), Y.ravel()])

for stage in range(4):
    
    k1 = geom_all_kernels[stage][0]
    k2 = ph_all_kernels[stage][0]


    data1 = k1(positions)
    data2 = k2(positions)

    entr = entropy(data1, data2)


    print(f'kl div between geom and ph for b mantles around NS ad comp for stage {stage} is {entr}')



















