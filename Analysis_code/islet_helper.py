import numpy as np
import math
import itertools as it
import pickle
import matplotlib as mpl
import os
import seaborn as sns
import networkx as nx
import pandas as pd
import scipy
from statistics import variance

from matplotlib import rc
rc('text',usetex=True)
rc('text.latex', preamble=r'\usepackage{xcolor}')

import matplotlib.pyplot as plt


# For Debbie's code, loop is pair of coordinates of edges of loop
# For PH, loop is set of edges of loop
# other_cell is [x,y] of the cell to check for containment

def PH_loop_has_cell(cell_locs, edge_arr, loop, other_cell):


    loop_coords = []

    for edge in loop:

        c1 = edge[0]
        c2 = edge[1]

        x1 = cell_locs[c1][0]
        y1 = cell_locs[c1][1]

        x2 = cell_locs[c2][0]
        y2 = cell_locs[c2][1]

        loop_coords.append([[x1, y1], [x2, y2]])

    
# loop_coords is list [ [[x1, y1], [x2,y2]], [...], ...     ]
def geom_loop_has_cell(loop_coords, check_cell):

    x, y = check_cell

    intersect = 0

    for edge in loop_coords:

        p1 = edge[0]
        p2 = edge[1]

        x1, y1 = p1
        x2, y2 = p2

        # If horizontal line, check if point is on the edge
        if y1 == y2:
            xmin = min(x1, x2)
            xmax = max(x1, x2)
            if (x > xmin and x < xmax) or (x == xmin) or (x == xmax):
                return 1

        # If horizontal line, check if point is on the edge
        if y1 == y2:
            xmin = min(x1, x2)
            xmax = max(x1, x2)
            if (x > xmin and x < xmax) or (x == xmin) or (x == xmax):
                return 1
            else:
                if x < xmin:
                    intersect += 1
        else:
            # Otherwise check intersection (horizontal ray algorithm)
            t = (y-y1)/(y2-y1)
            if (t > 0 and t < 1):
                x_intersect = (1-t)*x1 + t*x2
                if x_intersect >= x:
                    intersect += 1
            elif (t == 0):
                intersect += 0.5
            elif (t == 1):
                intersect += 0.5

    if (intersect%2==0):
        return 0
    else:
        return 1


def max_dist_loop_cell(loop_coords, cell):

    max_dist = 0

    x, y = cell

    for edge in loop_coords:

        x1, y1 = edge[0]
        x2, y2 = edge[1]

        ddist = (x1-x)**2 + (y1-y)**2

        if (ddist > max_dist):
            max_dist = ddist

        ddist = (x2-x)**2 + (y2-y)**2

        if (ddist > max_dist):
            max_dist = ddist


    return max_dist



def max_min_dist_loop_loop_one_way(loop_coords1, loop_coords2):

    max_dist = 0
    for pt1 in loop_coords1:

        x1, y1 = pt1
        min_dist = math.inf

        for pt2 in loop_coords2:
            x2, y2 = pt2

            ddist = (x1-x2)**2 + (y1-y2)**2
            if ddist < min_dist:
                min_dist = ddist

        if min_dist > max_dist:
            max_dist = min_dist



    return math.sqrt(max_dist)

# loop_coords are list []
def max_min_dist_loop_loop(loop_coords1, loop_coords2):

    dist1 = max_min_dist_loop_loop_one_way(loop_coords1, loop_coords2)
    dist2 = max_min_dist_loop_loop_one_way(loop_coords2, loop_coords1)

    max_dist = max(dist1, dist2)


    return max_dist

def get_geom_loops(geom_loop_file):

   loops = [] 

   vv = open(geom_loop_file, 'r')

   for line in vv:


       line = line.strip('\n')
       line = line.split(':')

       cells_in_loop = line[0]

       cells_in_loop = cells_in_loop.split(',')
       cells_in_loop = cells_in_loop[:-1]
       cells_in_loop = [int(x) for x in cells_in_loop]

       geom_loop = line[1]
       geom_loop = geom_loop.split(',')
       geom_loop = geom_loop[:-1]
       geom_loop = [float(x) for x in geom_loop]
       nn = 0
       geom_loop_v2 = []
       while (nn < len(geom_loop) - 4):
           x1 = geom_loop[nn]
           y1 = geom_loop[nn+1]
           e1 = [x1, y1]

           x2 = geom_loop[nn+2]
           y2 = geom_loop[nn+3]
           e2 = [x2, y2]

           edge = [[x1, y1], [x2, y2]]
           geom_loop_v2.append(edge)

           nn += 2

       this_loop = [cells_in_loop, geom_loop_v2]

       loops.append(this_loop)

   vv.close()

   return loops

def get_PH_loops(PH_loops_file, inv_map):

    loops = []

    PH = pickle.load(open(PH_loops_file, 'rb'))

    PH_loops = PH[0]
    PH_info = PH[1]


    for these_loops in PH_info:

        if (len(these_loops)) == 0:
            continue

        this_cells_in_loop = PH_info[these_loops]
        this_cells_in_loop = [inv_map[cell] for cell in this_cells_in_loop]

        for this_loop_idx in these_loops:

            this_loop = PH_loops[this_loop_idx]

            loops.append([this_cells_in_loop, this_loop])

    return loops


def plot_bfrac_total(all_info, fname, colors, stages, plt_rows, plt_cols):

    idx = 0
    
    fig, axs = plt.subplots(plt_rows, plt_cols, sharex=True, sharey=True)
    
    row = 0
    col = 0
    
    for stage in all_info:
    
        sub_total = []
        sub_b_frac = []
        sub_density = []
    
        for subject in all_info[stage]:
    
    
            for islet in all_info[stage][subject]:
    
                this_info = all_info[stage][subject][islet]
                n_ad = this_info['n_ad']
                n_b = this_info['n_b']
                area = this_info['area']
    
                total = n_ad + n_b
    
                if n_ad < 6 or n_b < 6:
                    continue
    
                sub_total.append(total)
                sub_b_frac.append(n_b/total)
                sub_density.append(total/area)
    
    
        #axs[row, col].scatter(sub_total, sub_b_frac, color=colors[idx], label=stages[idx], alpha=0.4)
        #axs[row, col].scatter(sub_total, sub_density, color=colors[idx], label=stages[idx], alpha=0.4)
    
        sub_total   = np.array(sub_total)
        sub_b_frac  = np.array(sub_b_frac)
    
        sub_density = np.array(sub_density)
    
        sub_total   = np.log(sub_total)
        sub_density = np.log(sub_density)
    
        if plt_rows > 1:
            axs[row, col].hist2d(sub_total\
                               # , sub_density\
                                , sub_b_frac\
                                , bins=50\
                                , norm=mpl.colors.LogNorm()\
                                , cmap=mpl.cm.gnuplot)
    
            axs[row, col].axhline(y=0.5, ls='--', color='black')
            axs[row, col].set_title(label=stages[idx])
        else:
            axs[col].hist2d(sub_total\
                          # , sub_density\
                           , sub_b_frac\
                           , bins=50\
                           , norm=mpl.colors.LogNorm()\
                           , cmap=mpl.cm.gnuplot)
    
            axs[col].axhline(y=0.5, ls='--', color='black')
            axs[col].set_title(label=stages[idx])
    
        col += 1
        if col == 2:
            row += 1
            col = 0
    
    
    
        idx += 1
    
    if plt_rows > 1:
        axs[1, 0].set_xlabel('total cells (log)')
        axs[1, 0].set_ylabel('beta cell fraction')
    else:
        axs[0].set_xlabel('total cells (log)')
        axs[0].set_ylabel('beta cell fraction')
    
    plt.savefig(fname)
    
    plt.cla()
    plt.clf()
    

def plot_bhavecyc_frac(all_info, fname, colors, stages):

    
    idx = 0
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
    
    row = 0
    col = 0
    
    for stage in all_info:
    
        sub_total = []
        sub_b_frac = []
        sub_b = []
        sub_ad = []
        sub_density = []

        b_in_mantle = []
        ad_in_mantle = []

        sub_b_in_mantle_frac = []
    
        for subject in all_info[stage]:
    
    
            for islet in all_info[stage][subject]:
    
                this_info = all_info[stage][subject][islet]
                n_ad = this_info['n_ad']
                n_b = this_info['n_b']
                area = this_info['area']

                gor = this_info['gor']


                if n_ad < 6 or n_b < 6:
                    continue

                n_b_in_mantle = 0

                total = n_ad + n_b

                if not gor:
                    continue
    
                sub_total.append(total)
                sub_b_frac.append(n_b/total)
                sub_b.append(n_b)
                sub_ad.append(n_ad)
                sub_density.append(total/area)

                if this_info['b_in_ad_geom'] == 0:
                    b_in_mantle.append(0)
                    sub_b_in_mantle_frac.append(0)
                else:
                    for loop in this_info['b_in_ad_geom']:
                        bb = loop[0]
                        n_b_in_mantle += len(bb)
        
                    sub_b_in_mantle_frac.append(n_b_in_mantle/n_b)
                    b_in_mantle.append(n_b_in_mantle)
    
    
        sub_b_in_mantle_frac  = np.array(sub_b_in_mantle_frac)

        b_in_mantle  = np.array(b_in_mantle)
        b_in_mantle  = np.log(1+b_in_mantle)

        sub_b  = np.array(sub_b)
        sub_b  = np.log(sub_b)

        sub_ad  = np.array(sub_ad)
        sub_ad  = np.log(sub_ad)

    
        axs[row, col].hist2d(sub_b_frac\
                           # , sub_density\
                           # , sub_b_in_mantle_frac\
                            , sub_b_in_mantle_frac\
                            , bins=50\
                            , norm=mpl.colors.LogNorm()\
                            , cmap=mpl.cm.gnuplot)
    
        axs[row, col].set_title(label=stages[idx])
    
        col += 1
        if col == 2:
            row += 1
            col = 0
    
    
    
        idx += 1
    
    #axs[0, 0].set_xlabel('total cells')
    #axs[0, 0].set_ylabel('beta cell fraction')
    #
    #plt.savefig(fname)
    #
    #plt.cla()
    #plt.clf()


    plt.show()




def plot_adhavecyc_frac(all_info, fname, colors, stages):

    
    idx = 0
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
    
    row = 0
    col = 0
    
    for stage in all_info:
    
        sub_total = []
        sub_b_frac = []
        sub_b = []
        sub_ad = []
        sub_density = []

        b_in_mantle = []
        ad_in_mantle = []

        sub_ad_in_mantle_frac = []
    
        for subject in all_info[stage]:
    
    
            for islet in all_info[stage][subject]:
    
                this_info = all_info[stage][subject][islet]
                n_ad = this_info['n_ad']
                n_b = this_info['n_b']
                area = this_info['area']

                gor = this_info['gor']


                if n_ad < 6 or n_b < 6:
                    continue

                n_ad_in_mantle = 0

                total = n_ad + n_b

                if not gor:
                    continue
    
                sub_total.append(total)
                sub_b_frac.append(n_b/total)
                sub_b.append(n_b)
                sub_ad.append(n_ad)
                sub_density.append(total/area)

                if this_info['ad_in_b_geom'] == 0:
                    ad_in_mantle.append(0)
                else:
                    for loop in this_info['ad_in_b_geom']:
                        bb = loop[0]
                        n_ad_in_mantle += len(bb)
        
                sub_ad_in_mantle_frac.append(n_ad_in_mantle/n_ad)
    
    
        sub_ad_in_mantle_frac  = np.array(sub_ad_in_mantle_frac)

        sub_b  = np.array(sub_b)
        sub_b  = np.log(sub_b)

        sub_ad  = np.array(sub_ad)
    
        sub_ad  = np.log(sub_ad)
    
        axs[row, col].hist2d(sub_ad\
                           # , sub_density\
                           # , sub_b_in_mantle_frac\
                            , sub_b\
                            , bins=50\
                            , norm=mpl.colors.LogNorm()\
                            , cmap=mpl.cm.gnuplot)
    
        axs[row, col].set_title(label=stages[idx])
    
        col += 1
        if col == 2:
            row += 1
            col = 0
    
    
    
        idx += 1
    
    #axs[0, 0].set_xlabel('total cells')
    #axs[0, 0].set_ylabel('beta cell fraction')
    #
    #plt.savefig(fname)
    #
    #plt.cla()
    #plt.clf()


    plt.show()




def plot_bhavecyc_islets(all_info, fname, colors, stages):

    
    idx = 0
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
    
    row = 0
    col = 0
    
    for stage in all_info:
    
        sub_total = []
        sub_b_frac = []
        sub_b = []
        sub_ad = []
        sub_density = []

        b_in_mantle = []
        ad_in_mantle = []

        sub_b_in_mantle_frac = []

        sub_ids = []

        sub_id = 0
    
        for subject in all_info[stage]:
    
            sub_id += 1
    
            for islet in all_info[stage][subject]:
    
                this_info = all_info[stage][subject][islet]
                n_ad = this_info['n_ad']
                n_b = this_info['n_b']
                area = this_info['area']

                gor = this_info['gor']


                if n_ad < 6 or n_b < 6:
                    continue

                n_b_in_mantle = 0

                total = n_ad + n_b

                if not gor:
                    continue

                sub_ids.append(sub_id)
    
                sub_total.append(total)
                sub_b_frac.append(n_b/total)
                sub_b.append(n_b)
                sub_ad.append(n_ad)
                sub_density.append(total/area)

                if this_info['b_in_ad_geom'] == 0:
                    b_in_mantle.append(0)
                    sub_b_in_mantle_frac.append(0)
                else:
                    for loop in this_info['b_in_ad_geom']:
                        bb = loop[0]
                        n_b_in_mantle += len(bb)
        
                    sub_b_in_mantle_frac.append(n_b_in_mantle/n_b)
                    b_in_mantle.append(n_b_in_mantle)
    
    
        sub_b_in_mantle_frac  = np.array(sub_b_in_mantle_frac)

        b_in_mantle  = np.array(b_in_mantle)
        b_in_mantle  = np.log(1+b_in_mantle)

        sub_b  = np.array(sub_b)
        sub_b  = np.log(sub_b)

        sub_ad  = np.array(sub_ad)
        sub_ad  = np.log(sub_ad)

    
        #axs[row, col].hist2d(sub_b_frac\
        #                   # , sub_density\
        #                   # , sub_b_in_mantle_frac\
        #                    , sub_b_in_mantle_frac\
        #                    , bins=50\
        #                    , norm=mpl.colors.LogNorm()\
        #                    , cmap=mpl.cm.gnuplot)

        axs[row, col].scatter(sub_ids, b_in_mantle)
    
        axs[row, col].set_title(label=stages[idx])
    
        col += 1
        if col == 2:
            row += 1
            col = 0
    
    
    
        idx += 1
    
    #axs[0, 0].set_xlabel('total cells')
    #axs[0, 0].set_ylabel('beta cell fraction')
    #
    #plt.savefig(fname)
    #
    #plt.cla()
    #plt.clf()


    plt.show()




def plot_bad_adb(all_info, fname, colors, stages, plt_rows, plt_cols):
    
    idx = 0
    
    fig, axs = plt.subplots(plt_rows, plt_cols, sharex=True, sharey=True)
    
    row = 0
    col = 0
    
    for stage in all_info:
    
        sub_total = []
        sub_b_frac = []


        sub_density = []

        sub_b = []
        sub_ad = []

        b_in_mantle = []
        ad_in_mantle = []

        sub_b_in_mantle_frac = []

        sub_ids = []

        sub_id = 0
    
        for subject in all_info[stage]:
    
            sub_id += 1
    
            for islet in all_info[stage][subject]:
    
                this_info = all_info[stage][subject][islet]
                n_ad = this_info['n_ad']
                n_b = this_info['n_b']
                area = this_info['area']

                gor = this_info['gor']


                if n_ad < 6 or n_b < 6:
                    continue



                if not gor:
                    continue

                total = n_ad + n_b

                n_ad_in_mantle = 0
                n_b_in_mantle = 0

                sub_ids.append(sub_id)
    
                sub_total.append(total)
                sub_b_frac.append(n_b/total)

                sub_b.append(n_b)
                sub_ad.append(n_ad)

                sub_density.append(total/area)

                if this_info['b_in_ad_geom'] == 0:
                    b_in_mantle.append(0)
                    sub_b_in_mantle_frac.append(0)
                else:
                    for loop in this_info['b_in_ad_geom']:
                        bb = loop[0]
                        n_b_in_mantle += len(bb)
        
                    sub_b_in_mantle_frac.append(n_b_in_mantle/n_b)
                    b_in_mantle.append(n_b_in_mantle)

                if this_info['ad_in_b_geom'] == 0:
                    ad_in_mantle.append(0)
                else:
                    for loop in this_info['ad_in_b_geom']:
                        bb = loop[0]
                        n_ad_in_mantle += len(bb)
                    ad_in_mantle.append(n_ad_in_mantle)
    
    
        ad_in_mantle  = np.array(ad_in_mantle)
        ad_in_mantle  = np.log(1+ad_in_mantle)

        b_in_mantle  = np.array(b_in_mantle)
        b_in_mantle  = np.log(1+b_in_mantle)

        sub_b  = np.array(sub_b)
        sub_b  = np.log(sub_b)

        sub_ad  = np.array(sub_ad)
        sub_ad  = np.log(sub_ad)

    
        #axs[row, col].hist2d(sub_b_frac\
        #                    , b_in_mantle_frac\
        #                    , bins=50\
        #                    , norm=mpl.colors.LogNorm()\
        #                    , cmap=mpl.cm.gnuplot)

        if plt_rows > 1:

            axs[row, col].scatter(sub_ad, b_in_mantle, alpha=0.2, color='blue', label='b cells')
            axs[row, col].scatter(sub_b, ad_in_mantle, alpha=0.2, color='red', label='ad cells')
            axs[row, col].set_title(label=stages[idx])

        else:
            axs[col].scatter(sub_ad, b_in_mantle, alpha=0.2, color='blue', label='b cells')
            axs[col].scatter(sub_b, ad_in_mantle, alpha=0.2, color='red', label='ad cells')
            axs[col].set_title(label=stages[idx])
    
        col += 1
        if col == 2:
            row += 1
            col = 0
    
    
    
        idx += 1

    plt.legend()
    
    #axs[0, 0].set_xlabel('total cells')
    #axs[0, 0].set_ylabel('beta cell fraction')
    #
    #plt.savefig(fname)
    #
    #plt.cla()
    #plt.clf()

    #rc('text',usetex=True)
    #rc('text.latex', preamble=r'\usepackage{xcolor}')

    #if plt_rows > 1:
    #    axs[1,0].set_xlabel(r'\textcolor{blue}{ad cells}'+r'\textcolor{red}{b cells} (log)')
    #    axs[1,0].set_ylabel(r'\textcolor{blue}{b cells}, \textcolor{red}{ad cells} in mantle (log)')
    #else:
    #    axs[0].set_xlabel(r'\textcolor{blue}{ad cells}, \textcolor{red}{b cells} (log)')
    #    axs[0].set_ylabel(r'\textcolor{blue}{b cells}, \textcolor{red}{ad cells} in mantle (log)')

    #plt.show()

    plt.savefig(fname)
    plt.cla()
    plt.clf()




def plot_bad_adb_stats(all_info, fname, colors, stages, plt_rows, plt_cols):
    
    idx = 0
    
    fig, axs = plt.subplots(plt_rows, plt_cols, sharex=True, sharey=True)
    
    row = 0
    col = 0
    
    for stage in all_info:
    
        sub_total = []
        sub_b_frac = []


        sub_density = []

        sub_b = []
        sub_ad = []
        
        xx_b_in_mantle = []
        xx_ad_in_mantle = []

        b_in_mantle = []
        ad_in_mantle = []

        sub_b_in_mantle_frac = []

        sub_ids = []

        sub_id = 0
    
        for subject in all_info[stage]:
    
            sub_id += 1
    
            for islet in all_info[stage][subject]:
    
                this_info = all_info[stage][subject][islet]
                n_ad = this_info['n_ad']
                n_b = this_info['n_b']
                area = this_info['area']

                gor = this_info['gor']


                if n_ad < 6 or n_b < 6:
                    continue



                if not gor:
                    continue

                total = n_ad + n_b

                n_ad_in_mantle = 0
                n_b_in_mantle = 0


                #sub_density.append(total/area)

                if this_info['b_in_ad_geom'] == 0:
                    b_in_mantle.append(0)
                    sub_b_in_mantle_frac.append(0)
                    sub_ids.append(sub_id)
    
                    sub_total.append(total)
                    sub_b_frac.append(n_b/total)

                    xx_b_in_mantle.append(n_ad)
                    sub_ad.append(n_ad)
                else:
                    for loop in this_info['b_in_ad_geom']:
                        bb = len(loop[0])
                        b_in_mantle.append(bb)
                        sub_ids.append(sub_id)
    
                        sub_total.append(total)
                        sub_b_frac.append(n_b/total)

                        xx_b_in_mantle.append(n_ad)


                        sub_b.append(n_b)
                        sub_ad.append(n_ad)

                if this_info['ad_in_b_geom'] == 0:
                    ad_in_mantle.append(0)
                    sub_ids.append(sub_id)
    
                    sub_total.append(total)
                    sub_b_frac.append(n_b/total)

                    xx_ad_in_mantle.append(n_b)

                    sub_b.append(n_b)
                    sub_ad.append(n_ad)
                else:
                    for loop in this_info['ad_in_b_geom']:
                        bb = len(loop[0])

                        #n_ad_in_mantle += len(bb)

                        ad_in_mantle.append(bb)
                        sub_ids.append(sub_id)

                        xx_ad_in_mantle.append(n_b)
    
                        sub_total.append(total)
                        sub_b_frac.append(n_b/total)

                        sub_b.append(n_b)
                        sub_ad.append(n_ad)
    
    
        ad_in_mantle  = np.array(ad_in_mantle)
        ad_in_mantle  = np.log(1+ad_in_mantle)

        b_in_mantle  = np.array(b_in_mantle)
        b_in_mantle  = np.log(1+b_in_mantle)

        sub_b  = np.array(sub_b)
        sub_b  = np.log(sub_b)

        sub_ad  = np.array(sub_ad)
        sub_ad  = np.log(sub_ad)

        xx_b_in_mantle = np.array(xx_b_in_mantle)
        xx_ad_in_mantle = np.array(xx_ad_in_mantle)

        xx_b_in_mantle = np.log(1+xx_b_in_mantle)
        xx_ad_in_mantle = np.log(1+xx_ad_in_mantle)
    
        #axs[row, col].hist2d(sub_b_frac\
        #                    , b_in_mantle_frac\
        #                    , bins=50\
        #                    , norm=mpl.colors.LogNorm()\
        #                    , cmap=mpl.cm.gnuplot)

        if plt_rows > 1:

            axs[row, col].scatter(xx_b_in_mantle, b_in_mantle, alpha=0.2, color='blue', label='b cells')
            axs[row, col].scatter(xx_ad_in_mantle, ad_in_mantle, alpha=0.2, color='red', label='ad cells')
            axs[row, col].set_title(label=stages[idx])

        else:
            axs[col].scatter(xx_b_in_mantle, b_in_mantle, alpha=0.2, color='blue', label='b cells')
            axs[col].scatter(xx_ad_in_mantle, ad_in_mantle, alpha=0.2, color='red', label='ad cells')
            axs[col].set_title(label=stages[idx])
    
        col += 1
        if col == 2:
            row += 1
            col = 0
    
    
    
        idx += 1

    plt.legend()
    
    #axs[0, 0].set_xlabel('total cells')
    #axs[0, 0].set_ylabel('beta cell fraction')
    #
    #plt.savefig(fname)
    #
    #plt.cla()
    #plt.clf()

    #rc('text',usetex=True)
    #rc('text.latex', preamble=r'\usepackage{xcolor}')

    #if plt_rows > 1:
    #    axs[1,0].set_xlabel(r'\textcolor{blue}{ad cells}'+r'\textcolor{red}{b cells} (log)')
    #    axs[1,0].set_ylabel(r'\textcolor{blue}{b cells}, \textcolor{red}{ad cells} in mantle (log)')
    #else:
    #    axs[0].set_xlabel(r'\textcolor{blue}{ad cells}, \textcolor{red}{b cells} (log)')
    #    axs[0].set_ylabel(r'\textcolor{blue}{b cells}, \textcolor{red}{ad cells} in mantle (log)')

    #plt.show()

    plt.savefig(fname)
    plt.cla()
    plt.clf()


def get_all_PD_info(info_file, save_name):

    #ff = open('all_islets_info.csv', 'r')
    ff = open(info_file, 'r')
    
    ff.readline()
    
    all_ad_files = []
    all_b_files = []
    
    
    all_H1_PD = dict()
    
    
    for line in ff:
    
        info = line.split(',')
    
        fprefix = info[3]
    
    
        PH_ad_file = 'PD_results/' + fprefix + '_ad.csv'

        if os.path.isfile(PH_ad_file):
            pts = np.loadtxt(PH_ad_file, delimiter=',')
            n_pts= len(pts)
    
            target = 'PD_results/' + fprefix + '_ad_H1_pers_data.txt'
    
            if os.path.isfile(target):
                this_PD = np.loadtxt(target, delimiter=',')
                if len(this_PD):
                    if this_PD.ndim == 1:
                        this_PD = np.reshape(this_PD, (1, 2))
                    max_pers = np.amax(this_PD[:,1] - this_PD[:,0])
                    all_H1_PD[fprefix+'_ad'] = {\
                                             'n_pts':n_pts\
                                            , 'PD':this_PD\
                                            , 'max_pers':max_pers\
                                            }
    
        PH_b_file = 'PD_results/' + fprefix + '_b.csv'
        if os.path.isfile(PH_b_file):
            pts = np.loadtxt(PH_b_file, delimiter=',')
            n_pts= len(pts)
    
            target = 'PD_results/' + fprefix + '_b_H1_pers_data.txt'
    
            if os.path.isfile(target):
                this_PD = np.loadtxt(target, delimiter=',')
                if len(this_PD):
                    if this_PD.ndim == 1:
                        this_PD = np.reshape(this_PD, (1, 2))
                    max_pers = np.amax(this_PD[:,1] - this_PD[:,0])
                    all_H1_PD[fprefix+'_b'] = {\
                                             'n_pts':n_pts\
                                            , 'PD':this_PD\
                                            , 'max_pers':max_pers\
                                            }
    
    
    ff.close()

    pickle.dump(all_H1_PD, open(save_name, 'wb'))
    
    
def plot_max_pers_distri_dev(PD_info, fname, typ):
    
    print('Inside plot')
    #plt.figure(figsize=(10,8))

    dev_stages = [\
            'Gestation/'\
            ,'1-35weeks/'\
            ,'12-24months/'
            ,'28MandAfter/'\
                    ]

    max_pers_dict = dict()

    stages = []
    max_pers = []

    G = nx.Graph()

    #progress = 0

    for key in PD_info:

        #print(progress)
        #progress += 1

        stage = key.split('_')[0][-1]

        G.add_node(stage)

        this_typ = key.split('_')[-1]

        if typ != 'all':
            if this_typ != typ:
                continue


        if stage not in max_pers_dict:
            max_pers_dict[stage] = []

        max_pers_dict[stage].append(PD_info[key]['max_pers'])

        stages.append(stage)
        max_pers.append(PD_info[key]['max_pers'])

    ptiles_dict = dict()

    for stage in range(4):
        vals = max_pers_dict[str(stage)]
        vals = np.array(vals)
        ptiles = np.percentile(vals, [5, 50, 95])
        ptiles_dict[stage] = ptiles

    d = {'Stages': stages, 'Max persistence': max_pers}
    df = pd.DataFrame(data=d)

    sns.set(font_scale=3)
    sns.set_style("whitegrid")

    sns.boxplot(data = df, x='Stages', y='Max persistence', palette=['lightblue', 'pink', 'lightgreen',
                                                                     'orange'])

    #plt.yscale('log', base=2)
    plt.tight_layout()

    #plt.title(title)
    plt.savefig(fname)
    #plt.show()



    plt.cla()
    plt.clf()


    print('Computing mann whitney')

    title = ''

    for s1, s2 in it.combinations(list(max_pers_dict.keys()), 2):

        res = get_mannwhitney(max_pers_dict[s1], max_pers_dict[s2], axis=0)
        print(s1, s2, res.pvalue)
        title += s1 + s2 + ': ' + f"{res.pvalue:.1E}" + ', '

        G.add_edge(s1, s2, weight=1-res.pvalue)

    dists = nx.get_edge_attributes(G, "weight")

    edge_labels = dict()
    edge_list = []
    edge_style = []
    edge_color = []
    edge_width = []

    for edge in dists:

        edge_list.append(edge)

        val = 1-dists[edge]
        #print(val)

        if val >= 0.05:
            style = '--'
            color = 'black'
            width = 1
        else:
            style = '-'
            if val < 0.001:
                width = 4
                color = 'red'
            elif val < 0.01:
                width = 2
                color = 'red'
            else:
                width = 1
                color = 'black'

        #edge_labels[edge] = str(round(1-dists[edge], 4))
        edge_width.append(width)
        edge_color.append(color)
        edge_style.append(style)
    
    pos = nx.kamada_kawai_layout(G)
    
    nx.draw_networkx_nodes(G\
                         , pos=pos\
                         , nodelist=['0', '1', '2', '3']\
                         , node_size=1000\
                         , node_color='#1A85FF'\
                         )

    #print(pos)
    nx.draw_networkx_labels(G\
                            , pos=pos\
                            , labels = {
                                '0':'0'\
                                ,'1':'1'\
                                ,'2':'2'\
                                ,'3':'3'\
                                }\
                            , font_color='white'\
                            , font_weight='heavy'\
                            , font_size=18\
                            )
    
    
    nx.draw_networkx_edges(G\
                         , pos=pos\
                         , edgelist = edge_list\
                         , edge_color = edge_color\
                         , width = edge_width\
                         , style = edge_style\
                         )
    
    #nx.draw_networkx_edge_labels(G\
    #                            , pos=pos\
    #                            , edge_labels=edge_labels\
    #                            )
    

    
    plt.grid(None)
    plt.gca().set_facecolor("white")
    #plt.savefig(title)
    #plt.show()
    #exit()
    plt.savefig('figures/max_pers_dev_pvals_'+typ+'.pdf')

    plt.cla()
    plt.clf()

    return ptiles_dict


def plot_all_pers_distri_dev(PD_info, fname, typ, clip_below, clip_above):
    

    dev_stages = [\
            'Gestation/'\
            ,'1-35weeks/'\
            ,'12-24months/'
            ,'28MandAfter/'\
                    ]

    all_pers = dict()

    #progress = 0

    for key in PD_info:

        #print(progress)
        #progress += 1

        stage = key.split('_')[0][-1]

        this_typ = key.split('_')[-1]

        if typ != 'all':
            if this_typ != typ:
                continue

        if stage not in all_pers:
            all_pers[stage] = []

        target = 'PD_results/' + key + '_H1_pers_data.txt'

        if not os.path.isfile(target):
            continue
        ff = np.loadtxt(target, delimiter=',')
        if len(ff) == 0:
            continue
        if ff.ndim == 1:
            ff = np.reshape(ff, (1, 2))

        this_pers = ff[:,1] - ff[:,0]

        for pers in this_pers:
            all_pers[stage].append(pers)


    #sns.kdeplot(data=max_pers, fill=True, clip=(0,trunc))
    #sns.displot(data=vals, kind='kde', hue=''fill=True)

    #sns.ecdfplot(data=max_pers['0'], label='0', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['1'], label='1', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['2'], label='2', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['3'], label='3', complementary=True, log_scale=2)

    sns.set(font_scale=1.35)

    sns.kdeplot(data=all_pers['0'], label=dev_stages[0], log_scale=2, clip=(clip_below, clip_above))
    sns.kdeplot(data=all_pers['1'], label=dev_stages[1], log_scale=2, clip=(clip_below, clip_above))
    sns.kdeplot(data=all_pers['2'], label=dev_stages[2], log_scale=2, clip=(clip_below, clip_above))
    sns.kdeplot(data=all_pers['3'], label=dev_stages[3], log_scale=2, clip=(clip_below, clip_above))


    #title =  '0:'+str(len(max_pers['0']))\
    #        ,', 1:'+str(len(max_pers['1']))\
    #        ,', 2:'+str(len(max_pers['2']))\
    #        ,', 3:'+str(len(max_pers['3']))

    plt.xlabel('persistence')

    plt.legend()

    plt.savefig(fname)

    plt.cla()
    plt.clf()

def plot_n_pers_distri_dev(PD_info, fname, typ, thresh):
    
    dev_stages = [\
            'Gestation/'\
            ,'1-35weeks/'\
            ,'12-24months/'
            ,'28MandAfter/'\
                    ]

    all_pers = dict()

    all_stages = []
    all_pers = []

    progress = 0

    for key in PD_info:

        #print(key)
        if progress % 1000 == 0:
            print(progress)

        progress += 1

        stage = key.split('_')[0][-1]

        this_typ = key.split('_')[-1]

        if this_typ == 'ad':
            islet_id = key[:-3]
        else:
            islet_id = key[:-2]

        fprefix = islet_id

        ###############################
        ###### Get G(r) thresh ########
        ###############################
        thresh_file = 'Thresholds/'+ fprefix + '.smooth2.thresh.minBetPeak2And3.dat'
        threshs_data = open(thresh_file, 'r')
        threshs = threshs_data.readline()
        threshs = threshs.strip('\n')
        threshs = threshs.split('\t')
        ad_thresh = float(threshs[1])
        b_thresh = float(threshs[2])
        threshs_data.close()
        ###############################

        if this_typ == 'ad':
            this_thresh = ad_thresh
        else:
            this_thresh = b_thresh

        #print(islet_id)
        #print(ad_thresh, b_thresh)

        if typ != 'all':
            if this_typ != typ:
                continue

        #if stage not in all_pers:
        #    all_pers[stage] = []

        target = 'PD_results/' + key + '_H1_pers_data.txt'


        if not os.path.isfile(target):
            continue

        ####################################
        ### CHECK PLOT ###
        ####################################
        #check_v_file = 'PD_results/' + key + '.csv'
        #check_v = np.loadtxt(check_v_file, delimiter=',')
        #plt.scatter(check_v[:,0], check_v[:,1])
        #plt.show()
        #print(check_v)
        #exit()
        ####################################

        ff = np.loadtxt(target, delimiter=',')
        if len(ff) == 0:
            continue
        if ff.ndim == 1:
            ff = np.reshape(ff, (1, 2))

        #print(ff)

        idxs = np.argwhere(ff[:,0] <= this_thresh).flatten()

        if not len(idxs):
            continue

        ff = ff[idxs]

        this_pers = ff[:,1] - ff[:,0]

        all_stages += len(this_pers) * [stage]
        all_pers += list(this_pers)

        #print(all_stages)
        #print(all_pers)
        #input('w')

        #print(this_pers)
        #exit()

        #count = len(np.argwhere(this_pers >= thresh).flatten())

        #all_pers[stage].append(count)

        #all_stages += len(this_pers)*[stage]
        #all_pers += list(this_pers)
        #print(all_stages)
        #print(all_pers)
        #exit()
        




    #sns.kdeplot(data=max_pers, fill=True, clip=(0,trunc))
    #sns.displot(data=vals, kind='kde', hue=''fill=True)

    #sns.ecdfplot(data=max_pers['0'], label='0', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['1'], label='1', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['2'], label='2', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['3'], label='3', complementary=True, log_scale=2)

    #sns.set(font_scale=1.35)

    #for stage in all_pers:
    #    vals = np.array(all_pers[stage])
    #    vals = np.log(1+vals)
    #    all_pers[stage] = vals

    #sns.kdeplot(data=all_pers['0'], label=dev_stages[0])
    #sns.kdeplot(data=all_pers['1'], label=dev_stages[1])
    #sns.kdeplot(data=all_pers['2'], label=dev_stages[2])
    #sns.kdeplot(data=all_pers['3'], label=dev_stages[3])


    d = {'Stages': all_stages, 'Persistence': all_pers}
    df = pd.DataFrame(data=d)

    sns.set(font_scale=1.5)

    sns.boxplot(data = df, x='Stages', y='Persistence')

    plt.show()

    #plt.savefig(fname, dpi=600)

    plt.cla()
    plt.clf()

    #plt.xscale('log', base = 2)
    ##title =  '0:'+str(len(max_pers['0']))\
    ##        ,', 1:'+str(len(max_pers['1']))\
    ##        ,', 2:'+str(len(max_pers['2']))\
    ##        ,', 3:'+str(len(max_pers['3']))

    #plt.xlabel('number of significant features')

    #plt.legend()

    #plt.savefig(fname)

    #plt.cla()
    #plt.clf()


def plot_box_significance(data, x, y, title, fname):

    sns.set(font_scale=1.5)

    sns.set_style("whitegrid")
    plt.tight_layout()

    sns.boxplot(data = data, x=x, y=y, palette=['lightblue', 'pink'])

    plt.yscale('log', base=2)

    plt.title(title)

    plt.savefig(fname, dpi=600)
    #plt.show()
    #exit()

    plt.cla()
    plt.clf()


def plot_pval_graph_dev(G, title):


    dists = nx.get_edge_attributes(G, "weight")
    edge_labels = dict()
    edge_list = []
    edge_style = []
    edge_color = []
    edge_width = []
    
    for edge in dists:

        edge_list.append(edge)

        val = 1-dists[edge]
        print(val)

        if val >= 0.05:
            style = '--'
            color = 'black'
            width = 1
        else:
            style = '-'
            if val < 0.001:
                width = 4
                color = 'red'
            elif val < 0.01:
                width = 2
                color = 'red'
            else:
                width = 1
                color = 'black'

        #edge_labels[edge] = str(round(1-dists[edge], 4))
        edge_width.append(width)
        edge_color.append(color)
        edge_style.append(style)
    
    pos = nx.kamada_kawai_layout(G)
    
    nx.draw_networkx_nodes(G\
                         , pos=pos\
                         , nodelist=['0', '1', '2', '3']\
                         , node_size=700\
                         , node_color='#1A85FF'\
                         )

    nx.draw_networkx_nodes(G\
                         , pos=pos\
                         , nodelist=['0_PH', '1_PH', '2_PH', '3_PH']\
                         , node_size=700\
                         , node_color='#D41159'\
                         #, linewidths = 8\
                         #, edgecolors='tab:orange'\
                         )
    
    
    #print(pos)
    nx.draw_networkx_labels(G\
                            , pos=pos\
                            , labels = {
                                 '0_PH':'0'\
                                ,'1_PH':'1'\
                                ,'2_PH':'2'\
                                ,'3_PH':'3'\
                                ,'0':'0'\
                                ,'1':'1'\
                                ,'2':'2'\
                                ,'3':'3'\
                                }
                            , font_color='white'\
                            , font_weight='heavy'\
                            , font_size=14\
                            )
    
    
    nx.draw_networkx_edges(G\
                         , pos=pos\
                         , edgelist = edge_list\
                         , edge_color = edge_color\
                         , width = edge_width\
                         , style = edge_style\
                         )
    
    #nx.draw_networkx_edge_labels(G\
    #                            , pos=pos\
    #                            , edge_labels=edge_labels\
    #                            )
    

    
    plt.grid(None)
    plt.gca().set_facecolor("white")
    plt.savefig(title)
    #plt.show()
    #exit()

    plt.cla()
    plt.clf()


def plot_pval_graph_diab(G, title):


    dists = nx.get_edge_attributes(G, "weight")
    edge_labels = dict()
    edge_list = []
    edge_style = []
    edge_color = []
    edge_width = []
    
    for edge in dists:

        edge_list.append(edge)

        val = 1-dists[edge]
        print(val)

        if val >= 0.05:
            style = '--'
            color = 'black'
            width = 1
        else:
            style = '-'
            if val < 0.001:
                width = 4
                color = 'red'
            elif val < 0.01:
                width = 2
                color = 'red'
            else:
                width = 1
                color = 'black'

        #edge_labels[edge] = str(round(1-dists[edge], 4))
        edge_width.append(width)
        edge_color.append(color)
        edge_style.append(style)
    
    pos = nx.kamada_kawai_layout(G)
    
    nx.draw_networkx_nodes(G\
                         , pos=pos\
                         , nodelist=['C', 'D']\
                         , node_size=700\
                         , node_color='#1A85FF'\
                         )
    
    nx.draw_networkx_nodes(G\
                         , pos=pos\
                         , nodelist=['C_PH', 'D_PH']\
                         , node_size=700\
                         , node_color='#D41159'\
                         )
    
    nx.draw_networkx_labels(G\
                            , pos=pos\
                            , labels = {
                                 'C_PH':'C'\
                                ,'D_PH':'D'\
                                ,'C':'C'\
                                ,'D':'D'\
                                }
                            , font_color='white'\
                            , font_weight='heavy'\
                            , font_size=14\
                            )
    
    
    nx.draw_networkx_edges(G\
                         , pos=pos\
                         , edgelist = edge_list\
                         , edge_color = edge_color\
                         , width = edge_width\
                         , style = edge_style\
                         )
    
    #nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)

    plt.grid(None)
    plt.gca().set_facecolor("white")
    plt.savefig(title)
    
    #plt.show()
    plt.cla()
    plt.clf()


def get_mannwhitney(d1, d2, axis):

    #d1 = d1[:,1]
    #d2 = d2[:,1]

    return scipy.stats.mannwhitneyu(d1, d2, alternative='two-sided', axis=axis )






def get_all_pers_dev(PD_info, fname):
    
    dev_stages = [\
            'Gestation/'\
            ,'1-35weeks/'\
            ,'12-24months/'
            ,'28MandAfter/'\
                    ]

    all_info = dict()

    #all_stages = []
    #all_pers = []
    #all_typ = []

    progress = 0

    for key in PD_info:

        #print(key)
        if progress % 1000 == 0:
            print(progress)

        progress += 1

        stage = key.split('_')[0][-1]

        this_typ = key.split('_')[-1]

        if this_typ == 'ad':
            islet_id = key[:-3]
        else:
            islet_id = key[:-2]

        fprefix = islet_id

        ################################
        ####### Get G(r) thresh ########
        ################################
        #thresh_file = 'Thresholds/'+ fprefix + '.smooth2.thresh.minBetPeak2And3.dat'
        #threshs_data = open(thresh_file, 'r')
        #threshs = threshs_data.readline()
        #threshs = threshs.strip('\n')
        #threshs = threshs.split('\t')
        #ad_thresh = float(threshs[1])
        #b_thresh = float(threshs[2])
        #threshs_data.close()

        #if this_typ == 'ad':
        #    this_thresh = ad_thresh
        #else:
        #    this_thresh = b_thresh
        ################################

        #print(islet_id)
        #print(ad_thresh, b_thresh)

        #if typ != 'all':
        #    if this_typ != typ:
        #        continue

        #if stage not in all_pers:
        #    all_pers[stage] = []

        target = 'PD_results/' + key + '_H1_pers_data.txt'


        if not os.path.isfile(target):
            continue

        ####################################
        ### CHECK PLOT ###
        ####################################
        #check_v_file = 'PD_results/' + key + '.csv'
        #check_v = np.loadtxt(check_v_file, delimiter=',')
        #plt.scatter(check_v[:,0], check_v[:,1])
        #plt.show()
        #print(check_v)
        #exit()
        ####################################

        ff = np.loadtxt(target, delimiter=',')
        if len(ff) == 0:
            continue
        if ff.ndim == 1:
            ff = np.reshape(ff, (1, 2))

        #print(ff)

        #idxs = np.argwhere(ff[:,0] <= this_thresh).flatten()

        #if not len(idxs):
        #    continue

        #ff = ff[idxs]

        this_pers = ff[:,1] - ff[:,0]

        #all_stages += len(this_pers) * [stage]
        #all_pers += list(this_pers)
        #all_typ += len(this_pers) * [this_typ]

        if stage not in all_info:
            all_info[stage] = dict()

        if this_typ not in all_info[stage]:
            all_info[stage][this_typ] = []

        #all_info[stage][this_typ] += list(this_pers)
        all_info[stage][this_typ] += list(ff)

        #all_info[stage][this_typ].append(max(list(this_pers)))

        #print(all_stages)
        #print(all_pers)
        #input('w')

        #print(this_pers)
        #exit()

        #count = len(np.argwhere(this_pers >= thresh).flatten())

        #all_pers[stage].append(count)

        #all_stages += len(this_pers)*[stage]
        #all_pers += list(this_pers)
        #print(all_stages)
        #print(all_pers)
        #exit()
        

    #all_info['all_stages'] = all_stages
    #all_info['all_pers'] = all_pers
    #all_info['all_typ'] = all_typ

    pickle.dump(all_info, open(fname, 'wb'))



    #sns.kdeplot(data=max_pers, fill=True, clip=(0,trunc))
    #sns.displot(data=vals, kind='kde', hue=''fill=True)

    #sns.ecdfplot(data=max_pers['0'], label='0', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['1'], label='1', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['2'], label='2', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['3'], label='3', complementary=True, log_scale=2)

    #sns.set(font_scale=1.35)

    #for stage in all_pers:
    #    vals = np.array(all_pers[stage])
    #    vals = np.log(1+vals)
    #    all_pers[stage] = vals

    #sns.kdeplot(data=all_pers['0'], label=dev_stages[0])
    #sns.kdeplot(data=all_pers['1'], label=dev_stages[1])
    #sns.kdeplot(data=all_pers['2'], label=dev_stages[2])
    #sns.kdeplot(data=all_pers['3'], label=dev_stages[3])


    #d = {'Stages': all_stages, 'Persistence': all_pers}
    #df = pd.DataFrame(data=d)

    #sns.set(font_scale=1.5)

    #sns.boxplot(data = df, x='Stages', y='Persistence')

    #plt.show()

    ##plt.savefig(fname, dpi=600)

    #plt.cla()
    #plt.clf()

    #plt.xscale('log', base = 2)
    ##title =  '0:'+str(len(max_pers['0']))\
    ##        ,', 1:'+str(len(max_pers['1']))\
    ##        ,', 2:'+str(len(max_pers['2']))\
    ##        ,', 3:'+str(len(max_pers['3']))

    #plt.xlabel('number of significant features')

    #plt.legend()

    #plt.savefig(fname)

    #plt.cla()
    #plt.clf()

def plot_all_pers_dev(fname):

    info = pickle.load(open(fname, 'rb'))

    fig, axs = plt.subplots(2, 4)

    row = 0
    col = 0

    for typ in ['ad', 'b']:

        for stage in info.keys():

            data = np.array(info[stage][typ])

            axs[row, col].hist2d(data[:,0]\
                    , data[:,1], bins=50, norm=mpl.colors.LogNorm(), cmap=mpl.cm.gnuplot)

            col += 1

        row += 1
        col = 0

    plt.show()
    exit()


    plt.xlabel('size diff', fontsize=14)
    plt.ylabel('bottleneck dist', fontsize=14)
    plt.colorbar()
    plt.savefig('figures/dev_11_PD_'+typ+'.pdf')
    plt.cla()
    plt.clf()

    stages = info.keys()
    types = ['ad', 'b']

    all_cats = list(it.product(stages, types))
    #print(all_cats)

    G = nx.Graph()

    for cat1, cat2 in it.combinations(all_cats, 2):

        stage1 = cat1[0]
        typ1 = cat1[1]
        vals1 = info[stage1][typ1]
        

        stage2 = cat2[0]
        typ2 = cat2[1]
        vals2 = info[stage2][typ2]

        res = get_mannwhitney(vals1, vals2, axis=0)

        G.add_edge(cat1, cat2, weight=1-res.pvalue)

    pos = nx.kamada_kawai_layout(G)
    

    #nx.draw_networkx_nodes(G\
    #                     , pos=pos\
    #                     , nodelist=['0_PH', '1_PH', '2_PH', '3_PH']\
    #                     , node_size=700\
    #                     , node_color='tab:orange'\
    #                     , alpha=0.5
    #                     )
    #
    #nx.draw_networkx_nodes(G\
    #                     , pos=pos\
    #                     , nodelist=['0', '1', '2', '3']\
    #                     , node_size=700\
    #                     , node_color='tab:blue'\
    #                     , alpha=0.5\
    #                     )
    
    nx.draw_networkx_labels(G\
                            , pos=pos\
                            , font_weight='heavy'\
                            )
    
    
    nx.draw_networkx_edges(G\
                         , pos=pos\
                         )

    dists = nx.get_edge_attributes(G, "weight")
    edge_labels  = dict()
    for edge in dists:
        edge_labels[edge] = str(round(1-dists[edge], 4))
    
    nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)

    #plt.savefig(title)
    
    plt.show()
    plt.cla()
    plt.clf()






def plot_max_pers_distri_diab(PD_info, fname, typ):
    
    print('Inside plot')
    plt.figure(figsize=(10,8))

    dev_stages = [\
            'Control'\
            ,'Diabetic'\
                    ]

    cats = ['C', 'D']

    max_pers_dict = dict()

    stages = []
    max_pers = []

    #progress = 0

    for key in PD_info:

        #print(progress)
        #progress += 1

        stage = key.split('_')[0][-1]

        this_typ = key.split('_')[-1]

        if typ != 'all':
            if this_typ != typ:
                continue


        if stage not in max_pers_dict:
            max_pers_dict[stage] = []

        max_pers_dict[stage].append(PD_info[key]['max_pers'])

        stages.append(stage)
        max_pers.append(PD_info[key]['max_pers'])

    ptiles_dict = dict()

    for cat in cats:
        vals = max_pers_dict[cat]
        vals = np.array(vals)
        ptiles = np.percentile(vals, [5, 50, 95])
        ptiles_dict[cat] = ptiles

    d = {'Stages': stages, 'Max persistence': max_pers}
    df = pd.DataFrame(data=d)

    sns.set(font_scale=3)
    
    sns.set_style("whitegrid")

    sns.boxplot(data = df, x='Stages', y='Max persistence', palette=['lightblue', 'pink'])
    plt.tight_layout()
    plt.yscale('log', base=2)

    print('Computing mann whitney')


    s1 = 'C'
    s2 = 'D'
    res = scipy.stats.mannwhitneyu(max_pers_dict[s1], max_pers_dict[s2], alternative='two-sided', axis=0 )
    print(s1, s2, res.pvalue)
    title = f"p-value = {res.pvalue:.1E}"


    #sns.kdeplot(data=max_pers, fill=True, clip=(0,trunc))
    #sns.displot(data=vals, kind='kde', hue=''fill=True)

    #sns.ecdfplot(data=max_pers['0'], label='0', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['1'], label='1', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['2'], label='2', complementary=True, log_scale=2)
    #sns.ecdfplot(data=max_pers['3'], label='3', complementary=True, log_scale=2)

    #sns.set(font_scale=1.35)

    #sns.kdeplot(data=max_pers['0'], label=dev_stages[0], log_scale=2, clip=(clip_below, clip_above))
    #sns.kdeplot(data=max_pers['1'], label=dev_stages[1], log_scale=2, clip=(clip_below, clip_above))
    #sns.kdeplot(data=max_pers['2'], label=dev_stages[2], log_scale=2, clip=(clip_below, clip_above))
    #sns.kdeplot(data=max_pers['3'], label=dev_stages[3], log_scale=2, clip=(clip_below, clip_above))






    #title =  '0:'+str(len(max_pers['0']))\
    #        ,', 1:'+str(len(max_pers['1']))\
    #        ,', 2:'+str(len(max_pers['2']))\
    #        ,', 3:'+str(len(max_pers['3']))

    #plt.xlabel('max persistence')

    
    #plt.legend()

    plt.title(title)
    plt.savefig(fname)
    #plt.show()

    plt.cla()
    plt.clf()

    return ptiles_dict


def plot_pval_graph(G, nodelist1, color1, nodelist2, color2, edge_labels, labels, fname):

    
    pos = nx.kamada_kawai_layout(G)
    #pos = nx.circular_layout(G)
    
    nx.draw_networkx_nodes(G\
                         , pos=pos\
                         , nodelist=nodelist1\
                         , node_size=700\
                         , node_color=color1\
                         )

    nx.draw_networkx_nodes(G\
                         , pos=pos\
                         , nodelist=nodelist2\
                         , node_size=700\
                         , node_color=color2\
                         , node_shape= 's'
                         #, linewidths = 8\
                         #, edgecolors='tab:orange'\
                         )
    
    
    #print(pos)
    nx.draw_networkx_labels(G\
                            , pos=pos\
                            , labels = labels
                            #{
                            #     '0_PH':'0'\
                            #    ,'1_PH':'1'\
                            #    ,'2_PH':'2'\
                            #    ,'3_PH':'3'\
                            #    ,'0':'0'\
                            #    ,'1':'1'\
                            #    ,'2':'2'\
                            #    ,'3':'3'\
                            #    }
                            , font_color='white'\
                            , font_weight='heavy'\
                            , font_size=14\
                            )
    
    
    nx.draw_networkx_edges(G\
                         , pos=pos\
                         )
    
    #nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)

    plt.savefig(fname)
    
    #plt.show()
    #exit()

    plt.cla()
    plt.clf()



def get_pairwise_pval_graph(info):


    G = nx.Graph()

    for s1, s2 in it.combinations(list(info.keys()), 2):

        #print(s1, s2)
    
        val1 = info[s1]
        val2 = info[s2]

        res = get_mannwhitney(val1, val2, axis=0)

        if res.pvalue < 0.05:
            G.add_edge(s1, s2, weight = 1-res.pvalue)
        else:
            G.add_node(s1)
            G.add_node(s2)

    return G


def get_pairwise_energydist_graph(info, fname, nodelist1, color1, nodelist2, color2, labels):


    G = nx.Graph()

    energy_dist_info = dict()

    types = ['CC', 'CD', 'DD']

    all_types = []
    all_vals = []
    dist_info = dict()

    for s1, s2 in it.combinations(list(info.keys()), 2):
    
        val1 = info[s1]
        val2 = info[s2]


        typ1 = s1[0]
        typ2 = s2[0]

        if typ1 != typ2:
            this_type = 'CD'
            continue
        elif typ1 == 'C':
            this_type = 'CC'
        else:
            this_type = 'DD'

        val = scipy.stats.energy_distance(val1, val2)

        if this_type not in dist_info:
            dist_info[this_type] = [val]
        else:
            dist_info[this_type].append(val)

        all_types.append(this_type)
        all_vals.append(val)

        if s1 not in energy_dist_info:
            energy_dist_info[s1] = [val]
        else:
            energy_dist_info[s1].append(val)

        if s2 not in energy_dist_info:
            energy_dist_info[s2] = [val]
        else:
            energy_dist_info[s2].append(val)

        G.add_edge(s1, s2, weight = val)
        
    d = {'Comparison': all_types, 'Energy dist': all_vals}
    df = pd.DataFrame(data=d)

    this_title = ''

    for t1, t2 in it.combinations(dist_info.keys(), 2):

        #res = get_mannwhitney(dist_info[t1], dist_info[t2], axis=0)
        res = scipy.stats.kstest(dist_info[t1], dist_info[t2])

        print(t1, t2, res.pvalue)

        this_title = 'pval(' + t1 + ',' + t2 + '):' + str(round(res.pvalue, 4))

        print(variance(dist_info[t1]), variance(dist_info[t2]))

    sns.boxplot(data = df, x='Comparison', y='Energy dist')
    sns.stripplot(data = df, x='Comparison', y='Energy dist', color='black')

    plt.title(this_title)

    #plt.show()

    plt.savefig(fname)
    plt.cla()
    plt.clf()

    energy_dist_pval_graph = nx.Graph()
    for s1, s2 in it.combinations(list(info.keys()), 2):
        
        typ1 = s1[0]
        typ2 = s2[0]

        if typ1 != typ2:
            continue

        val1 = energy_dist_info[s1]
        val2 = energy_dist_info[s2]

        res = get_mannwhitney(val1, val2, axis=0)

        if res.pvalue < 0.05:
            energy_dist_pval_graph.add_edge(s1, s2, weight=1-res.pvalue)
        else:
            energy_dist_pval_graph.add_node(s1)
            energy_dist_pval_graph.add_node(s2)




    dists = nx.get_edge_attributes(energy_dist_pval_graph, "weight")
    edge_labels = dict()
    for edge in dists:
        edge_labels[edge] = str(round(1-dists[edge], 4))

    plot_pval_graph(energy_dist_pval_graph\
                    , nodelist1, color1\
                    , nodelist2, color2\
                    , edge_labels, labels\
                    ,fname[:-4]+'_pvalgraph.pdf')


    return G





