#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
#
##############################################################################
# Script used in the paper:
# Evaluation of TanDEM-X DEMs on selected Brazilian sites: comparison with SRTM, ASTER GDEM and ALOS AW3D30
# by
# Carlos H. Grohmann - 2017
# guano (at) usp (dot) br
# Institute of Energy and Environment - University of Sao Paulo
##############################################################################


# This file has the pre-analysis
# test for SRTM how many random points are needed to characterize each study area. 

# 1 - Extract elevation for sets of randon points [50,100,250,500,1000,2000,5000,10000] and then
# calculate stats for each. 

# 2 - MonteCarlo-like analysis:
# 0 - for a series of n random points [50,100,250,500,1000,2000,5000,10000]:
#     1 - get X sets of n random points (rand_n_01, rand_n_02, rand_n_03,...) - sorted
#     2 - calculate correlation between first set and all others
#     3 - put data in a table and plot the results
#     4 - make plot off all values (X=n_points, Y=correlation)


import sys, os, csv
import itertools
import numpy as np
import scipy.stats as ss
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import gc

import grass.script as grass
import grass.script.array as garray

# import homemade functions
funcDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/scripts/'
os.chdir(funcDir)
from functions import *

#----------------------------------------------------
#----------------------------------------------------
# helper func: delete all *random maps
def clean_rand():
    grass.run_command('g.remove', type='vector', pattern='*random*', flags='f')


#----------------------------------------------------
# Stats for random points
#----------------------------------------------------
# list of srtm dems
dem_list = ['srtm30_Araca','srtm30_Barcelos','srtm30_Pantanal',\
    'srtm30_SerraMar','srtm30_RioClaro','srtm30_Iporanga','srtm30_SantaCatarina']

# one list for all dataframes
df_all = []

# get stats for each srtm tile (for each area), and join stats in dataframe
for dem in dem_list:
    print ('')
    print ('')
    print ('----------  '+dem+'  ----------')
    grass.run_command('g.region', raster=dem, flags='a')
    # dataframe to hold stats
    df = pd.DataFrame()
    df['stats'] = ['min','max','range','mean','median','stddev','skew','kurt','p25','p75']
    # get raster data at full resolution (30m)
    reference = flatParam(dem)
    # stats for reference data 
    df[dem] = calc_stats(reference)
    # random points + stats
    for points in [50,100,250,500,1000,2000,5000,10000]:
        col_name = 'rand_' + str(points)
        elev = sample_dems(dem, points, r_seed=1234, ow_vector=True, ow_what=True)
        df[col_name] = calc_stats(elev)
    # set index col
    df_indexed = df.set_index('stats')
    df_all.append(df_indexed)
    print ('----------  end  ----------')
    print ('')
    print ('')

# export tables of stats to latex
for df in df_all:
    print df.to_latex(float_format=lambda x:'%10.3f' % x)
    print('\n\n')


#----------------------------------------------------
# MonteCarlo-like analysis
#----------------------------------------------------
workDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/rand_mc/'
os.chdir(workDir)

n_random = 50
dem = 'srtm30_Araca'
npoints_list = [50,100,250,500,1000,2500,5000,10000]

# list of srtm dems
dem_list = ['srtm30_Araca','srtm30_Barcelos','srtm30_Pantanal',\
    'srtm30_SerraMar','srtm30_RioClaro','srtm30_Iporanga','srtm30_SantaCatarina']

# run monte carlo sampling and save data as csv
for dem,points in itertools.product(dem_list, npoints_list):
    df = pd.DataFrame()
    file_out = dem + '_rand_MC_' + str(points)
    for run in range(n_random):
        col_name = 'rand_' + str(points) + '_' + str(run).zfill(2)
        elev = sample_dems_mc(dem, points, mc_run=run, ow_vector=True, ow_what=True)
        df[col_name] = np.sort(elev)
        df.to_csv(path_or_buf=file_out+'.csv', na_rep='NaN')


# reads data from csv files and calculates correlation
for dem in dem_list:
    avg_corr = []
    df_corr = pd.DataFrame()
    for points in npoints_list:
        csv_file = dem + '_rand_MC_' + str(points) + '.csv'
        print(csv_file)
        df = pd.read_csv(csv_file, index_col=0)
        # correlation of first column[0] with  all the others [1:]. 
        # No need to define column by name
        corr = df.corr().iloc[0,1:]
        avg_corr.append(corr.mean())
        # plot correlation values for this set of random points
        x_ax = np.empty(n_random -1)
        x_ax.fill(points)
        plt.plot(x_ax, corr, 'o')
        # add correlations (df.series) to dataframe
        df_corr[str(points)+'pts']=corr.values
    # export df to latex
    print dem, points
    print df_corr.to_latex(float_format=lambda x:'%4.4f' % x)

# curve fitting for mean of correlation values
# Initial guess for curve fitting parameters
p0 = [1, 1, 1, 1]
# observations
x = npoints_list
y_meas = avg_corr
# least-squares fitting
plsq = leastsq(residuals, p0, args=(y_meas, x)) 
equation = 'y = ((A-D)/(1.0+((x/C)**B))) + D' 
A = plsq[0][0]
B = plsq[0][1]
C = plsq[0][2]
D = plsq[0][3]

# sequence of values for curve plotting
xx=np.arange(0,10500,25)
# plot fitted curve
plt.plot(xx, peval(xx,plsq[0]), color='0.5', lw=0.9)
plt.plot(x, y_meas, 'o', color='0.5', ms=5)

# final plot decorations - values of (ymax,ymin) for each area were set manually
ymin_dict = {'srtm30_Araca':0.747,'srtm30_Barcelos':0.91,'srtm30_Pantanal':0.962,\
    'srtm30_SerraMar':0.857,'srtm30_RioClaro':0.937,'srtm30_Iporanga':0.907,'srtm30_SantaCatarina':0.797}

x_ticks = npoints_list
x_labels = ['50','100','250','500','1000','2500','5000','10000']
plt.xlim((40,12500))
ymax = 1.003
ymin = ymin_dict[dem]
yrange = ymax - ymin
plt.ylim((ymin,ymax))
plt.xscale('log')
plt.xticks(x_ticks, x_labels)
plt.xlabel('Number of random points')
plt.ylabel("Pearson's correlation coefficient")
plt.text(550,ymin+yrange*0.60,'4PL curve fit', fontsize=11)
plt.text(550,ymin+yrange*0.55, equation, fontsize=11)
plt.text(550,ymin+yrange*0.50, 'A = ' +str(A), fontsize=11)
plt.text(550,ymin+yrange*0.45, 'B = ' +str(B), fontsize=11)
plt.text(550,ymin+yrange*0.40, 'C = ' +str(C), fontsize=11)
plt.text(550,ymin+yrange*0.35, 'D = ' +str(D), fontsize=11)
plt.tight_layout()
# plt.show()

# save final figure
plt.savefig(dem + '_rand_MC.svg')
plt.cla()
plt.clf()


#----------------------------------------------------
# adds mean values of correlations to dataframe so we can export the table to latex
dem_corr_4pd = []
cols = ['50','100','250','500','1000','2500','5000','10000']

for dem in dem_list:
    dem_corr = []
    # df_corr = pd.DataFrame()
    for points in npoints_list:
        npoints_names.append('avg_corr_'+str(points))
        csv_file = dem + '_rand_MC_' + str(points) + '.csv'
        print(csv_file)
        df = pd.read_csv(csv_file, index_col=0)
        # correlation of first column[0] with  all the others [1:]. 
        # No need to define column by name
        corr = df.corr().iloc[0,1:]
        dem_corr.append(corr.mean())
    row_name = dem.split('_')[1]
    dem_corr_4pd.append((row_name,dem_corr))

pd_corrs = pd.DataFrame.from_items(dem_corr_4pd, columns=cols, orient='index')

print(pd_corrs.to_latex(float_format=lambda x:'%3.6f' % x))


#----------------------------------------------------
# adds standard deviations values of correlations to dataframe so we can export the table to latex
dem_std_4pd = []
cols = ['50','100','250','500','1000','2500','5000','10000']

for dem in dem_list:
    dem_std = []
    # df_corr = pd.DataFrame()
    for points in npoints_list:
        # npoints_names.append('avg_corr_'+str(points))
        csv_file = dem + '_rand_MC_' + str(points) + '.csv'
        print(csv_file)
        df = pd.read_csv(csv_file, index_col=0)
        # correlation of first column[0] with  all the others [1:]. 
        # No need to define column by name
        corr = df.corr().iloc[0,1:]
        dem_std.append(corr.std())
    row_name = dem.split('_')[1]
    dem_std_4pd.append((row_name,dem_std))

pd_corrs = pd.DataFrame.from_items(dem_std_4pd, columns=cols, orient='index')

print(pd_corrs.to_latex(float_format=lambda x:'%3.6f' % x))


