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


# This file has the main analysis
# it uses the functions defined in functions.py

import sys, os, itertools
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns

try:
    import grass.script as grass
    import grass.script.array as garray
except:
    pass

#----------------------------------------------------
# import homemade functions
funcDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/scripts/'
# funcDir = '/Users/guano/Dropbox/USP/projetosPesquisa/TanDEM-X/scripts'
os.chdir(funcDir)
from functions import *


#----------------------------------------------------
# WGS84 x EGM96
# stats of DEMs before and after converting elevation to WGS84 ellipsoid
#----------------------------------------------------

#----------------------------------------------------
# set working directory
work_dir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/rand_10k_csv/'
os.chdir(work_dir)

#----------------------------------------------------
# list of tiles, by geographical area (egm96 + wgs84)
araca =         ['tdx12_Araca',         'tdx30_Araca',         'srtm30_Araca',         'srtm30_wgs84_Araca',         'aster30_Araca',         'aster30_wgs84_Araca',         'aw3d30_Araca',         'aw3d30_wgs84_Araca',          ]
barcelos =      ['tdx12_Barcelos',      'tdx30_Barcelos',      'srtm30_Barcelos',      'srtm30_wgs84_Barcelos',      'aster30_Barcelos',      'aster30_wgs84_Barcelos',      'aw3d30_Barcelos',      'aw3d30_wgs84_Barcelos',       ]
pantanal =      ['tdx12_Pantanal',      'tdx30_Pantanal',      'srtm30_Pantanal',      'srtm30_wgs84_Pantanal',      'aster30_Pantanal',      'aster30_wgs84_Pantanal',      'aw3d30_Pantanal',      'aw3d30_wgs84_Pantanal',       ]
serramar =      ['tdx12_SerraMar',      'tdx30_SerraMar',      'srtm30_SerraMar',      'srtm30_wgs84_SerraMar',      'aster30_SerraMar',      'aster30_wgs84_SerraMar',      'aw3d30_SerraMar',      'aw3d30_wgs84_SerraMar',       ]
rioclaro =      ['tdx12_RioClaro',      'tdx30_RioClaro',      'srtm30_RioClaro',      'srtm30_wgs84_RioClaro',      'aster30_RioClaro',      'aster30_wgs84_RioClaro',      'aw3d30_RioClaro',      'aw3d30_wgs84_RioClaro',       ]
iporanga =      ['tdx12_Iporanga',      'tdx30_Iporanga',      'srtm30_Iporanga',      'srtm30_wgs84_Iporanga',      'aster30_Iporanga',      'aster30_wgs84_Iporanga',      'aw3d30_Iporanga',      'aw3d30_wgs84_Iporanga',       ]
santacatarina = ['tdx12_SantaCatarina', 'tdx30_SantaCatarina', 'srtm30_SantaCatarina', 'srtm30_wgs84_SantaCatarina', 'aster30_SantaCatarina', 'aster30_wgs84_SantaCatarina', 'aw3d30_SantaCatarina', 'aw3d30_wgs84_SantaCatarina',  ]
# all areas
areas=[araca,barcelos,pantanal,iporanga,rioclaro,santacatarina,serramar]


#----------------------------------------------------
# Sample DEMs -  10,000 random points
points = 10000

#----------------------------------------------------
# sample dems and save data as csv files - one file per area
for area in areas:
    print(area)
    df = pd.DataFrame()
    file_out = area[0].split('_')[-1] + '_randon_' + str(points) + '.csv'
    for dem in area:
        col_name = dem#.split('_')[0]
        df['X'],df['Y'],df[col_name] = sample_dems(dem, points, r_seed=1234, ow_vector=False, ow_what=False, coords=True)
    df.to_csv(path_or_buf=file_out, na_rep='NaN')

#----------------------------------------------------
# reads data from csv files and calculates stats - display as latex tables
for area in areas:
    df_stats = pd.DataFrame()
    df_stats['stats'] = ['min','max','range','mean','median','stddev','skew','kurt','p25','p75']
    csv_file = area[0].split('_')[-1] + '_randon_' + str(points) + '.csv'
    print(csv_file)
    df = pd.read_csv(work_dir + csv_file, index_col=0)
    # runs stats only in columns with elev values
    for column in df.iloc[:,2:]:
        elev = df[column].values
        df_stats[column] = calc_stats(elev)
    # export stats to latex
    print('\n\n')
    # print df_stats.set_index('stats').to_latex(float_format=lambda x:'%4.4f' % x)
    print df_stats.transpose().to_latex(float_format=lambda x:'%4.2f' % x)


#----------------------------------------------------
#----------------------------------------------------
# WGS84
# actual analysis starts here. all elevations are WGS84
#----------------------------------------------------
#----------------------------------------------------

# set working directory
work_dir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/rand_10k_csv/'
os.chdir(work_dir)

#----------------------------------------------------
# list of tiles, by geographical area (all elevations are wgs84)
araca = ['tdx12_Araca', 'tdx30_Araca', 'srtm30_wgs84_Araca', 'aster30_wgs84_Araca', 'aw3d30_wgs84_Araca']
barcelos = ['tdx12_Barcelos', 'tdx30_Barcelos', 'srtm30_wgs84_Barcelos', 'aster30_wgs84_Barcelos', 'aw3d30_wgs84_Barcelos']
pantanal = ['tdx12_Pantanal', 'tdx30_Pantanal', 'srtm30_wgs84_Pantanal', 'aster30_wgs84_Pantanal', 'aw3d30_wgs84_Pantanal']
serramar = ['tdx12_SerraMar', 'tdx30_SerraMar', 'srtm30_wgs84_SerraMar', 'aster30_wgs84_SerraMar', 'aw3d30_wgs84_SerraMar']
rioclaro = ['tdx12_RioClaro', 'tdx30_RioClaro', 'srtm30_wgs84_RioClaro', 'aster30_wgs84_RioClaro', 'aw3d30_wgs84_RioClaro']
iporanga = ['tdx12_Iporanga', 'tdx30_Iporanga', 'srtm30_wgs84_Iporanga', 'aster30_wgs84_Iporanga', 'aw3d30_wgs84_Iporanga']
santacatarina = ['tdx12_SantaCatarina', 'tdx30_SantaCatarina', 'srtm30_wgs84_SantaCatarina', 'aster30_wgs84_SantaCatarina', 'aw3d30_wgs84_SantaCatarina']
# all areas
areas=[araca,barcelos,pantanal,iporanga,rioclaro,santacatarina,serramar]

#----------------------------------------------------
# calculate histograms, one area at a time, so we can set up limits and bins individually 
do_histogram_np(araca, work_dir, bin_width=2, x_lim=(-10,150), density_plot=True, points=points)
do_histogram_np(barcelos, work_dir, bin_width=2, x_lim=(-10,120), density_plot=True, points=points)
do_histogram_np(pantanal, work_dir, bin_width=2, x_lim=(50,200), density_plot=True, points=points)
do_histogram_np(serramar, work_dir, bin_width=20, x_lim=(-10,2700), density_plot=True, points=points)
do_histogram_np(rioclaro, work_dir, bin_width=10, x_lim=(400,1100), density_plot=True, points=points)
do_histogram_np(iporanga, work_dir, bin_width=20, x_lim=(-10,1250), density_plot=True, points=points)
do_histogram_np(santacatarina, work_dir, bin_width=2, x_lim=(-2,150), density_plot=True, points=points)


# full histograms (limits = min-max) for srtm, for each area (they will be insets on the hists for all dems)
do_histogram_full('srtm30_Araca', work_dir, bin_width=10, points=points)
do_histogram_full('srtm30_Barcelos', work_dir, bin_width=1, points=points)
do_histogram_full('srtm30_Pantanal', work_dir, bin_width=1, points=points)
do_histogram_full('srtm30_SerraMar', work_dir, bin_width=10, points=points)
do_histogram_full('srtm30_RioClaro', work_dir, bin_width=5, points=points)
do_histogram_full('srtm30_Iporanga', work_dir, bin_width=10, points=points)
do_histogram_full('srtm30_SantaCatarina', work_dir, bin_width=10, points=points)


#----------------------------------------------------
# Mean slope versus elevation
#----------------------------------------------------
# R
# sf <- readRAST6(c("elevation.dem", "slope"))
# df <- as(sf, "data.frame")
# df1 <- df[order(df$elevation.dem),]
# df1$elev.factor <- cut(df1$elevation.dem, breaks=seq(1051, 1840, 10),
#    ordered_result=TRUE)
# df2 <- tapply(df1$slope, df1$elev.factor, mean, na.rm=TRUE) -- pandas pivot_table
# str(df2)
# plot(df2, seq(1051, 1840, 10), type="l")

work_dir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/rand_10k_csv/'
os.chdir(work_dir)

# calculate slope
for area in areas:
    for dem in area:
        slope = dem.split('_')[0]+'_slope_'+area[0].split('_')[1]
        grass.run_command('g.region', raster=dem, flags='pa')
        grass.run_command("r.slope.aspect", elevation=dem, slope=slope, overwrite=True)

# sample DEMs - one CSV file per DEM
points = 100000
for area in areas:
    for dem in area:
        slope = dem.split('_')[0]+'_slope_'+area[0].split('_')[1]
        df = pd.DataFrame()
        file_out = dem + '_elev_w_slope_'+str(points)+'.csv'   
        # sample rasters
        print dem, slope
        df[dem] = sample_dems(dem, points, r_seed=1234, ow_vector=False, ow_what=True)
        df[slope] = sample_dems(slope, points, r_seed=1234, ow_vector=True, ow_what=True, flags='i')
        df.to_csv(path_or_buf=file_out, na_rep='NaN')

# Stats of slope
# reads data from csv files and calculates stats - display as latex tables
points = 10000
for area in areas:
    df_stats = pd.DataFrame()
    df_stats['stats'] = ['min','max','range','mean','median','stddev','skew','kurt','p25','p75']    
    for dem in area:
        slope = dem.split('_')[0]+'_slope_'+area[0].split('_')[1]
        csv_file = dem + '_elev_w_slope_'+str(points)+'pts.csv' 
        print(csv_file)
        df = pd.read_csv(csv_file, index_col=0)
        elev = df[slope].values
        df_stats[slope] = calc_stats(elev)
        # export stats to latex
        print('\n\n')
    print df_stats.transpose().to_latex(float_format=lambda x:'%4.2f' % x)

# calculate histograms, one area at a time, so we can set up limits and bins individually 
do_histogram_slope(araca         , work_dir, bin_width=1, x_lim=(0,90), density_plot=True, points=points, n_colors=5)
do_histogram_slope(barcelos      , work_dir, bin_width=1, x_lim=(0,90), density_plot=True, points=points, n_colors=5)
do_histogram_slope(pantanal      , work_dir, bin_width=1, x_lim=(0,90), density_plot=True, points=points, n_colors=5)
do_histogram_slope(serramar      , work_dir, bin_width=1, x_lim=(0,90), density_plot=True, points=points, n_colors=5)
do_histogram_slope(rioclaro      , work_dir, bin_width=1, x_lim=(0,90), density_plot=True, points=points, n_colors=5)
do_histogram_slope(iporanga      , work_dir, bin_width=1, x_lim=(0,90), density_plot=True, points=points, n_colors=5)
do_histogram_slope(santacatarina , work_dir, bin_width=1, x_lim=(0,90), density_plot=True, points=points, n_colors=5)

# plot mean-slope x elevation, WITHOUT adjusted slope limits
mean_slope_elev(area_list=araca         ,elev_interval=50 , zmin=0  , zmax=1750, smin=0, smax=90, points=points, work_dir=work_dir, suffix='full')
mean_slope_elev(area_list=barcelos      ,elev_interval=2  , zmin=0  , zmax=75  , smin=0, smax=90, points=points, work_dir=work_dir, suffix='full')
mean_slope_elev(area_list=pantanal      ,elev_interval=2  , zmin=90 , zmax=170 , smin=0, smax=90, points=points, work_dir=work_dir, suffix='full')
mean_slope_elev(area_list=iporanga      ,elev_interval=50 , zmin=0  , zmax=1250, smin=0, smax=90, points=points, work_dir=work_dir, suffix='full')
mean_slope_elev(area_list=rioclaro      ,elev_interval=50 , zmin=400, zmax=1050, smin=0, smax=90, points=points, work_dir=work_dir, suffix='full')
mean_slope_elev(area_list=santacatarina ,elev_interval=50 , zmin=0  , zmax=1300, smin=0, smax=90, points=points, work_dir=work_dir, suffix='full')
mean_slope_elev(area_list=serramar      ,elev_interval=100, zmin=0  , zmax=2700, smin=0, smax=90, points=points, work_dir=work_dir, suffix='full')

# plot mean-slope x elevation, with adjusted slope limits
mean_slope_elev(area_list=araca         ,elev_interval=50 , zmin=0  , zmax=1750, smin=0, smax=50, points=points, work_dir=work_dir, suffix='tight')
mean_slope_elev(area_list=barcelos      ,elev_interval=2  , zmin=0  , zmax=75  , smin=0, smax=45, points=points, work_dir=work_dir, suffix='tight')
mean_slope_elev(area_list=pantanal      ,elev_interval=2  , zmin=70 , zmax=170 , smin=0, smax=15, points=points, work_dir=work_dir, suffix='tight')
mean_slope_elev(area_list=iporanga      ,elev_interval=50 , zmin=0  , zmax=1250, smin=0, smax=55, points=points, work_dir=work_dir, suffix='tight')
mean_slope_elev(area_list=rioclaro      ,elev_interval=50 , zmin=400, zmax=1050, smin=0, smax=15, points=points, work_dir=work_dir, suffix='tight')
mean_slope_elev(area_list=santacatarina ,elev_interval=50 , zmin=0  , zmax=1300, smin=0, smax=25, points=points, work_dir=work_dir, suffix='tight')
mean_slope_elev(area_list=serramar      ,elev_interval=100, zmin=0  , zmax=2700, smin=0, smax=40, points=points, work_dir=work_dir, suffix='tight')

# export slope as pngs
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        slope = dem.split('_')[0]+'_slope_'+area_name
        grass.run_command('g.region', raster=slope, flags='pa', res='0:0:01')
        grass.run_command('d.mon', start='png', output=slope+'.png', resolution = '2', \
            height=500, width=500, overwrite=True)
        grass.run_command('d.rast', map=slope)
        grass.run_command('d.grid', size='0.25', text_color='black', flags='c')
        grass.run_command('d.legend', raster=slope, flags='sb', at='4,25,2,4', font='Helvetica', fontsize=14, bgcolor='240:240:240')
        grass.run_command('d.mon', stop='png')


#----------------------------------------------------
# Errors compared to SRTM (rmse, mae, me, std)
#----------------------------------------------------

# calculate error metrics
rmse_list, mae_list, me_list, std_list, dems, n_nans, area_list = calc_error_metrics(areas,work_dir,points)

# make list of lists and put all in a DataFrame
err4pd = [('DEM', dems),
         ('area', area_list),
         ('ME', me_list),
         ('MAE', mae_list),                  
         ('RMSE', rmse_list),
         ('STDE', std_list),
         ('n_NANs', n_nans),]
pd_err = pd.DataFrame.from_items(err4pd)

# export metrics to latex
print('\n\n')
print pd_err.sort_values(by='DEM').to_latex(float_format=lambda x:'%4.2f' % x)

# export stats of each metric to latex
# aster
print pd_err[pd_err['DEM'].str.startswith('aster')].describe().to_latex(float_format=lambda x:'%4.3f' % x)

# alos
print pd_err[pd_err['DEM'].str.startswith('aw3d')].describe().to_latex(float_format=lambda x:'%4.3f' % x)

# tdx12
print pd_err[pd_err['DEM'].str.startswith('tdx12')].describe().to_latex(float_format=lambda x:'%4.3f' % x)

# tdx30
print pd_err[pd_err['DEM'].str.startswith('tdx30')].describe().to_latex(float_format=lambda x:'%4.3f' % x)


#----------------------------------------------------
# Stats of differences - SRTM 30m minus all other DEMs
#----------------------------------------------------

# calculate differences (with data from the CSV files)
for area in areas:
    area_name = area[0].split('_')[-1]
    csv_file = area_name + '_randon_' + str(points) + '.csv'
    print(csv_file)
    # get elev data from csv files 
    df = pd.read_csv(work_dir + csv_file, index_col=0)
    # calculate differences from SRTM 30m
    df['srtm30_wgs84-tdx12'] = df['srtm30_wgs84_'+area_name]-df['tdx12_'+area_name]
    df['srtm30_wgs84-tdx30'] = df['srtm30_wgs84_'+area_name]-df['tdx30_'+area_name]
    df['srtm30_wgs84-aster30_wgs84'] = df['srtm30_wgs84_'+area_name]-df['aster30_wgs84_'+area_name]
    df['srtm30_wgs84-aw3d30_wgs84'] = df['srtm30_wgs84_'+area_name]-df['aw3d30_wgs84_'+area_name]
    # save dataframe with diffs as new csv file
    file_out = area[0].split('_')[-1] + '_randon_' + str(points) + '_with_diffs_from_srtm' + '.csv'
    df.to_csv(path_or_buf=work_dir+file_out, na_rep='NaN')

# make new datafram with all columns of differences, so we can have histograms grouped by DEM, not per area
df_diffs_per_dem = pd.DataFrame()
for area in areas:
    area_name = area[0].split('_')[-1]
    csv_file = area_name + '_randon_' + str(points) + '_with_diffs_from_srtm' + '.csv'
    print(csv_file)
    # get elev data from csv files 
    df = pd.read_csv(work_dir + csv_file, index_col=0)
    # headers = df.dtypes.index
    headers_diff = ['srtm30_wgs84-tdx12','srtm30_wgs84-tdx30','srtm30_wgs84-aster30_wgs84','srtm30_wgs84-aw3d30_wgs84']
    for column in df[headers_diff]:
        dem_name = column.split('-')[-1]
        if dem_name.endswith('wgs84'):
            dem_name = dem_name.split('_')[0]
        col_name = area_name + '_' + dem_name
        df_diffs_per_dem[col_name] = df[column].values

# export dataframe
df_diffs_per_dem.to_csv(path_or_buf='all_dems_diffs_from_srtm.csv', na_rep='NaN')


#----------------------------------------------------
# histograms of differences - per dem
csv_file = 'all_dems_diffs_from_srtm.csv'
df = pd.read_csv(work_dir + csv_file, index_col=0)
headers = df.dtypes.index

for dem in ['tdx12','tdx30','aster30','aw3d30']:
    file_out = dem + '_diffs_from_srtm_hist.pdf'
    sns.set_palette(sns.color_palette("muted", n_colors=7, desat=.5))
    for column in df[headers]:
        if column.endswith(dem):
            elev = df[column].values
            elev = elev[~np.isnan(elev)]
            # define bins
            g_max = int(np.ceil(np.max(elev)))
            g_min = int(np.min(elev))
            nbins = g_max - g_min
            # plot histogram
            hist, edges = np.histogram(elev, bins=nbins, density=True)
            plt.plot(edges[:-1], hist, label=column)
    # plot decorations
    plt.title(dem + ' - deviation from SRTM 30m')
    plt.xlabel('Difference (m)')
    plt.ylabel('Normalized probability density function')
    plt.legend()
    plt.tight_layout()
    plt.savefig(file_out)
    print 'histogram OK'
    # plt.show()
    plt.clf()
    plt.cla()


#----------------------------------------------------
# stats of differences from SRTM30 - per dem
for dem in ['tdx12','tdx30','aster30','aw3d30']:
    df_stats = pd.DataFrame()
    df_stats['stats'] = ['min','max','range','mean','median','stddev','skew','kurt','p25','p75']
    for column in df[headers]:
        if column.endswith(dem):
            elev = df[column].values
            df_stats[column] = calc_stats(elev)
    # export stats to latex
    print('\n\n')
    # print df_stats.set_index('stats').to_latex(float_format=lambda x:'%4.4f' % x)
    print df_stats.transpose().to_latex(float_format=lambda x:'%4.2f' % x)


#----------------------------------------------------
# calculate maps of differences
for area in areas:
    area_name = area[0].split('_')[-1]
    srtm = 'srtm30_wgs84_' + area_name
    for dem in area:
        diff_map = 'srtm_diff_' + dem.split('_')[0] + '_wgs84_' + area_name
        grass.run_command('g.region', raster=dem, flags='pa')
        grass.mapcalc("${out} = ${srtm} - ${dem}",
            out=diff_map,
            srtm=srtm,
            dem=dem,
            overwrite = True)

# fix  colors, flags for r.colors: e=equalize hist, a=absolute log scale
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        diff_map = 'srtm_diff_' + dem.split('_')[0] + '_wgs84_' + area_name
        grass.run_command('g.region', raster=dem, flags='pa')
        grass.run_command('r.colors.stddev', map=diff_map, flags='z')
        # grass.run_command('r.colors', map=diff_map, color='differences', flags='a')

# make a list of the diffs maps, to exclude srtm30-srtm30
diffs_list = []
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        diff_map = 'srtm_diff_' + dem.split('_')[0] + '_wgs84_' + area_name
        diffs_list.append(diff_map)

diffs_list = [dmap for dmap in diffs_list if not dmap.startswith('srtm_diff_srtm')]

# export diffs maps and histograms as pngs
for diff_map in diffs_list:
    grass.run_command('g.region', raster=diff_map, flags='pa', res='0:0:01')
    # diff map
    grass.run_command('d.mon', start='png', output=diff_map+'.png', resolution = '2', \
        height=500, width=500, overwrite=True)
    grass.run_command('d.rast', map=diff_map)
    grass.run_command('d.grid', size='0.25', text_color='black', flags='c')
    grass.run_command('d.legend', raster=diff_map, flags='tsbd', at='4,25,6,8', font='Helvetica', fontsize=12, bgcolor='240:240:240')
    grass.run_command('d.mon', stop='png')
    # histogram
    grass.run_command('d.mon', start='png', output=diff_map+'_hist.png', resolution = '2', \
        height=350, width=600, overwrite=True)
    grass.run_command('d.histogram', map=diff_map)
    grass.run_command('d.mon', stop='png')


#----------------------------------------------------
# scatter plots of original elevations (not difference)
# manually delete srtm x srtm
for area in areas:
    area_name = area[0].split('_')[-1]
    csv_file = area_name + '_randon_' + str(points) + '.csv'
    df = pd.read_csv(work_dir + csv_file, index_col=0)
    for dem in area:
        file_svg = dem + '_scatter.pdf'
        col_name = dem#.split('_')[0]
        y_min_srtm = min(df['srtm30_wgs84_'+area_name])
        y_min_col = min(df[col_name])
        y_min_plt = min(y_min_srtm,y_min_col)
        plt.plot(df['srtm30_wgs84_'+area_name],df[col_name],'o')
        plt.title(area_name + ' - SRTM 30m x ' + col_name)
        plt.xlabel('SRTM 30m elevation')
        plt.ylabel(col_name + ' elevation')
        xmin = plt.gca().get_xlim()[0]
        xmax = plt.gca().get_xlim()[1]
        plt.plot((xmin,xmax),(xmin,xmax))
        plt.ylim(xmin,xmax)
        plt.xlim(xmin,xmax)
        # plt.legend()
        plt.tight_layout()
        plt.savefig(file_svg)
        print 'scatterplot OK'
        # plt.show()
        plt.clf()
        plt.cla()


#----------------------------------------------------
# Timelapse and topo profiles - this is all pretty much manual...
#----------------------------------------------------
# set working directory
work_dir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/pngs_diff_inset/'
os.chdir(work_dir)

# Barcelos
area = 'Barcelos'
n='0:27s'
s='0:41s'
w='62:48w'
e='62:34w'
make_diff_map_png(tdx='tdx30',area=area,n=n,s=s,w=w,e=e,shade_map=False,grid_size='0.05',suffix='')
coords = '62:45:43w,0:35:47s,62:35:25w,0:30:30s'
fig_aspect = 200
make_diff_profile(area,n,s,w,e,coords,fig_aspect)

# Iporanga
area = 'Iporanga'
n='24:40s'
s='24:45s'
w='48:10:30w'
e='48:05:30w'
make_diff_map_png(tdx='tdx30',area=area,n=n,s=s,w=w,e=e,shade_map=False,grid_size='0.05',suffix='')
coords = '48:08:45w,24:42:40s,48:06:40w,24:41:20s'
fig_aspect = 10
make_diff_profile(area,n,s,w,e,coords,fig_aspect)

# SerraMar - gridded artifact from difference in resolution
# looks a lot like projection errors from lat/long to utm...
area = 'SerraMar'
n='22:19s'
s='22:21s'
w='44:40w'
e='44:38w'
make_diff_map_png(tdx='tdx30',area=area,n=n,s=s,w=w,e=e,shade_map=True,grid_size='0.01',suffix='_proj_error_30m')
make_diff_map_png(tdx='tdx12',area=area,n=n,s=s,w=w,e=e,shade_map=True,grid_size='0.01',suffix='_proj_error_12m')


#----------------------------------------------------
# Shades of Tucum river erosion area
#----------------------------------------------------
# set working directory
work_dir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/rand_10k_csv/'
os.chdir(work_dir)


grid_size='0.01'
dems = ['aster30','aw3d30_wgs84','srtm30','tdx30','tdx12']
for dem in dems:
    shade_map = dem+'_RioClaro_shade_315_20'
    grass.run_command('g.region', raster=shade_map, flags='pa')
    grass.run_command('g.region', region='tucum', flags='pa')
    grass.run_command('d.mon', start='png', output=dem+'_shade_tucum.png', resolution = '3', height=500, width=500, overwrite=True)
    grass.run_command('d.rast', map=shade_map)
    grass.run_command('d.grid', size=grid_size, text_color='black', flags='c')
    grass.run_command('d.mon', stop='png')
