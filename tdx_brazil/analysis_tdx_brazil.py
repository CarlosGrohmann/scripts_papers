
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
# funcDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/_submission_2/scripts_rev1/'
os.chdir(funcDir)
from functions import *


#----------------------------------------------------
# WGS84 x EGM96
# stats of DEMs before and after converting elevation to WGS84 ellipsoid
# conversion to WGS84 is in the import_tdx_brazil.py file 
#----------------------------------------------------

#----------------------------------------------------
# set working directory
work_dir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/stats_plots/'
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
# read data from GRASS rasters and calculate stats - display as latex tables
for area in areas:
    df_stats = pd.DataFrame()
    df_stats['stats'] = ['min','max','range','mean','median','stddev','skew','kurt','p25','p75']    
    for dem in area:
        dem_array = raster_as_1d_array(dem)
        df_stats[dem] = calc_stats(dem_array)
    # export stats to latex
    print('\n\n')
    print df_stats.transpose().to_latex(float_format=lambda x:'%4.2f' % x)


#----------------------------------------------------
#----------------------------------------------------
# WGS84
# actual analysis starts here. all elevations are WGS84
#----------------------------------------------------
#----------------------------------------------------

# set working directory
work_dir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/stats_plots/'
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
# calculate histograms of elevation, one area at a time, so we can set up limits and bins individually 
os.chdir(work_dir) # output dir
do_histogram_np(araca, bin_width=2, x_lim=(-10,150), density_plot=True)
do_histogram_np(barcelos, bin_width=2, x_lim=(-10,120), density_plot=True)
do_histogram_np(pantanal, bin_width=2, x_lim=(50,200), density_plot=True)
do_histogram_np(serramar, bin_width=20, x_lim=(-10,2700), density_plot=True)
do_histogram_np(rioclaro, bin_width=10, x_lim=(400,1100), density_plot=True)
do_histogram_np(iporanga, bin_width=20, x_lim=(-10,1250), density_plot=True)
do_histogram_np(santacatarina, bin_width=2, x_lim=(-2,150), density_plot=True)

# full histograms (limits = min-max) for srtm, for each area (they will be insets on the hists for all dems)
os.chdir(work_dir) # output dir
do_histogram_full('srtm30_wgs84_Araca', bin_width=10)
do_histogram_full('srtm30_wgs84_Barcelos', bin_width=1)
do_histogram_full('srtm30_wgs84_Pantanal', bin_width=1)
do_histogram_full('srtm30_wgs84_SerraMar', bin_width=10)
do_histogram_full('srtm30_wgs84_RioClaro', bin_width=5)
do_histogram_full('srtm30_wgs84_Iporanga', bin_width=10)
do_histogram_full('srtm30_wgs84_SantaCatarina', bin_width=10)


#----------------------------------------------------
# slope
# calculate slope in GRASS
for area in areas:
    for dem in area:
        slope = dem + '_slope'
        grass.run_command('g.region', raster=dem, flags='pa')
        grass.run_command('r.slope.aspect', elevation=dem, slope=slope, overwrite=True)

# export slope maps as pngs
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        slope = dem + '_slope'
        grass.run_command('g.region', raster=slope, flags='pa', res='0:0:01')
        grass.run_command('d.mon', start='cairo', output=slope+'.png', resolution='3', \
            height=500, width=500, overwrite=True)
        grass.run_command('d.rast', map=slope)
        grass.run_command('d.grid', size='0.25', text_color='black', flags='c')
        grass.run_command('d.legend', raster=slope, flags='sb', at='4,25,2,4', \
            font='Helvetica', fontsize=14, bgcolor='240:240:240', range=(0,90))
        grass.run_command('d.mon', stop='cairo')

# stats of slope
# reads data from GRASS rasters and calculates stats - display as latex tables
for area in areas:
    df_stats = pd.DataFrame()
    df_stats['stats'] = ['min','max','range','mean','median','stddev','skew','kurt','p25','p75']    
    for dem in area:
        slope = dem + '_slope'
        slope_array = raster_as_1d_array(slope)
        df_stats[slope] = calc_stats(slope_array)
    # export stats to latex
    print('\n\n')
    print df_stats.transpose().to_latex(float_format=lambda x:'%4.2f' % x)

# calculate histograms of slope 
os.chdir(work_dir) # output dir
for area_list in areas:
    do_histogram_slope(area_list, bin_width=1, x_lim=(0,90), density_plot=True, n_colors=len(area_list))

# Mean slope versus elevation
os.chdir(work_dir) # output dir
mean_slope_elev(area_list=araca         ,elev_interval=50 , ymin=0  , ymax=1800)
mean_slope_elev(area_list=barcelos      ,elev_interval=2  , ymin=0  , ymax=180)
mean_slope_elev(area_list=pantanal      ,elev_interval=5  , ymin=0  , ymax=200)
mean_slope_elev(area_list=iporanga      ,elev_interval=50 , ymin=0  , ymax=1350)
mean_slope_elev(area_list=rioclaro      ,elev_interval=50 , ymin=300, ymax=1100)
mean_slope_elev(area_list=santacatarina ,elev_interval=50 , ymin=0  , ymax=1350)
mean_slope_elev(area_list=serramar      ,elev_interval=100, ymin=0  , ymax=2800)


#----------------------------------------------------
# comparison of contour lines
# dict for contour levels
cont_dict = {'Araca':(20,0,1750), 'Barcelos':(2,0,150), 'Pantanal':(2,20,200), \
        'SerraMar':(20,0,2850), 'RioClaro':(10,370,1100), 'Iporanga':(20,0,1300), 'SantaCatarina':(5,0,1300),}

# calculate contours
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        grass.run_command('g.region', raster=dem, flags='pa')
        vect_contour = dem + '_contours'
        interval = cont_dict[area_name][0]
        min_c = cont_dict[area_name][1]
        max_c = cont_dict[area_name][2]
        grass.run_command('r.contour', input=dem, output=vect_contour, step=interval, minlevel=min_c, maxlevel=max_c, overwrite=True)

# calculate length of contours, contours per level and export to CSV
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        vect_contour = dem + '_contours'
        grass.run_command('v.db.addcolumn', map=vect_contour, columns="length double, count int", overwrite=True)
        grass.run_command('v.to.db', map=vect_contour, type='line', option='length', columns='length', overwrite=True)
        grass.run_command('v.to.db', map=vect_contour, type='line', option='count', columns='count', overwrite=True)
        grass.run_command('v.db.select', map=vect_contour, \
            columns='level,length,count', separator='comma', file=vect_contour+'.csv', overwrite=True)

# read CSV, add to pandas dataframe and plot stats and tables
# cumulative sum of length per interval
for area in areas:
    area_name = area[0].split('_')[-1]
    file_out = area_name + '_cumsum_contour_length.pdf'
    sns.set_palette(sns.color_palette("muted", n_colors=len(area), desat=.5))
    for dem in area:
        vect_contour = dem + '_contours'
        df_vect = pd.read_csv(work_dir+vect_contour+'.csv')
        plt.plot(df_vect['level'], df_vect['length'].cumsum(), label=dem)
    # plot decorations
    plt.title(area_name + ' cumsum - lenght')
    plt.xlabel('Elevation (m)')
    plt.ylabel('Cumulative Sum of Contour Lenght (m)')
    # plt.xlim((g_min,g_max))
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(file_out)
    print 'CumSum Plot OK'
    # plt.show()
    plt.clf()
    plt.cla()

# cumulative sum of number of countour curves per interval
for area in areas:
    area_name = area[0].split('_')[-1]
    file_out = area_name + '_cumsum_contour_count.pdf'
    sns.set_palette(sns.color_palette("muted", n_colors=len(area), desat=.5))
    for dem in area:
        vect_contour = dem + '_contours'
        df_vect = pd.read_csv(work_dir+vect_contour+'.csv')
        plt.plot(df_vect['level'], df_vect['count'].cumsum(), label=dem)
    # plot decorations
    plt.title(area_name + ' cumsum - contour count')
    plt.xlabel('Elevation (m)')
    plt.ylabel('Cumulative Sum of Contour Count (m)')
    # plt.xlim((g_min,g_max))
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(file_out)
    print 'CumSum Plot OK'
    # plt.show()
    plt.clf()
    plt.cla()

# sum length and count per vector per area, plot summary
dict_contour = {}
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        dem_type = dem.split('_')[0]
        vect_contour = dem + '_contours'
        df_vect = pd.read_csv(work_dir+vect_contour+'.csv')
        dict_contour[dem] = [area_name,dem_type,df_vect['length'].sum(),df_vect['count'].sum()]

# add to DataFrame
df = pd.DataFrame.from_dict(dict_contour,orient='index')
df.columns = ['area','type','length','count']
# export stats to latex
print('\n\n')
print df.sort_values(by=['area','length']).to_latex(float_format=lambda x:'%4.2f' % x)

# plot total length versus area, for each dem
sns.stripplot(x='area', y='length', hue='type', data=df.sort_values(by=['area']), \
    palette={'tdx12':'g','tdx30':'m','srtm30':'r','aster30':'b','aw3d30':'c'}, size=8)
# plot decorations
plt.title('Total contour length')
plt.xlabel('Study area')
plt.ylabel('Total contours length (m x10^9)')
plt.ylim(0, 1.2e9)
plt.legend()
plt.tight_layout()
plt.savefig('total_contours_length.svg')
print 'plot OK'
plt.clf()
plt.cla()


# plot total count versus area, for each dem
sns.stripplot(x='area', y='count', hue='type', data=df.sort_values(by=['area']), \
    palette={'tdx12':'g','tdx30':'m','srtm30':'r','aster30':'b','aw3d30':'c'}, size=8)
# plot decorations
plt.title('Total contour count')
plt.xlabel('Study area')
plt.ylabel('Total contours count')
# plt.ylim(0, 1.2e9)
plt.legend()
plt.tight_layout()
plt.savefig('total_contours_count.svg')
print 'plot OK'
plt.clf()
plt.cla()


#----------------------------------------------------
# Deviations from TDX12 (rmse, mae, me, std)
# 1 - interpolate all rasters to 12m using bicubic interpolation
for area in areas:
    grass.run_command('g.region', raster=area[0], flags='pa')
    for dem in area:
        if not dem.startswith('tdx12'):
            dem_resamp = dem + '_12m_bicubic'
            print('resampling ' + dem + ' to ' + dem_resamp + '...')
            grass.run_command('r.resamp.interp', input=dem, output=dem_resamp, method='bicubic', overwrite=True)

# 2 - calculate error metrics
rmse_list, me_list, std_list, le90_list, le95_list, le99_list, dem_list, area_list, min_list, max_list = calc_error_metrics(areas)

# make list of lists and put all in a DataFrame
err4pd = [('DEM', dem_list),
     ('AREA', area_list),
     ('MIN', min_list),
     ('MAX', max_list),
     ('ME', me_list),                 
     ('RMSE', rmse_list),
     ('STDE', std_list),
     ('LE90', le90_list),
     ('LE95', le95_list),
     ('LE99', le99_list),]
pd_err = pd.DataFrame.from_items(err4pd)

# export all metrics to latex
print pd_err.sort_values(by='DEM').to_latex(float_format=lambda x:'%4.2f' % x)

# export stats of each metric to latex
# srtm
print pd_err[pd_err['DEM'].str.startswith('srtm')].describe().to_latex(float_format=lambda x:'%4.3f' % x)
# aster
print pd_err[pd_err['DEM'].str.startswith('aster')].describe().to_latex(float_format=lambda x:'%4.3f' % x)
# alos
print pd_err[pd_err['DEM'].str.startswith('aw3d')].describe().to_latex(float_format=lambda x:'%4.3f' % x)
# tdx30
print pd_err[pd_err['DEM'].str.startswith('tdx30')].describe().to_latex(float_format=lambda x:'%4.3f' % x)


#----------------------------------------------------
# DEMs of differences - TDX12 minus all other DEMs (resampled to 12m)
# calculate maps of differences
for area in areas:
    area_name = area[0].split('_')[-1]
    grass.run_command('g.region', raster=area[0], flags='pa')
    tdx12 = 'tdx12_' + area_name
    for dem in area:
        dem_resamp = dem + '_12m_bicubic'
        if not dem.startswith('tdx12'):
            print('calculating DoD of TDX12 with ' + dem_resamp)
            diff_map = 'tdx12_diff_' + dem.split('_')[0] + '_wgs84_' + area_name
            grass.mapcalc("${out} = ${tdx} - ${dem}",
                out=diff_map,
                tdx=tdx12,
                dem=dem_resamp,
                overwrite = True)

# stats of DoDs
# reads data from GRASS rasters and calculates stats - display as latex tables
for area in areas:
    area_name = area[0].split('_')[-1]
    df_stats = pd.DataFrame()
    df_stats['stats'] = ['min','max','range','mean','median','stddev','skew','kurt','p25','p75']    
    for dem in area:
        if not dem.startswith('tdx12'):
            diff_map = 'tdx12_diff_' + dem.split('_')[0] + '_wgs84_' + area_name
            diff_map_array = raster_as_1d_array(diff_map)
            df_stats[diff_map] = calc_stats(diff_map_array)
    # export stats to latex
    print('\n\n')
    print df_stats.transpose().to_latex(float_format=lambda x:'%4.2f' % x)

# lists of all diff maps
diff_list = grass.list_grouped('raster', pattern='*diff*')['tandemX_brasil']

# fix colortable (color gradient from blue-white-red beteween -15/+15, then darker colors to min-max)
rule = '''0% 0:30:110 \n-15 blue \n0 white \n15 red \n100% 125:25:15'''
for diff_map in diff_list:
    grass.write_command('r.colors', map=diff_map, rules='-', stdin=rule)
    
# export diffs maps as pngs
for diff_map in diff_list:
    area_name = diff_map.split('_')[-1]
    grass.run_command('g.region', raster=diff_map, flags='pa', res='0:0:01')
    # diff map
    grass.run_command('d.mon', start='cairo', output=diff_map+'.png', resolution='3', \
        height=500, width=500, overwrite=True)
    grass.run_command('d.rast', map=diff_map)
    grass.run_command('d.grid', size='0.25', text_color='black', flags='c')
    grass.run_command('d.legend', raster=diff_map, flags='tsbd', at='4,25,6,8', \
        font='Helvetica', fontsize=12, bgcolor='240:240:240', range=(-30,30))
    grass.run_command('d.mon', stop='cairo')

# export histograms as pngs
for diff_map in diff_list:
    grass.run_command('g.region', raster=diff_map, flags='pa', res='0:0:01')
    # histogram
    grass.run_command('d.mon', start='cairo', output=diff_map+'_hist.png', resolution='3', \
        height=350, width=600, overwrite=True)
    grass.run_command('d.histogram', map=diff_map)
    grass.run_command('d.mon', stop='cairo')


#----------------------------------------------------
# scatter plots of aspect x differences
# calculate aspect in GRASS for resampled DEMs
for area in areas:
    for dem in area:
        if not dem.startswith('tdx12'):
            dem12 = dem + '_12m_bicubic'
            asp = dem12 + '_aspect'
            grass.run_command('g.region', raster=dem12, flags='pa')
            grass.run_command('r.slope.aspect', elevation=dem12, aspect=asp, overwrite=True)

# convert GRASS default aspect cartesian angles to azimuth angles
# see: Grohmann 2004. Caomputers and Geosciences 30:1055-1067
for area in areas:
    for dem in area:
        if not dem.startswith('tdx12'):
            dem12 = dem + '_12m_bicubic'
            asp = dem12 + '_aspect'
            asp_comp = asp + '_compass'
            grass.run_command('g.region', raster=dem12, flags='pa')
            grass.mapcalc("${out} = if(isnull(${asp_map}), null(), if((${asp_map} < 90), 90-${asp_map}, 360+90-${asp_map}))",
                out=asp_comp,
                asp_map=asp,
                overwrite = True)

# scatter plots of aspect x differences
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area:
        if not dem.startswith('tdx12'):
            dem_name = dem.split('_')[0]
            dem_diff = 'tdx12_diff_' + dem_name + '_wgs84_' + area_name
            if dem.startswith('tdx30'):
                asp = dem_name + '_' + area_name + '_12m_bicubic_aspect_compass'
            else:
                asp = dem_name + '_wgs84_' + area_name + '_12m_bicubic_aspect_compass'
            file_out = area_name + '_' + asp + '_scatter_dod.png'
            do_scatter_aspect(asp_name=asp,dod_name=dem_diff,file_out=file_out)

# local scatterplots - Serra do Mar
n='22:38s'
s='22:58s'
w='44:50w'
e='44:30w'
grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:00:00.4', flags='pa')

# tdx30
dem_diff = 'tdx12_diff_tdx30_wgs84_SerraMar'
asp = 'tdx30_SerraMar_12m_bicubic_aspect_compass'
file_out = 'tdx30_scatter_dod_aspect_Bocaina.png'
do_scatter_aspect(asp_name=asp, dod_name=dem_diff, file_out=file_out)

# srtm
dem_diff = 'tdx12_diff_srtm30_wgs84_SerraMar'
asp = 'srtm30_wgs84_SerraMar_12m_bicubic_aspect_compass'
file_out = 'srtm30_scatter_dod_aspect_Bocaina.png'
do_scatter_aspect(asp_name=asp, dod_name=dem_diff, file_out=file_out)

# used these params in do_scatter_aspect:
# plt.xlim(-20,380)
# plt.ylim(-50,50)
# plt.savefig(file_out, dpi=(600))

# make png of DoD for TDX30 in this sub-area
area_name = 'SerraMar'
dem_diff = 'tdx12_diff_tdx30_wgs84_SerraMar'
shade = 'tdx30_SerraMar_shade_315_20'
grass.run_command('d.mon', start='cairo', output=dem_diff+'_bocaina.png', resolution='3', \
    height=600, width=600, overwrite=True)
grass.run_command('d.shade', shade=shade, color=dem_diff, brighten=70)
grass.run_command('d.grid', size='0.1', text_color='black', flags='c', fontsize=18)
grass.run_command('d.legend', raster=dem_diff, flags='tsbd', at='4,25,6,8', \
    font='Helvetica', fontsize=14, bgcolor='240:240:240', range=(-50,50))
grass.run_command('d.mon', stop='cairo')

#----------------------------------------------------
#profile to check moving-window issues 
n='22:38s'
s='22:58s'
w='44:50w'
e='44:30w'

tdx12 = 'tdx12_SerraMar'
tdx30 = 'tdx30_SerraMar'
srtm = 'srtm30_wgs84_SerraMar'

# set size of moving-window here so we can change it later for all DEMs
rolling_size = 15


# run moving-window in raster
# grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:00:00.4', flags='pa')
# grass.run_command('r.neighbors', input=tdx12, output=tdx12+'_mean_15x15', size=15, flags='a', overwrite=True)
# (-44.583511547,-22.7918946833,-44.583511547,-22.8379865417) north-south

# get profile data, put it in a dataframe and create rolling window
# tdx12
grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:00:00.4', flags='pa')
p_tdx12 = grass.read_command('r.profile', input=tdx12 , output='-', coordinates=(-44.5991428842,-22.8482315113,-44.5852093365,-22.863772776), null='*') 
x_tdx12 = [float(i) for i in p_tdx12.split()[0::2]]
y_tdx12 = [float(i) for i in p_tdx12.split()[1::2]]
df = pd.DataFrame(list(zip(x_tdx12, y_tdx12)), columns=['X','Y'])
df_tdx12 = df.dropna()
# r_tdx12 = df_tdx12.rolling(rolling_size, center=False)

# tdx12 mean
# p_tdx12 = grass.read_command('r.profile', input=tdx12+'_mean_15x15' , output='-', coordinates=(-44.583511547,-22.7918946833,-44.583511547,-22.8379865417), null='*') 
# x_tdx12 = [float(i) for i in p_tdx12.split()[0::2]]
# y_tdx12 = [float(i) for i in p_tdx12.split()[1::2]]
# df = pd.DataFrame(list(zip(x_tdx12, y_tdx12)), columns=['X','Y'])
# df_tdx12_mw15 = df.dropna()

# tdx12 - inverted
# grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:00:00.4', flags='pa')
# p_tdx12 = grass.read_command('r.profile', input=tdx12 , output='-', coordinates=(-44.583511547,-22.7918946833,-44.583511547,-22.8379865417), null='*') 
# x_tdx12 = [float(i) for i in p_tdx12.split()[0::2]]
# y_tdx12 = [float(i) for i in p_tdx12.split()[1::2]]
# dfi = pd.DataFrame(list(zip(x_tdx12[::-1], y_tdx12[::-1])), columns=['X','Y'])
# df_tdx12i = dfi.dropna()
# r_tdx12i = df_tdx12i.rolling(rolling_size, center=False)

# tdx30
grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:00:01', flags='pa')
p_tdx30 = grass.read_command('r.profile', input=tdx30 , output='-', coordinates=(-44.5991428842,-22.8482315113,-44.5852093365,-22.863772776), null='*') 
x_tdx30 = [float(i) for i in p_tdx30.split()[0::2]]
y_tdx30 = [float(i) for i in p_tdx30.split()[1::2]]
df = pd.DataFrame(list(zip(x_tdx30, y_tdx30)), columns=['X','Y'])
df_tdx30 = df.dropna()

# srtm
grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:00:01', flags='pa')
p_srtm = grass.read_command('r.profile', input=srtm , output='-', coordinates=(-44.5991428842,-22.8482315113,-44.5852093365,-22.863772776), null='*') 
x_srtm = [float(i) for i in p_srtm.split()[0::2]]
y_srtm = [float(i) for i in p_srtm.split()[1::2]]
df = pd.DataFrame(list(zip(x_srtm, y_srtm)), columns=['X','Y'])
df_srtm = df.dropna()


# save plot
fig = plt.figure(figsize=(10,4))
# fig = plt.figure(figsize=(4,4)) #zoom
plt.plot(df_tdx12['X'], df_tdx12['Y'], label='TanDEM-X 12m')
# plt.plot(df_tdx12['X'], r_tdx12.mean()['Y'], label='TDX12m - moving-window: ' + str(rolling_size))
# plt.plot(df_tdx12i['X'], r_tdx12i.mean()['Y'], label='TDX12m - inverted moving-window: ' + str(rolling_size))
plt.plot(df_tdx30['X'], df_tdx30['Y'], label='TanDEM-X 30m')
plt.plot(df_srtm['X'], df_srtm['Y'], label='SRTM30m')
# plt.xlim(0,5150)
# plt.ylim(1025,1250)
# zoom
# plt.xlim(3000,3700)
# plt.ylim(1075,1225)
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Elevation (m)')
fig = 'profile_tdx12_tdx30_srtm_NW-SE.svg'
# fig = 'profile_tdx12_tdx30_mw5_zoom.svg'
plt.savefig(fig, dpi=(600))
plt.cla()
plt.clf()
plt.close(fig)



#----------------------------------------------------
# scatter plots of original elevations (not difference)
# remember to reset area list to the one with original (wgs84) dems, not the one with diffs

for area in areas:
    area_name = area[0].split('_')[-1]
    tdx = 'tdx12_' + area_name
    for dem in area:
        if not dem.startswith('tdx12'):
            dem_name = dem.split('_')[0]
            dem = dem + '_12m_bicubic'
            file_out = area_name + '_scatter_tdx12_' + dem_name + '.png'
            do_scatter_elev(tdx,dem,file_out,area_name,dem_name)

# ----------------------------------------------------
#----------------------------------------------------
# Export PNGs of 'zoom' areas
#----------------------------------------------------

# export PNGs pf contours with shaded relief
# set working directory
workDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/stats_plots/'
os.chdir(workDir)

# iporanga (Cajati Mining ara)
area_name = 'Iporanga'
cont_param = (20,0,300)
n='24:40s'
s='24:45s'
w='48:10:30w'
e='48:05:30w'
# make_zoom_contour(area_name,cont_param,n,s,w,e)
# make_shade_contour(n,s,w,e,tc='black',gc='black')
make_shadecolor_zoom(area_name,n,s,w,e,tc='white',gc='white')

# Santa catarina (coastal area + dune field)
area_name = 'SantaCatarina'
cont_param = (5,0,350)
n='27:59s'
s='28:01s'
w='48:39w'
e='48:37w'
# make_zoom_contour(area_name,cont_param,n,s,w,e)
# make_shade_contour(n,s,w,e,tc='black',gc='black')
make_shadecolor_zoom(area_name,n,s,w,e,tc='white',gc='white')

# barcelos
# area_name = 'Barcelos'
# cont_param = (5,0,60)
# n='0:27s'
# s='0:41s'
# w='62:48w'
# e='62:34w'
# # make_zoom_contour(area_name,cont_param,n,s,w,e)
# make_shade_contour(n,s,w,e,tc='black',gc='black')

# tucum
area_name = 'RioClaro'
cont_param = (5,450,650)
n='22:32:36S'
s='22:34:29S'
w='47:54:19W'
e='47:52:20W'
# make_zoom_contour(area_name,cont_param,n,s,w,e)
# make_shade_contour(n,s,w,e,tc='black',gc='black')
make_shadecolor_zoom(area_name,n,s,w,e,tc='white',gc='white')

# pantanal
area_name = 'Pantanal'
cont_param = (2,80,180)
n='18:55:50S'
s='18:58:50S'
w='56:27:00W'
e='56:24:00W'
# make_zoom_contour(area_name,cont_param,n,s,w,e)
# make_shade_contour(n,s,w,e,tc='white',gc='white')
make_shadecolor_zoom(area_name,n,s,w,e,tc='white',gc='white')

# Serra do Mar - pousada Sankay
# area_name = 'SerraMar'
# cont_param = (5,0,700)
# n='23:05:40S'
# s='23:07:40s'
# w='44:16w'
# e='44:14w'
# make_zoom_contour(area_name,cont_param,n,s,w,e)
# make_shade_contour(n,s,w,e,tc='white',gc='white')

# Serra do Mar - Alluvial fans
# area_name = 'SerraMar'
# cont_param = (5,500,1100)
# n='22:28S'
# s='22:32s'
# w='44:55w'
# e='44:51w'
# make_zoom_contour(area_name,cont_param,n,s,w,e)
# make_shade_contour(n,s,w,e,tc='white',gc='white')

# Serra do Mar - Void Filling?
area_name = 'SerraMar'
cont_param = (10,0,1800)
n='23:13:15S'
s='23:16:00s'
w='44:54:45w'
e='44:52:00w'
make_zoom_contour(area_name,cont_param,n,s,w,e)
make_shade_contour(n,s,w,e,tc='white',gc='white')
make_shadecolor_zoom(area_name,n,s,w,e,tc='white',gc='white')


#----------------------------------------------------
# Timelapse and topo profiles - this is all pretty much manual...
# Barcelos
area = 'Barcelos'
n='0:27s'
s='0:41s'
w='62:48w'
e='62:34w'
make_diff_map_png(dem='srtm30',area=area,n=n,s=s,w=w,e=e,shade_map=False,grid_size='0.05',suffix='')
coords = '62:45:43w,0:35:47s,62:35:25w,0:30:30s'
fig_aspect = 200
make_diff_profile(area,n,s,w,e,coords,fig_aspect)

# Iporanga
area = 'Iporanga'
n='24:40s'
s='24:45s'
w='48:10:30w'
e='48:05:30w'
make_diff_map_png(dem='srtm30',area=area,n=n,s=s,w=w,e=e,shade_map=False,grid_size='0.05',suffix='')
coords = '48:08:45w,24:42:40s,48:06:40w,24:41:20s'
fig_aspect = 10
make_diff_profile(area,n,s,w,e,coords,fig_aspect)


#----------------------------------------------------
# Shades of Tucum river erosion area (for fig. 2)
dems = ['aster30_wgs84','aw3d30_wgs84','srtm30_wgs84','tdx30','tdx12']
for dem in dems:
    shade_map = dem + '_RioClaro_shade_315_20'
    # grass.run_command('g.region', raster=shade_map, flags='pa')
    grass.run_command('g.region', region='tucum', flags='pa')
    grass.run_command('d.mon', start='cairo', output=dem+'_shade_tucum.png', resolution='3', height=500, width=500, overwrite=True)
    grass.run_command('d.rast', map=shade_map)
    grass.run_command('d.grid', size='0.025', text_color='black', fontsize=40, flags='c')
    grass.run_command('d.mon', stop='cairo')


#----------------------------------------------------
# Details of slope maps
# Barcelos
area_name = 'Barcelos'
n='0:16:24s'
s='0:20:24s'
w='62:09w'
e='62:05w'
make_shade_slope_zoom(area_name,n,s,w,e,tc='white',gc='white')

# pantanal
area_name = 'Pantanal'
n='18:55:50S'
s='18:58:50S'
w='56:27:00W'
e='56:24:00W'
make_shade_slope_zoom(area_name,n,s,w,e,tc='white',gc='white')


#----------------------------------------------------
# Water mask issues
# Araca
area_name = 'Barcelos'
make_shade_wam(area_name,tc='black',gc='black')

area_name = 'Itatiaia'
make_shade_wam(area_name,tc='black',gc='black')

area_name = 'Paraty'
make_shade_wam(area_name,tc='black',gc='black')
