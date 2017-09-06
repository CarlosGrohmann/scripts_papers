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


# This file has functions used to analyse raster and vector data
# these functions are imported by the analysis script


import sys, os, csv
import itertools
import numpy as np
import scipy.stats as ss
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import gc

try:
    import grass.script as grass
    import grass.script.array as garray
except:
    pass
    
#----------------------------------------------------
def rm_mask(kind='cell'):
    ''' removes any existing mask, defaults to raster '''
    exists = bool(grass.find_file('MASK', element=kind)['file'])
    if exists:
        grass.run_command('r.mask', flags='r') 


#----------------------------------------------------
def flat_param(parameter, aslist=False):
    ''' return GRASS raster as numpy array or flat list removing the nulls '''
    print ('----------   flatParam   ----------')
    param = garray.array()
    param.read(parameter, null=np.nan)
    # total_cells = param.size # total number of cells
    param_nan = param[~np.isnan(param)]
    # null_cells = total_cells - param_nan.size # number of null cells
    if aslist == True:
        parflat = param_nan.tolist()
    else:
        parflat = np.array(param_nan)
    param = None
    param_nan = None
    print ('----------  flatParam OK ----------')
    gc.collect()
    return parflat


#----------------------------------------------------
def round5(x):
    ''' round to nearest 5 ''' 
    rounded = int(round(x/5.0)*5.0)
    return rounded

#----------------------------------------------------
def round_arb(x,y):
    ''' round to nearest (arbitrary) value ''' 
    rounded = int(round(x/float(y))*float(y))
    return rounded


#----------------------------------------------------
def rmse(predictions, targets):
    ''' calculate RMSE from 2 values or lists or arrays ''' 
    return np.sqrt(((predictions - targets) ** 2).mean())


#----------------------------------------------------
def calc_stats(dem):
    ''' get stats from a list of values 
    (min, max, range, mean, median, stddev, p25, p75, skew, kurt) 
    removes NaNs from array'''
    print ('----------   calc_stats   ----------')
    dem = dem[~np.isnan(dem)]
    dem_max = np.max(dem)
    dem_min = np.min(dem)
    dem_range = dem_max - dem_min
    dem_mean = np.mean(dem)
    dem_median = np.median(dem)
    dem_stddev = np.std(dem)
    dem_p25 = np.percentile(dem, 25)
    dem_p75 = np.percentile(dem, 75)
    dem_skew = ss.skew(dem)
    dem_kurt = ss.kurtosis(dem)
    print ('---------- calc_stats OK  ----------')
    return [dem_min,dem_max,dem_range,dem_mean,dem_median,dem_stddev,dem_skew,dem_kurt,dem_p25,dem_p75]


#----------------------------------------------------
def sample_dems(dem, n_points, r_seed, ow_vector, ow_what, coords=False, flags=None):
    ''' create random points for selected tile and sample the elevation values 
    # returns list of stats or np.array of elevation values
    dem = raster map (elevation)
    n_points = number of random points
    r_seed = random seed
    ow_vector = should vector maps be overwritten?
    ow_what = re-run r.what.rast ? 
    '''
    print ('----------   sample_dems   ----------')
    grass.run_command('g.region', raster=dem, flags='pa')
    # coastal areas need to use a vector map to constrain random points only to terrain
    # this mask was created in the pre-processing step (import_data.py)
    poly_name = dem.split('_')[-1] + '_poly'
    find_mask = grass.find_file(poly_name, element = 'vector')
    if find_mask['name'] == "":
        vect_mask = None
    else:
        vect_mask = poly_name
    # create random points - check if map exists
    vector_name = dem + '_random_' + str(n_points) #dem.split('_')[-1] + '_random_' + str(n_points)
    find_vect = grass.find_file(vector_name, element = 'vector')
    # if there is no vector map OR if ow_vector is True, create vector map
    if find_vect['name'] == "" or ow_vector==True:
        grass.run_command('v.random', output=vector_name, npoints=n_points, seed=r_seed, column='z', restrict=vect_mask, overwrite=ow_vector)
    # add column with number of points as name, check if already exists
    column_existing = grass.vector_columns(vector_name).keys()
    rand_col = dem#.split('_')[0]
    if rand_col not in column_existing:
        grass.run_command('v.db.addcolumn', map=vector_name, columns=rand_col+' double precision', quiet=True)
    # sample raster map - force if overwrite vector is true
    if ow_vector or ow_what:
        grass.run_command('v.what.rast', map=vector_name, raster=dem, column=rand_col, quiet=True)
    # export as ascii and read into python
    xyz = grass.read_command('v.out.ascii', input=vector_name, type='point', format='point', columns=rand_col, overwrite=True)
    elev_list = [float(attr.split('|')[3]) if attr.split('|')[3] != '' else None for attr in xyz.split('\n')[:-1]]
    elev = np.asarray(elev_list, dtype=np.float64) 
    if coords==True:
        x_list = [float(attr.split('|')[0]) if attr.split('|')[0] != '' else None for attr in xyz.split('\n')[:-1]]
        y_list = [float(attr.split('|')[1]) if attr.split('|')[1] != '' else None for attr in xyz.split('\n')[:-1]]
    else:
        x_list = None
        y_list = None
    print ('---------- sample_dems OK  ----------')
    return x_list, y_list, elev


#----------------------------------------------------
def sample_slopes(dem, slope, n_points, r_seed, ow_vector, ow_what):
    ''' create random points for selected tile and sample the elevation values 
    # returns list of stats or np.array of elevation values
    dem = raster map (elevation)
    n_points = number of random points
    r_seed = random seed
    ow_vector = should vector maps be overwritten?
    ow_what = re-run r.what.rast ? 
    '''
    print ('----------   sample_dems   ----------')
    grass.run_command('g.region', raster=dem, flags='a')
    # coastal areas need to use a vector map to constrain random points only to terrain
    # this mask was created in the pre-processing step (import_data.py)
    poly_name = dem.split('_')[-1] + '_poly'
    find_mask = grass.find_file(poly_name, element = 'vector')
    if find_mask['name'] == "":
        vect_mask = None
    else:
        vect_mask = poly_name
    # create random points - check if map exists
    vector_name = dem + '_random_' + str(n_points) #dem.split('_')[-1] + '_random_' + str(n_points)
    find_vect = grass.find_file(vector_name, element = 'vector')
    # if there is no vector map OR if ow_vector is True, create vector map
    if find_vect['name'] == "" or ow_vector==True:
        grass.run_command('v.random', output=vector_name, npoints=n_points, seed=r_seed, column='z', restrict=vect_mask, overwrite=ow_vector)
    # add column with number of points as name, check if already exists
    column_existing = grass.vector_columns(vector_name).keys()
    rand_col = dem.split('_')[0]
    if rand_col not in column_existing:
        grass.run_command('v.db.addcolumn', map=vector_name, columns=rand_col+' double precision', quiet=True)
    # sample raster map - force if overwrite vector is true
    if ow_vector or ow_what:
        grass.run_command('v.what.rast', map=vector_name, raster=dem, column=rand_col, quiet=True)
    # export as ascii and read into python
    xyz = grass.read_command('v.out.ascii', input=vector_name, type='point', format='point', columns=rand_col, overwrite=True)
    elev_list = [float(attr.split('|')[3]) if attr.split('|')[3] != '' else None for attr in xyz.split('\n')[:-1]]
    elev = np.asarray(elev_list, dtype=np.float64) 
    print ('---------- sample_dems OK  ----------')
    return elev


#----------------------------------------------------
def sample_dems_mc(dem, n_points, mc_run, ow_vector, ow_what):
    ''' create random points for selected tile and sample the elevation values 
    simplified version of sample_dems tailored for MonteCarlo-like analysis
    dem = raster map (elevation)
    n_points = number of random points
    mc_run = number of times a set of random pints will be created 
    ow_vector = should vector maps be overwritten?
    ow_what = re-run r.what.rast ? 

    note: to keep random points really randon and yet ensure reproducibility,
    random seed is set to the value of mc_run * 42'''
    print ('----------   sample_dems_mc   ----------')
    grass.run_command('g.region', raster=dem, flags='a')
    # chack if vector map of random points alreay exists
    vector_name = dem.split('_')[1] + '_random_' + str(n_points) + '_' + str(mc_run).zfill(2)
    find_elev = grass.find_file(vector_name, element = 'vector')
    grass.run_command('v.random', output=vector_name, npoints=n_points, seed=mc_run*42, column='z', quiet=True, overwrite=ow_vector)   
    # add column with number of points as name, check if already exists
    column_existing = grass.vector_columns(vector_name).keys()
    rand_col = 'rand_' + str(n_points) + '_' + str(mc_run)
    if rand_col not in column_existing:
        grass.run_command('v.db.addcolumn', map=vector_name, columns=rand_col+' double precision', quiet=True)
    # sample raster map - force if overwrite vector is true
    if ow_vector or ow_what:
        grass.run_command('v.what.rast', map=vector_name, raster=dem, column=rand_col, quiet=True)
    # export as ascii and read into python
    xyz = grass.read_command('v.out.ascii', input=vector_name, type='point', format='point', columns=rand_col, overwrite=True)
    elev_list = [float(attr.split('|')[3]) if attr.split('|')[3] != '' else None for attr in xyz.split('\n')[:-1]]
    elev = np.asarray(elev_list, dtype=np.float64)
    print ('---------- sample_dems_mc OK  ----------')
    return elev

        
#----------------------------------------------------
# fits a 4PL curve to mean of correlation values
# plots and funcs from http://people.duke.edu/~ccc14/pcfb/analysis.html
def logistic4(x, A, B, C, D):
    ''' 4PL logistic equation ''' 
    return ((A-D)/(1.0+((x/C)**B))) + D

def residuals(p, y, x):
    ''' Deviations of data from fitted 4PL curve ''' 
    A,B,C,D = p
    err = y-logistic4(x, A, B, C, D)
    return err

def peval(x, p):
    ''' Evaluated value at x with current parameters ''' 
    A,B,C,D = p
    return logistic4(x, A, B, C, D)


#----------------------------------------------------
# calculates a histogram for each DEM tile and saves it to an svg file
# uses matplotlib histogram
def do_histogram_np(dem_list,work_dir,bin_width,x_lim,density_plot, points):
    ''' calculates histogram on DEM tiles '''
    os.chdir(work_dir) # output dir
    area_name = dem_list[0].split('_')[1]
    fileOutName = area_name + '_histogram.pdf' # output file
    csv_file = area_name + '_randon_' + str(points) + '.csv' # csv file with elev for all dems
    df = pd.read_csv(work_dir + csv_file, index_col=0)
    sns.set_palette(sns.color_palette("Set1", n_colors=len(df.columns)+2, desat=.5))
    # get data for each dem (only columns with elev values)
    for column in df.iloc[:,2:]:
        elev = df[column].values
        elev = elev[~np.isnan(elev)]
        # define bins
        g_max = int(np.ceil(np.max(elev)))
        g_min = int(np.min(elev))
        nbins = round5((g_max - g_min)/bin_width)
        # plot histogram
        # plt.hist(elev, bins=nbins, normed=density_plot, histtype='step', label=dem)
        hist, edges = np.histogram(elev, bins=nbins, density=density_plot)
        plt.plot(edges[:-1], hist, label=column)
    # plot decorations
    plt.title(area_name + ' (bin width = ' + str(bin_width) +'m)')
    plt.xlabel('Elevation')
    if density_plot == True:
        plt.ylabel('Normalized probability density function')
    else:
        plt.ylabel('Cell count')
    plt.xlim(x_lim)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fileOutName)
    print 'histogram OK'
    # plt.show()
    plt.clf()
    plt.cla()


#----------------------------------------------------
# calculates a histogram for each DEM tile and saves it to an svg file
# uses matplotlib histogram
# uses min-max as limits, one dem per run (does not accept lists of dems)
def do_histogram_full(dem,work_dir,bin_width, points):
    ''' calculates histogram on DEM tiles '''
    os.chdir(work_dir) # output dir
    area_name = dem.split('_')[1]
    fileOutName = area_name + '_histogram_full_srtm.pdf' # output file
    csv_file = area_name + '_randon_' + str(points) + '.csv' # csv file with elev for all dems
    df = pd.read_csv(work_dir + csv_file, index_col=0)
    # get data for each dem
    column = dem
    elev = df[column].values
    elev = elev[~np.isnan(elev)]
    # define bins
    g_max = int(np.ceil(np.max(elev)))
    g_min = int(np.min(elev))
    nbins = round5((g_max - g_min)/bin_width)
    hist, edges = np.histogram(elev, bins=nbins, density=True)
    plt.plot(edges[:-1], hist, label=column)
    # plot decorations
    plt.title(area_name + ' (bin width = ' + str(bin_width) +'m)')
    plt.xlabel('Elevation')
    plt.ylabel('Normalized probability density function')
    plt.xlim((g_min,g_max))
    plt.legend()
    plt.tight_layout()
    plt.savefig(fileOutName)
    print 'histogram OK'
    # plt.show()
    plt.clf()
    plt.cla()


#----------------------------------------------------
# calculates histograms for slope
def do_histogram_slope(dem_list,work_dir,bin_width,x_lim,density_plot,points,n_colors):
    ''' calculates histogram on DEM tiles '''
    os.chdir(work_dir) # input/output dir
    area_name = dem_list[0].split('_')[1]
    fileOutName = area_name + '_histogram_slope.pdf' # output file
    for dem in dem_list:
        slope_rast = dem.split('_')[0]+'_slope_'+area_name
        csv_file = dem + '_elev_w_slope_'+str(points)+'pts.csv'
        df = pd.read_csv(csv_file, index_col=0)
        slope = df[slope_rast].values
        sns.set_palette(sns.color_palette("Set1", n_colors=n_colors, desat=.5))
        slope = slope[~np.isnan(slope)]
        # define bins
        g_max = int(np.ceil(np.max(slope)))
        g_min = int(np.min(slope))
        nbins = round5((g_max - g_min)/bin_width)
        # plot histogram
        # plt.hist(slope, bins=nbins, normed=density_plot, histtype='step', label=dem)
        hist, edges = np.histogram(slope, bins=nbins, density=density_plot)
        plt.plot(edges[:-1], hist, label=dem)
    # plot decorations
    plt.title(area_name + ' slope (bin width = ' + str(bin_width) +'m)')
    plt.xlabel('Slope (degrees)')
    if density_plot == True:
        plt.ylabel('Normalized probability density function')
    else:
        plt.ylabel('Cell count')
    plt.xlim(x_lim)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fileOutName)
    print 'histogram OK'
    # plt.show()
    plt.clf()
    plt.cla()


#----------------------------------------------------
# plot mean_slope per elevation
def mean_slope_elev(area_list,elev_interval,zmin,zmax,smin,smax,points,work_dir,suffix):
    os.chdir(work_dir)
    area_name = area_list[0].split('_')[-1]
    file_svg = area_name + '_mean_slope_elev_'+suffix+'.pdf'
    sns.set_palette(sns.color_palette("muted", n_colors=5, desat=.5))
    for dem in area_list:
        df = pd.DataFrame()
        csv_file = dem + '_elev_w_slope_'+str(points)+'pts.csv' 
        slope = dem.split('_')[0]+'_slope_'+area_name
        df = pd.read_csv(csv_file)#, header=None, names=['z','s'], index_col=None,sep='\s+')
        # sort dem, define breaks
        df_sort = df.sort_values(by=dem)
        # elev_min = int(round(df[dem].min()/elev_interval)*elev_interval)
        # elev_max = int(round(df[dem].max()/elev_interval)*elev_interval)
        bins = range(zmin,zmax+1,int(elev_interval))
        # aggregate slope into breaks
        df_sort['cut'] = pd.cut(df_sort[dem], bins=bins)
        pvt = pd.pivot_table(df_sort, columns='cut', values=slope, aggfunc=np.mean)
        # plot 
        # x = pvt.values
        # y = bins[:-1]
        df2 = pd.DataFrame()
        df2['x'] = pvt.values
        df2['y'] = bins[:-1]
        df2 = df2[df2.x != np.nan]
        plt.plot(df2['x'],df2['y'], label=dem)
    # plot decorations
    plt.title(area_name + ' - mean slope x elevation ')
    plt.xlabel('Mean slope - interval:' + str(elev_interval) + ' m')
    plt.ylabel('Elevation')
    # zdelta = 0.1 * (zmax-zmin)
    # sdelta = 0.1 * (smax-smin)
    plt.ylim(zmin, zmax)
    plt.xlim(0,smax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(file_svg)
    print 'plot OK'
    # plt.show()
    plt.clf()
    plt.cla()


#----------------------------------------------------
# error metrics
# error metrics are just a copy from:
# http://www.statsmodels.org/dev/_modules/statsmodels/tools/eval_measures.html
# except mean error, that was adapted from mean squared error
def err_me(x1, x2, axis=0):
    """mean error"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    return np.mean((x1-x2), axis=axis)

def err_mse(x1, x2, axis=0):
    """mean squared error"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    return np.mean((x1-x2)**2, axis=axis)

def err_rmse(x1, x2, axis=0):
    """root mean squared error"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    return np.sqrt(err_mse(x1, x2, axis=axis))

def err_mae(x1, x2, axis=0):
    """mean absolute error"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    return np.mean(np.abs(x1-x2), axis=axis)

def err_std(x1, x2, axis=0):
    """standard deviation of error"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    return np.std(x1-x2, axis=axis)


#----------------------------------------------------
# calculate error metric (rmse, mae, me, std)
def calc_error_metrics(areas,work_dir,points):
    dems = []
    rmse_list = []
    mae_list = []
    me_list = []
    std_list = []
    n_nans = []
    area_list = []
    for area in areas:
        area_name = area[0].split('_')[-1]
        csv_file = area_name + '_randon_' + str(points) + '.csv'
        print(csv_file)
        # get elev data from csv files 
        df = pd.read_csv(work_dir + csv_file, index_col=0)
        # caclulate RMSE about tdx12
        for column in ['tdx12_'+area_name,'tdx30_'+area_name,'aster30_wgs84_'+area_name,'aw3d30_wgs84_'+area_name]:
            # set values of column srtm30 to be NaN if they are NaN at column in the loop
            df.loc[df['srtm30_wgs84_'+area_name].isnull(), column] = np.nan
            # set values of column in the loop to be NaN if they are NaN at column srtm30
            df.loc[df[column].isnull(), 'srtm30_wgs84_'+area_name] = np.nan
            # remove NaNs
            srtm30 = df['srtm30_wgs84_'+area_name].values
            srtm30 = srtm30[~np.isnan(srtm30)]
            dem2 = df[column].values
            dem2 = dem2[~np.isnan(dem2)]
            # calculate errors
            rmse_val = err_rmse(srtm30,dem2)
            mae_val = err_mae(srtm30,dem2)
            me_val = err_me(srtm30,dem2)
            std_val = err_std(srtm30,dem2)
            # append data to lists
            rmse_list.append(rmse_val)
            mae_list.append(mae_val)
            me_list.append(me_val)
            std_list.append(std_val)
            dems.append(column)
            n_nans.append(len(dem2))
            area_list.append(area_name)
    return rmse_list, mae_list, me_list, std_list, dems, n_nans, area_list


#----------------------------------------------------
#  make a png of a diff map and plot topo profile
def make_diff_map_png(tdx,area,n,s,w,e,shade_map,grid_size,suffix):
    ''' make a png of a diff map and plot topo profile
    tdx: 'tdx12' or 'tdx30'
    area: name of study area (as in 'Barcelos')
    n,s,w,e: limits of GRASS region
    shade_map: boolean
    grid_size: interval of grid on output map
    suffix: so we can have outputs name differently
    '''
    diff_map = 'srtm_diff_'+tdx+'_wgs84_'+area
    diff_map_inset = diff_map+'_inset'+suffix
    shade = tdx+'_'+area+'_shade_315_20'
    grass.run_command('g.region', raster=diff_map, flags='a')
    grass.run_command('g.region', n=n, s=s, w=w, e=e)
    # make new map (for new color table)
    print('clipping with r.mapcalc...')
    grass.mapcalc("${out} = ${orig}",
        out=diff_map_inset,
        orig=diff_map,
        overwrite = True)
    # export map as png
    if shade_map==True:
        grass.run_command('d.mon', start='png', output=diff_map_inset+'_shade.png', resolution = '2', height=500, width=500, overwrite=True)
        grass.run_command('d.shade', shade=shade, color=diff_map)
    else:
        grass.run_command('d.mon', start='png', output=diff_map_inset+'.png', resolution = '2', height=500, width=500, overwrite=True)
        grass.run_command('d.rast', map=diff_map_inset)
    grass.run_command('d.grid', size=grid_size, text_color='black', flags='c')
    grass.run_command('d.legend', raster=diff_map_inset, flags='tsbd', at='4,25,92,94', font='Helvetica', fontsize=12, bgcolor='240:240:240')
    grass.run_command('d.mon', stop='png')
    print('export as png OK')
    # export simple PNG + world file
    grass.run_command('r.out.png',input=diff_map_inset, output=diff_map_inset+'_wld.png', flags='w', overwrite=True)


def make_diff_profile(area,n,s,w,e,coords,fig_aspect):
    df_tdx = None
    df_srtm = None
    grass.run_command('g.region', n=n, s=s, w=w, e=e, res='0:0:01')
    # make profile
    csv_file_tdx = 'tdx30_'+area+'_profile.csv'
    csv_file_srtm = 'srtm30_'+area+'_profile.csv'
    grass.run_command('r.profile', input='tdx30_'+area, output=csv_file_tdx, coordinates=coords, overwrite=True)
    grass.run_command('r.profile', input='srtm30_wgs84_'+area, output=csv_file_srtm, coordinates=coords, overwrite=True)
    # read csv and plot profile
    df_tdx = pd.read_csv(work_dir + csv_file_tdx, sep=' ', header=None, names=['nan','dist','elev'],na_values='*')
    df_srtm = pd.read_csv(work_dir + csv_file_srtm, sep=' ', header=None, names=['nan','dist','elev'],na_values='*')
    file_svg = area + '_profile.pdf'
    plt.plot(df_tdx['dist'],df_tdx['elev'],'k-',label='TanDEM-X 30M')
    plt.plot(df_srtm['dist'],df_srtm['elev'],'b-',label='SRTM 30M')
    xmax = round_arb(df_tdx['dist'].max(),500)
    plt.xlim(0,xmax)
    plt.axes().set_aspect(fig_aspect)
    plt.title('Topographic profiles - ' + area)
    plt.xlabel('Distance along profile')
    plt.ylabel('Elevation')
    plt.legend()
    # plt.show()
    # plt.tight_layout()
    plt.savefig(file_svg, format='pdf', bbox_inches = 'tight')
    plt.clf()   
    plt.cla()
    print('export plot as PDF OK')



