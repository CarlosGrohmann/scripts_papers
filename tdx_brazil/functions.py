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
def raster_as_1d_array(raster):
    ''' return GRASS raster as numpy array - keep null values '''
    print ('----------   raster_as_1d_array   ----------')
    print (raster)
    grass.run_command('g.region', raster=raster, flags='pa')
    raster_array = garray.array()
    raster_array.read(raster, null=np.nan)
    print ('----------  raster_as_1d_array OK ----------')
    return raster_array.flatten(order='C')


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


def calc_stats(raster):
    ''' get stats from a list of values 
    (min, max, range, mean, median, stddev, p25, p75, skew, kurt) 
    removes NaNs from array'''
    print ('----------   calc_stats   ----------')
    rast = raster[~np.isnan(raster)]
    rast_max = np.max(rast).item()
    rast_min = np.min(rast).item()
    rast_range = rast_max - rast_min
    rast_mean = np.mean(rast).item()
    rast_median = np.median(rast).item()
    rast_stddev = np.std(rast).item()
    rast_p25 = np.percentile(rast, 25)
    rast_p75 = np.percentile(rast, 75)
    rast_skew = ss.skew(rast)
    rast_kurt = ss.kurtosis(rast)
    print ('---------- calc_stats OK  ----------')
    return [rast_min,rast_max,rast_range,rast_mean,rast_median,rast_stddev,rast_skew,rast_kurt,rast_p25,rast_p75]


#----------------------------------------------------
# calculates a histogram for each DEM tile and saves it to an svg file
# uses matplotlib histogram
def do_histogram_np(dem_list,bin_width,x_lim,density_plot):
    ''' calculates histogram on DEM tiles '''
    area_name = dem_list[0].split('_')[1]
    fileOutName = area_name + '_histogram.pdf' # output file
    sns.set_palette(sns.color_palette("Set1", n_colors=len(dem_list)+2, desat=.5))
        # get data for each dem
    for dem in dem_list:
        dem_1d = raster_as_1d_array(dem)
        elev = dem_1d[~np.isnan(dem_1d)]
        # define bins
        g_max = int(np.ceil(np.max(elev)))
        g_min = int(np.min(elev))
        nbins = round5((g_max - g_min)/bin_width)
        # plot histogram
        # plt.hist(elev, bins=nbins, normed=density_plot, histtype='step', label=dem)
        hist, edges = np.histogram(elev, bins=nbins, density=density_plot)
        plt.plot(edges[:-1], hist, label=dem)
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
def do_histogram_full(dem,bin_width):
    ''' calculates histogram on DEM tiles '''
    area_name = dem.split('_')[-1]
    fileOutName = area_name + '_histogram_full_srtm.pdf' # output file
    dem_1d = raster_as_1d_array(dem)
    elev = dem_1d[~np.isnan(dem_1d)]
    # define bins
    g_max = int(np.ceil(np.max(elev)))
    g_min = int(np.min(elev))
    nbins = round5((g_max - g_min)/bin_width)
    hist, edges = np.histogram(elev, bins=nbins, density=True)
    plt.plot(edges[:-1], hist, label=dem)
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
def do_histogram_slope(dem_list,bin_width,x_lim,density_plot,n_colors):
    ''' calculates histogram on DEM tiles '''
    area_name = dem_list[0].split('_')[1]
    fileOutName = area_name + '_histogram_slope.pdf' # output file
    for dem in dem_list:
        slope_rast = dem + '_slope'
        slope_1d = raster_as_1d_array(slope_rast)
        slope = slope_1d[~np.isnan(slope_1d)]
        sns.set_palette(sns.color_palette("Set1", n_colors=n_colors, desat=.5))
        # slope = slope[~np.isnan(slope)]
        # define bins
        g_max = int(np.ceil(np.max(slope)))
        g_min = int(np.min(slope))
        nbins = round5((g_max - g_min)/bin_width)
        # plot histogram
        # plt.hist(slope, bins=nbins, normed=density_plot, histtype='step', label=dem)
        hist, edges = np.histogram(slope, bins=nbins, density=density_plot)
        plt.plot(edges[:-1], hist, label=slope_rast)
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
def mean_slope_elev(area_list,elev_interval,ymin,ymax):
    ''' plot of mean slope per elevation interval'''
    area_name = area_list[0].split('_')[-1]
    file_svg = area_name + '_mean_slope_elev.pdf'
    sns.set_palette(sns.color_palette("muted", n_colors=len(area_list), desat=.5))
    for dem in area_list:
        # read raster data and put that into a DataFrame
        slope_rast = dem + '_slope'
        df = pd.DataFrame()
        df['elev'] = raster_as_1d_array(dem)
        df['slope'] = raster_as_1d_array(slope_rast)
        # bins
        zmin = int(round(df['elev'].min()/elev_interval)*elev_interval)
        zmax = int(round(df['elev'].max()/elev_interval)*elev_interval)
        bins = range(zmin,zmax+1,int(elev_interval))
        # aggregate slope into breaks
        df['cut'] = pd.cut(df['elev'], bins=bins)
        pvt = pd.pivot_table(df, columns='cut', values='slope', aggfunc=np.mean)
        # plot 
        x = pvt.values[0]
        y = [l.mid for l in list(pvt)]
        print('x = %s, y = %s' %(len(x),len(y)))
        plt.plot(x,y, label=dem)
        df = None
    # plot decorations
    plt.title(area_name + ' - mean slope x elevation ')
    plt.xlabel('Mean slope - interval:' + str(elev_interval) + ' m')
    plt.ylabel('Elevation')
    plt.ylim(ymin, ymax)
    plt.xlim(0,90)
    plt.legend()
    plt.tight_layout()
    plt.savefig(file_svg)
    print 'plot OK'
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

def err_min(x1, x2, axis=0):
    """mean error"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    return np.min((x1-x2), axis=axis)

def err_max(x1, x2, axis=0):
    """mean error"""
    x1 = np.asanyarray(x1)
    x2 = np.asanyarray(x2)
    return np.max((x1-x2), axis=axis)

#----------------------------------------------------
# calculate error metric (rmse, mae, me, std)
def calc_error_metrics(areas):
    ''' error metrics (rmse, mae, me, std) between two raster maps'''
    area_list = []
    dem_list = []
    rmse_list = []
    me_list = []
    std_list = []
    le90_list = []
    le95_list = []
    le99_list = []
    min_list = []
    max_list = []
    for area in areas:
        area_name = area[0].split('_')[-1]
        suffix = area_name + '_12m_bicubic'
        df = pd.DataFrame()
        for dem in area:
            if dem.startswith('tdx12'):
                df[dem] = raster_as_1d_array(dem)
            else:
                dem = dem + '_12m_bicubic'
                df[dem] = raster_as_1d_array(dem)
            print(dem + ' added to DataFrame')
        dfna = df.dropna(axis=0, how='any')
        for column in ['tdx30_'+suffix,'srtm30_wgs84_'+suffix,'aster30_wgs84_'+suffix,'aw3d30_wgs84_'+suffix]:
            # calculate errors
            rmse = err_rmse(dfna['tdx12_'+area_name],dfna[column])
            me = err_me(dfna['tdx12_'+area_name],dfna[column])
            std = err_std(dfna['tdx12_'+area_name],dfna[column])
            le90 = 1.6449 * std
            le95 = 1.9600 * std
            le99 = 3.0000 * std
            min_err = err_min(dfna['tdx12_'+area_name],dfna[column])
            max_err = err_max(dfna['tdx12_'+area_name],dfna[column])
            # append data to lists
            rmse_list.append(rmse)
            me_list.append(me)
            std_list.append(std)
            le90_list.append(le90)
            le95_list.append(le95)
            le99_list.append(le99)
            dem_list.append(column)
            area_list.append(area_name)
            min_list.append(min_err)
            max_list.append(max_err)
        # clean memory
        df = None
        dfna = None
        print('\n\n')
    return rmse_list, me_list, std_list, le90_list, le95_list, le99_list, dem_list, area_list, min_list, max_list


#----------------------------------------------------
# calculates a histogram for each DEM tile and saves it to an svg file
# uses matplotlib histogram
def do_histogram_dod(dem_list,bin_width,density_plot):
    ''' calculates histogram on DEM tiles '''
    area_name = dem_list[0].split('_')[1]
    fileOutName = area_name + '_histograms_DoD.pdf' # output file
    sns.set_palette(sns.color_palette("Set1", n_colors=len(dem_list)+2, desat=.5))
        # get data for each dem
    for dem in dem_list:
        dem_1d = raster_as_1d_array(dem)
        elev = dem_1d[~np.isnan(dem_1d)]
        # define bins
        g_max = int(np.ceil(np.max(elev)))
        g_min = int(np.min(elev))
        nbins = round(g_max - g_min)
        # plot histogram
        # plt.hist(elev, bins=nbins, normed=density_plot, histtype='step', label=dem)
        hist, edges = np.histogram(elev, bins=nbins, density=density_plot)
        plt.plot(edges[:-1], hist, label=dem)
    # plot decorations
    plt.title(area_name + ' (bin width = ' + str(bin_width) +'m)')
    plt.xlabel('Difference in Elevation from TanDEM-X 12m')
    if density_plot == True:
        plt.ylabel('Normalized probability density function')
    else:
        plt.ylabel('Cell count')
    plt.legend()
    plt.tight_layout()
    plt.savefig(fileOutName)
    print 'histogram OK'
    # plt.show()
    plt.clf()
    plt.cla()



#----------------------------------------------------
# scatter plots of original elevations
def do_scatter_elev(elev1,elev2,file_out,area_name,dem_name):
    ''' scatter plots of original elevations '''
    elev_x = raster_as_1d_array(elev1)
    elev_y = raster_as_1d_array(elev2)
    plt.plot(elev_x,elev_y,'+')
    plt.title(area_name + ' - TDX 12m x ' + dem_name)
    plt.xlabel('TDX 12m elevation')
    plt.ylabel(dem_name + ' elevation')
    xmin = plt.gca().get_xlim()[0]
    xmax = plt.gca().get_xlim()[1]
    plt.plot((xmin,xmax),(xmin,xmax))
    # plt.ylim(xmin,xmax)
    # plt.xlim(xmin,xmax)
    # plt.legend()
    plt.tight_layout()
    plt.savefig(file_out)
    print 'scatterplot OK'
    # plt.show()
    plt.clf()
    plt.cla()


#----------------------------------------------------
# scatter plots of aspect x DoD
def do_scatter_aspect(asp_name,dod_name,file_out):
    ''' scatter plots of original elevations '''
    elev_x = raster_as_1d_array(asp_name)
    elev_y = raster_as_1d_array(dod_name)
    fig = plt.figure(figsize=(10,6))
    plt.plot(elev_x,elev_y,'+',color='k',alpha=.2,ms=2)
    plt.title('DoD x aspect')
    plt.xlabel(asp_name)
    plt.ylabel(dod_name)
    plt.xlim(-20,380)
    plt.ylim(-100,100)
    plt.savefig(file_out, dpi=(600))
    print 'scatterplot OK'
    # plt.show()
    plt.clf()
    plt.cla()
    plt.close(fig)


# scatter plots of aspect x DoD with reunning mean
# unused, left here as reference and example
# def do_scatter_aspect_mean(asp_name,dod_name,file_out):
#     ''' scatter plots of original elevations '''
#     # elev_x = raster_as_1d_array(asp_name)
#     # elev_y = raster_as_1d_array(dod_name)
#     bins = np.arange(0,361,5)
#     df = pd.DataFrame({'X':raster_as_1d_array(asp_name), 'Y':raster_as_1d_array(dod_name)})
#     data_cut = pd.cut(df.X,bins)
#     grp = df.groupby(by=data_cut)
#     ret = grp.aggregate(np.median)
#     plt.figure(figsize=(10,7.5))
#     plt.plot(df.X,df.Y,'+',color='k',alpha=.2,ms=2)
#     plt.plot(ret.X,ret.Y,'r-',lw=2,alpha=.8)
#     # plt.plot(elev_x,elev_y,'+')
#     plt.title('DoD x aspect')
#     plt.xlabel(asp_name)
#     plt.ylabel(dod_name)
#     plt.xlim(0,360)
#     plt.ylim(-300,300)
#     plt.savefig(file_out)
#     print 'scatterplot OK'
#     # plt.show()
#     plt.clf()
#     plt.cla()



#----------------------------------------------------
#  make a png of a diff map and plot topo profile
def make_diff_map_png(dem,area,n,s,w,e,shade_map,grid_size,suffix):
    ''' make a png of a diff map and plot topo profile
    tdx: 'tdx12' or 'tdx30'
    area: name of study area (as in 'Barcelos')
    n,s,w,e: limits of GRASS region
    shade_map: boolean
    grid_size: interval of grid on output map
    suffix: so we can have outputs name differently
    '''
    diff_map = 'tdx12_diff_'+dem+'_wgs84_'+area
    diff_map_inset = diff_map+'_inset'+suffix
    shade = 'tdx12_'+area+'_shade_315_20'
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
        grass.run_command('d.mon', start='cairo', output=diff_map_inset+'_shade.png', resolution='3', height=500, width=500, overwrite=True)
        grass.run_command('d.shade', shade=shade, color=diff_map)
    else:
        grass.run_command('d.mon', start='cairo', output=diff_map_inset+'.png', resolution='3', height=500, width=500, overwrite=True)
        grass.run_command('d.rast', map=diff_map_inset)
    grass.run_command('d.grid', size=grid_size, text_color='black', flags='c')
    grass.run_command('d.legend', raster=diff_map_inset, flags='tsbd', at='4,25,92,94', font='Helvetica', fontsize=12, bgcolor='240:240:240')
    grass.run_command('d.mon', stop='cairo')
    print('export as png OK')
    # export simple PNG + world file
    grass.run_command('r.out.png',input=diff_map_inset, output=diff_map_inset+'_wld.png', flags='w', overwrite=True)


def make_diff_profile(area,n,s,w,e,coords,fig_aspect):
    df_tdx = None
    df_srtm = None
    grass.run_command('g.region', n=n, s=s, w=w, e=e, res='0:0:01')
    # make profile
    csv_file_tdx = 'tdx12_'+area+'_profile.csv'
    csv_file_srtm = 'srtm30_'+area+'_profile.csv'
    grass.run_command('r.profile', input='tdx12_'+area, output=csv_file_tdx, coordinates=coords, overwrite=True)
    grass.run_command('r.profile', input='srtm30_wgs84_'+area, output=csv_file_srtm, coordinates=coords, overwrite=True)
    # read csv and plot profile
    df_tdx = pd.read_csv(work_dir + csv_file_tdx, sep=' ', header=None, names=['nan','dist','elev'],na_values='*')
    df_srtm = pd.read_csv(work_dir + csv_file_srtm, sep=' ', header=None, names=['nan','dist','elev'],na_values='*')
    file_svg = area + '_diff_profile.pdf'
    plt.plot(df_tdx['dist'],df_tdx['elev'],'k-',label='TanDEM-X 12m')
    plt.plot(df_srtm['dist'],df_srtm['elev'],'b-',label='SRTM 30m')
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




#----------------------------------------------------
# make contours only in a small area
def make_zoom_contour(area_name,cont_param,n,s,w,e):
    ''' make contours only in a small area'''
    dems = ['aster30_wgs84_','aw3d30_wgs84_','srtm30_wgs84_','tdx30_','tdx12_']
    for dem in dems:
        dem_in = dem + area_name
        grass.run_command('g.region',raster=dem_in, flags='pa')
        dem_zoom = dem_in + '_zoom'
        print('clipping with r.mapcalc...')
        grass.run_command('g.region',n=n,s=s,w=w,e=e, flags='pa')
        grass.mapcalc("${out} = ${orig}",
            out=dem_zoom,
            orig=dem_in,
            overwrite = True)
        vect_contour = dem + area_name + '_contours_zoom'
        interval = cont_param[0]
        min_c = cont_param[1]
        max_c = cont_param[2]
        grass.run_command('r.contour', input=dem_zoom, output=vect_contour, \
            step=interval, minlevel=min_c, maxlevel=max_c, overwrite=True)

# export PNGs of contours with shaded relief
def make_shade_contour(n,s,w,e,tc,gc):
    ''' make pngs of shaded relief maps overlaid by vector contours'''
    dems = ['aster30_wgs84_','aw3d30_wgs84_','srtm30_wgs84_','tdx30_','tdx12_']
    for dem in dems:
        out_shade = dem + area_name + '_shade_315_20'
        vect_contour = dem + area_name + '_contours_zoom'
        grass.run_command('g.region',n=n,s=s,w=w,e=e, flags='pa')
        grass.run_command('d.mon', start='cairo', output=dem+area_name+'_shade_contours.png', \
            resolution='3', height=500, width=500, overwrite=True)
        grass.run_command('d.rast', map=out_shade)
        grass.run_command('d.vect', map=vect_contour, type='line')
        grass.run_command('d.grid', size='0.025', text_color=tc, color=gc, fontsize=16, flags='c')
        grass.run_command('d.mon', stop='cairo')

# export PNGs of colored shaded relief (tdx12)
def make_shadecolor_zoom(area_name,n,s,w,e,tc,gc):
    ''' export PNGs of colored shaded relief (tdx12)'''
    dem = 'tdx12_' + area_name + '_zoom'
    dem_shade = 'tdx12_' + area_name + '_shade_315_20'
    grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:0:00.4', flags='pa')
    grass.run_command('d.mon', start='cairo', output=dem+area_name+'_colorshade.png', \
        resolution='3', height=500, width=500, overwrite=True)
    grass.run_command('d.shade', shade=dem_shade, color=dem)
    grass.run_command('d.grid', size='0.025', text_color=tc, color=gc, fontsize=16, flags='c')
    grass.run_command('d.legend', raster=dem, flags='sb', at='4,25,2,4', font='Helvetica', \
        fontsize=13, bgcolor='240:240:240', title='(m)')
    grass.run_command('d.mon', stop='cairo')

# export PNGs of colored shaded relief (tdx12)
def make_shade_slope_zoom(area_name,n,s,w,e,tc,gc):
    ''' export PNGs of colored shaded relief (tdx12)'''
    dem = 'tdx12_' + area_name
    dem_shade = 'tdx12_' + area_name + '_shade_315_20'
    dem_slope = 'tdx12_' + area_name + '_slope'
    grass.run_command('g.region',n=n,s=s,w=w,e=e, res='0:0:00.4', flags='pa')
    grass.run_command('d.mon', start='cairo', output=dem+'_shade_slope.png', \
        resolution='3', height=500, width=500, overwrite=True)
    grass.run_command('d.rast', map=dem_shade)
    grass.run_command('d.rast', map=dem_slope, values='15-90')
    grass.run_command('d.grid', size='0.025', text_color=tc, color=gc, fontsize=16, flags='c')
    grass.run_command('d.legend', raster=dem_slope, flags='sb', at='4,25,2,4', font='Helvetica', \
        fontsize=13, bgcolor='240:240:240', title='(m)', range=(15,90))
    grass.run_command('d.mon', stop='cairo')


# export PNGs of shaded relief overlaid by WAM mask
def make_shade_wam(area_name,tc,gc):
    ''' export PNGs of colored shaded relief (tdx12)'''
    dem = 'tdx12_' + area_name
    dem_shade = 'tdx12_' + area_name + '_shade_315_20'
    dem_wam = 'tdx12_wam_' + area_name + '_cats'
    grass.run_command('g.region', raster=dem, flags='pa')
    grass.run_command('d.mon', start='cairo', output=dem+'_wam_colorcats.png', \
        resolution='3', height=600, width=600, overwrite=True)
    # grass.run_command('d.shade', shade=dem_shade, color=dem_wam)
    grass.run_command('d.rast', map=dem_shade)
    grass.run_command('d.rast', map=dem_wam)
    grass.run_command('d.grid', size='0.25', text_color=tc, color=gc, fontsize=24, flags='c')
    grass.run_command('d.mon', stop='cairo')
