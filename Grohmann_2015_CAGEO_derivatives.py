#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
#

# Testing behaviour of DEM derivatives
# What is best: to calculate derivatives (slope, aspect, curvatures)
# in high resolution and then resample to small resolution
# or resample the DEM to small resolution and then calculate derivatives?

# test in three distinct landscapes: 
# flat (Amazon), hilly (southeastern Brazil) and mountains (Andes)

# base data: SRTM 03" v.4 from CGIAR

# Carlos H. Grohmann - 2012/2013/2014
# guano (at) usp (dot) br


#----------------------------------------------------
# imports
import sys, os, csv
import math as m
import numpy as np
import scipy.stats as ss
import grass.script as grass
import grass.script.array as garray
import grass.script.vector as gvect
import matplotlib.pyplot as plt
import gc

from rpy2.robjects import r
from rpy2.robjects import FloatVector
from rpy2.robjects import IntVector
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as rpyn
circular = importr('circular')
graphics = importr('graphics')
grdevices = importr('grDevices')

gc.enable()

#----------------------------------------------------
# set working directory
workDir = 'Volumes/HDD/Users/guano/Dropbox/derivadas_dem'
os.chdir(workDir)

# set overwrite flag
ow = True

# three test areas: flat (amazon), hilly (qf) and mountains (andes)
region = ['amazon','qf','andes']


# aux func
def flatParam(parameter, aslist=False):
    print 'flatParam: ' + parameter
    param = garray.array()
    param.read(parameter, null=np.nan)
    param = param[~np.isnan(param)]
    if aslist == True:
        parflat = param.tolist()
    else:
        parflat = np.array(param)
    param = None
    gc.collect()
    return parflat




# g.region region=amazon -pa
# projection: 3 (Latitude-Longitude)
# zone:       0
# datum:      wgs84
# ellipsoid:  wgs84
# north:      1:35:48S
# south:      5:25:12S
# west:       62:40:03W
# east:       57:04:48W
# nsres:      0:00:03
# ewres:      0:00:03
# rows:       4588
# cols:       6705
# cells:      30762540

# g.region region=qf -pa
# projection: 3 (Latitude-Longitude)
# zone:       0
# datum:      wgs84
# ellipsoid:  wgs84
# north:      19:00:48S
# south:      21:08:48S
# west:       45:40:57W
# east:       42:27:57W
# nsres:      0:00:03
# ewres:      0:00:03
# rows:       2560
# cols:       3860
# cells:      9881600

# g.region region=andes -pa
# projection: 3 (Latitude-Longitude)
# zone:       0
# datum:      wgs84
# ellipsoid:  wgs84
# north:      32:07S
# south:      34:46S
# west:       71:43:30W
# east:       67:57:30W
# nsres:      0:00:03
# ewres:      0:00:03
# rows:       3180
# cols:       4520
# cells:      14373600


#----------------------------------------------------
# 1 - import SRTM data
# working with CGIAR SRTM v.4 data, so no need to fill voids

# unzip \*.zip

# get list of SRTM files
dataDir = '/Volumes/HDD/Users/guano/Documents/artigos/derivadas_dem/srtm/'
fileList = os.listdir(dataDir)

# import files
for dem in fileList:
    if dem.endswith('tif'):
        outName = dem.split('.')[0]
        inName = dataDir + dem
        grass.run_command('r.in.gdal', input=inName, output=outName, overwrite=ow)


#----------------------------------------------------
# 2 - patch raster maps, clip data to region
# and resample with bilinear to avoid artifacts in derivatives
# (slight change in resolution necessary)
# amazon
grass.run_command('g.region', region='amazon', res='0:0:04')
patchList = 'srtm_23_13,srtm_24_13'
grass.run_command('r.patch', input=patchList, output='patch_amazon', overwrite=ow)
grass.run_command('r.resamp.interp', input='patch_amazon', output='amazon_bilinear', overwrite=ow)

# qf
grass.run_command('g.region', region='qf', res='0:0:04')
patchList = 'srtm_28_16,srtm_28_17'
grass.run_command('r.patch', input=patchList, output='patch_qf', overwrite=ow)
grass.run_command('r.resamp.interp', input='patch_qf', output='qf_bilinear', overwrite=ow)

# andes
grass.run_command('g.region', region='andes', res='0:0:04')
patchList = 'srtm_22_19,srtm_23_19'
grass.run_command('r.patch', input=patchList, output='patch_andes', overwrite=ow)
grass.run_command('r.resamp.interp', input='patch_andes', output='andes_bilinear', overwrite=ow)


#----------------------------------------------------
# 3 - low-res shaded reliefs for illustration purposes
# amazon
grass.run_command('g.region', region='amazon', res='0:0:10')
grass.run_command('r.shaded.relief', map='patch_amazon', \
    shadedmap='shade_amazon_045_25_z10_10sec', \
    units='meters', zmult=10, altitude=25, azimuth=45, overwrite=ow)

# qf
grass.run_command('g.region', region='qf', res='0:0:10')
grass.run_command('r.shaded.relief', map='patch_qf', \
    shadedmap='shade_qf_045_25_z5_10sec', \
    units='meters', zmult=5, altitude=25, azimuth=45, overwrite=ow)

# andes
grass.run_command('g.region', region='andes', res='0:0:06')
grass.run_command('r.shaded.relief', map='patch_andes', \
    shadedmap='shade_andes_045_25_z2_06sec', \
    units='meters', zmult=2, altitude=25, azimuth=45, overwrite=ow)


#----------------------------------------------------
#----------------------------------------------------
# let's begin. create all maps first, make stats later

region = ['amazon','qf','andes']
# parameters = [slope, aspect]
resList = ['0:0:10','0:0:15','0:0:20','0:0:25','0:0:30','0:0:35', \
        '0:0:40','0:0:45','0:0:50','0:0:55','0:01']


#----------------------------------------------------
# SRTM 03" resolution (original)

# 4 - base derivatives (slope and aspect)
for reg in region:
    grass.run_command('g.region', region=reg)
    dem = reg + '_bilinear'
    slope = reg + '_slope'
    aspect = reg + '_aspect'
    grass.run_command('r.slope.aspect', elevation=dem, slope=slope, \
        aspect=aspect, min_slp_allowed=0.5, overwrite=ow)
    # change aspect orientation to compass
    aspectComp = aspect + '_compass'
    grass.mapcalc("${out} = if(${rast1} < 90 ,90 - ${rast1}, 450 - ${rast1})", \
            out = aspectComp,
            rast1 = aspect)

#----------------------------------------------------
#----------------------------------------------------
# 5 - average base derivatives (param-average = PA)
# vector mean of aspect
region = ['amazon','qf','andes']
for reg in region:
    dem = reg + '_bilinear'
    slope = reg + '_slope'
    aspect = reg + '_aspect_compass'
    cosX = reg + '_cosine_x'
    cosY = reg + '_cosine_y'
    grass.run_command('g.region', region=reg, res='0:0:03', flags='pa')
    grass.mapcalc("${out} = sin(${rast1})",
        out = cosX,
        rast1 = aspect,
        overwrite=ow)
    grass.mapcalc("${out} = cos(${rast1})",
        out = cosY,
        rast1 = aspect,
        overwrite=ow)
    for res in resList:
        grass.run_command('g.region', region=reg, res=res, flags='pa')
        # slope
        avg = slope +'_avg_' + res[-2:] + '_PA'
        grass.run_command('r.resamp.stats', input=slope, output=avg, method='average', overwrite=ow) 
        # aspect - vector mean
        # calculate SUM of direction cosines
        grass.run_command('r.resamp.stats', input=cosX, output='sumX', method='sum', overwrite=ow) 
        grass.run_command('r.resamp.stats', input=cosY, output='sumY', method='sum', overwrite=ow) 
        # mean vector
        outAspect = aspect + '_VectAvg_' + res[-2:] + '_PA'
        grass.mapcalc("${out} = atan(${rast2}, ${rast1})",
        out = outAspect,
        rast1 = 'sumX',
        rast2 = 'sumY',
        overwrite=ow)
    grass.run_command('g.remove', rast=('sumX','sumY'), quiet=True)
    gc.collect()

print 'step 5: param-average (PA) done'


#----------------------------------------------------
#----------------------------------------------------
# Resampled SRTM

# 6 - average DEM then take derivatives (average-param = AP)
region = ['amazon','qf','andes']
for reg in region:
    dem = reg + '_bilinear'
    for res in resList:
        grass.run_command('g.region', region=reg, res=res, flags='pa')
        demAvg = reg + '_avg_' + res[-2:]
        grass.run_command('r.resamp.stats', input=dem, output=demAvg, method='average', overwrite=ow) 
        # calculate parameters for resampled DEM
        slope = reg + '_avg_slope_' + res[-2:] + '_AP'
        aspect = reg + '_avg_aspect_' + res[-2:] + '_AP'
        grass.run_command('r.slope.aspect', elevation=demAvg, slope=slope, \
        aspect=aspect, overwrite=ow)
        # convert aspect to compass-oriented
        aspectComp = reg + '_avg_aspect_compass_' + res[-2:] + '_AP'
        grass.mapcalc("${out} = if(${rast1} < 90 ,90 - ${rast1}, 450 - ${rast1})", \
                out = aspectComp,
                rast1 = aspect,
                overwrite=ow)


print 'step 6: average-param (AP) done'




#----------------------------------------------------
#----------------------------------------------------
# 7 - Summary stats 
# 7.1 - SRTM
# files for results
os.chdir('/Volumes/HDD/Users/guano/Dropbox/artigos/derivadas_dem/stats/descr_stats/')
fileOut = open('estatisticas_descritivas_SRTM.txt', 'w') 
fileOut.write('estatisticas_descritivas \n') 

region = ['andes','amazon','qf']
for reg in region:
        grass.run_command('g.region', region=reg)
        # original DEM
        print 'original DEM'
        demOrig = reg + '_bilinear'
        parfOrig = flatParam(demOrig)
        fileOut.write(demOrig + '\n')
        fileOut.write('min: max: mean: median: stddev: skewness: kurtosis: \n')
            fileOut.write('%s %s %s %s %s %s %s \n') % (np.min(parfOrig), \
                np.max(parfOrig), np.mean(parfOrig), \
                np.median(parfOrig), np.std(parfOrig), \
                ss.stats.skew(parfOrig), ss.stats.kurtosis(parfOrig))
        fileOut.write('\n')
        fileOut.write('\n')
        gc.collect()
        # averaged DEMs
        print 'averaged DEMs'
        for res in resList:
            grass.run_command('g.region', region=reg, res=res, flags='pa')
            demAvg = reg + '_avg_' + res[-2:]
            parfOrig = flatParam(demAvg)
            fileOut.write(demAvg + '\n')
            fileOut.write('min: max: mean: median: stddev: skewness: kurtosis: \n')
            fileOut.write('%s %s %s %s %s %s %s \n' % (np.min(parfOrig), \
                np.max(parfOrig), np.mean(parfOrig), \
                np.median(parfOrig), np.std(parfOrig), \
                ss.stats.skew(parfOrig), ss.stats.kurtosis(parfOrig)))
            fileOut.write('\n')
            fileOut.write('\n')
            gc.collect()
# 

fileOut.close()
print 'step 7.1: summary stats SRTM done'


#----------------------------------------------------
# 7.2 - summary stats - original parameter
# files for results
os.chdir('/Volumes/HDD/Users/guano/Dropbox/artigos/derivadas_dem/stats/descr_stats/')
fileOut = open('estatisticas_descritivas_ParamOrig.txt', 'w') 
fileOut.write('estatisticas_descritivas \n') 

region = ['andes','amazon','qf']
for reg in region:
        grass.run_command('g.region', region=reg)
        for param in ['slope', 'aspect_compass']:
            paramOrig = reg + '_' + param
            print paramOrig
            parfOrig = flatParam(paramOrig)
            fileOut.write(paramOrig + '\n')
            fileOut.write('min: max: mean: median: stddev: skewness: kurtosis: \n')
            fileOut.write('%s %s %s %s %s %s %s \n' % (np.min(parfOrig), \
                np.max(parfOrig), np.mean(parfOrig), \
                np.median(parfOrig), np.std(parfOrig), \
                ss.stats.skew(parfOrig), ss.stats.kurtosis(parfOrig)))
            fileOut.write('\n')
            fileOut.write('\n')
            gc.collect()
# 

fileOut.close()
print 'step 7.2: summary stats original parameters done'


#----------------------------------------------------
# 7.3 - summary stats (average-param = AP)
# files for results
os.chdir('/Volumes/HDD/Users/guano/Dropbox/artigos/derivadas_dem/stats/descr_stats/')
fileOut = open('estatisticas_descritivas_AvgParam.txt', 'w') 
fileOut.write('estatisticas_descritivas \n') 

region = ['andes','amazon','qf']
for reg in region:
    for param in ['slope', 'aspect']:
        for res in resList:
            grass.run_command('g.region', region=reg, res=res, flags='pa')
            if param == 'aspect':
                paramAP = reg + '_avg_aspect_compass_' + res[-2:] + '_AP'
            else:
                paramAP = reg + '_avg_slope_' + res[-2:] + '_AP'
            print paramAP
            parfAP = flatParam(paramAP)
            fileOut.write(paramAP + '\n')
            fileOut.write('min: max: mean: median: stddev: skewness: kurtosis: \n')
            fileOut.write('%s %s %s %s %s %s %s \n' % (np.min(parfAP), \
                np.max(parfAP), np.mean(parfAP), \
                np.median(parfAP), np.std(parfAP), \
                ss.stats.skew(parfAP), ss.stats.kurtosis(parfAP)))
            fileOut.write('\n')
            fileOut.write('\n')
            gc.collect()
# 

fileOut.close()
print 'step 7.3: summary stats average-param AP done'


#----------------------------------------------------
# 7.4 - summary stats (param-average = PA)
# files for results
os.chdir('/Volumes/HDD/Users/guano/Dropbox/artigos/derivadas_dem/stats/descr_stats/')
fileOut = open('estatisticas_descritivas_ParamAvg.txt', 'w') 
fileOut.write('estatisticas_descritivas \n') 

region = ['andes','amazon','qf']
for reg in region:
    for param in ['slope', 'aspect']:
        for res in resList:
            grass.run_command('g.region', region=reg, res=res, flags='pa')
            # for method in ['avg', 'median', 'mode']:
            if param == 'aspect':
                paramPA = reg + '_aspect_compass_VectAvg_' + res[-2:] + '_PA'
            else:
                paramPA = reg + '_slope_avg_' + res[-2:] + '_PA'
            print paramPA
            parfPA = flatParam(paramPA)
            fileOut.write(paramPA + '\n')
            fileOut.write('min: max: mean: median: stddev: skewness: kurtosis: \n')
            fileOut.write('%s %s %s %s %s %s %s \n' % (np.min(parfPA), \
                np.max(parfPA), np.mean(parfPA), \
                np.median(parfPA), np.std(parfPA), \
                ss.stats.skew(parfPA), ss.stats.kurtosis(parfPA)))
            fileOut.write('\n')
            fileOut.write('\n')
            gc.collect()
# 

fileOut.close()
print 'step 7.4: summary stats param-average PA done'



#----------------------------------------------------
#----------------------------------------------------
# 8 - Analysis

#----------------------------------------------------
#----------------------------------------------------
# density plots
os.chdir('/Volumes/HDD/Users/guano/Dropbox/artigos/derivadas_dem/stats/')

# 8.1 - SRTM elevation
region = ['andes','amazon','qf']
for reg in region:
        grass.run_command('g.region', region=reg)
        fig = reg + '_SRTM_density.svg'
        print 'fig will be: ' + fig
        # original DEM        
        demOrig = reg + '_bilinear'
        parfOrig = flatParam(demOrig)
        print 'flatParam for ' + demOrig + ' OK'
        densityOrig = ss.kde.gaussian_kde(parfOrig)
        gc.collect()
        print 'KDE for ' + demOrig + ' OK'
        # determine steps for density plot (to save time)
        if reg == 'andes': 
            step = 30.
        elif reg == 'amazon':
            step = 1.
        elif reg == 'qf':
            step = 20.
        x = np.arange(np.min(parfOrig), np.max(parfOrig), step)
        yOrig = densityOrig(x)
        plt.ylabel('density')
        plt.xlabel('elevation')
        plt.plot(x, yOrig, label=demOrig)
        print 'plot for ' + demOrig + ' OK'
        gc.collect()
        # averaged DEMs
        for res in resList[::-1]:
            grass.run_command('g.region', region=reg, res=res, flags='pa')
            demAvg = reg + '_avg_' + res[-2:]
            parfAvg = flatParam(demAvg)
            print 'flatParam for ' + demAvg + ' OK'
            densityAvg = ss.kde.gaussian_kde(parfAvg)
            gc.collect()
            print 'KDE for ' + demAvg + ' OK'
            x = np.arange(np.min(parfAvg), np.max(parfAvg), step)
            yAvg = densityAvg(x)
            plt.plot(x, yAvg, label=demAvg)
            print 'plot for ' + demAvg + ' OK'
            gc.collect()
        plt.legend()
        plt.savefig(fig)
        plt.clf()



#----------------------------------------------------
# paramOrig
# 8.2 - density plots: original slope
region = ['andes','amazon','qf']
for reg in region:
    grass.run_command('g.region', region=reg)
    fig = reg + '_slope_density_original.svg'
    print 'fig will be: ' + fig
    paramOrig = reg + '_slope'
    parfOrig = flatParam(paramOrig)
    print 'flatParam for ' + paramOrig + ' OK'
    densityOrig = ss.kde.gaussian_kde(parfOrig)
    gc.collect()
    print 'KDE for ' + paramOrig + ' OK'
    x = np.arange(0., np.max(parfOrig), .5)
    yOrig = densityOrig(x)
    plt.ylabel('density')
    plt.xlabel('slope')
    plt.plot(x, yOrig, label=paramOrig)
    print 'plot for ' + paramOrig + ' OK'
    gc.collect()
    plt.legend()
    plt.savefig(fig)
    plt.clf()



# 8.3 - density plots: original aspect
# use circular function from R
# NOTE: 0:0:03 resolution was too much for the MBP with 8GB RAM,
# the andes region could be calculated in 3 hours but the amazon
# was too much. re-calculated with 0:0:06
# vector mean of aspect
# calculate direction cosines
# u = sin(alpha) # x
# v = cos(alpha) # y
region = ['andes','amazon','qf']
for reg in region:
    grass.run_command('g.region', region=reg, flags='pa', res='0:0:06')
    paramOrig = reg + '_aspect_compass'
    parfOrig = flatParam(paramOrig, aslist=True)
    print 'flatParam for ' + paramOrig + ' OK'
    # create R object for circular density
    parfOrig_r = IntVector(parfOrig)
    parfOrig_rc = r.circular(parfOrig_r, units='degrees', template='geographics')
    # optimal bandwith determined based on a VonMises distribution
    # using circular.bw_nrd_circular()
    if reg == 'andes':
        bandwith = 2.0
    elif reg == 'qf':
        bandwith = 3.5
    elif reg == 'amazon':
        bandwith = 1.0
    # calculate circular density
    print 'calculating density for ' + paramOrig + ', please wait...'
    dens = r.density(parfOrig_rc, bw=bandwith, kernel='vonmises')
    print 'density for ' + paramOrig + ' OK'
    # get density values back to python
    densX = rpyn.ri2numpy(dens[1])
    densY = rpyn.ri2numpy(dens[2])
    gc.collect()
    # densityOrig = ss.kde.gaussian_kde(parfOrig)
    # x = np.arange(0., np.max(parfOrig), .1)
    # yOrig = densityOrig(x)
    # plot
    # cartesian plot
    fig = reg + '_aspect_compass_original_density_circular.svg'
    print 'fig will be: ' + fig
    plt.ylabel('density')
    plt.xlabel('aspect_compass')
    plt.plot(densX, densY, label=paramOrig)
    plt.xticks((-270, -180, -90, 0, 90), ('E','S','W','N','E'))
    print 'cartesian plot for ' + paramOrig + ' OK'
    gc.collect()
    plt.legend()
    plt.savefig(fig)
    plt.clf()
    # polar plot
    fig = reg + '_aspect_compass_original_density_circular_polar.svg'
    print 'fig will be: ' + fig
    ax = plt.subplot(111, polar=True)
    ax.plot(np.radians(densX), densY, label=paramOrig)
    print 'cartesian plot for ' + paramOrig + ' OK'
    plt.legend()
    plt.savefig(fig)
    plt.clf()




# cleanup temp mapcalc maps
grass.run_command('g.remove', rast=('sumX','sumY'), quiet=True)
gc.collect()



#----------------------------------------------------
# paramPA - SLOPE
# 8.4 - density plots (param-average = PA)
region = ['andes','amazon','qf']
for reg in region:
    for res in resList[::-1]:
        grass.run_command('g.region', region=reg, res=res, flags='pa')
        paramPA = reg + '_slope_avg_' + res[-2:] + '_PA'
        fig = reg + '_slope_avg_PA_density.svg'
        print 'fig will be: ' + fig
        parfPA = flatParam(paramPA)
        print 'flatParam for ' + paramPA + ' OK'
        densityPA = ss.kde.gaussian_kde(parfPA)
        print 'KDE for ' + paramPA + ' OK'
        gc.collect()
        x = np.arange(0., np.max(parfPA), .2)
        yPA = densityPA(x)
        plt.ylabel('density')
        plt.xlabel('slope')
        plt.plot(x, yPA, label=paramPA)
        print 'plot for ' + paramPA + ' OK'
        gc.collect()
    plt.legend()
    plt.savefig(fig)
    plt.clf()
# 



#----------------------------------------------------
# paramPA - ASPECT
# 8.5 - density plots (param-average = PA)
region = ['andes','amazon','qf']
for reg in region:
    for res in resList[::-1]:
        grass.run_command('g.region', region=reg, res=res, flags='a')
        paramPA = reg + '_aspect_compass_VectAvg_' + res[-2:] + '_PA'
        parfPA = flatParam(paramPA)
        print 'flatParam for ' + paramPA + ' OK'
        # create R object for circular density
        parfPA_r = IntVector(parfPA)
        parfPA_rc = r.circular(parfPA_r, units='degrees', template='geographics')
        # optimal bandwith determined based on a VonMises distribution 
        # bandwith = circular.bw_nrd_circular(parfPA_rc)
        # print '%s %s %s \n' % (reg, res, bandwith)
        # calculate circular density
        # optimal bandwith determined based on a VonMises distribution
        # using circular.bw_nrd_circular()
        if reg == 'andes':
            bandwith = 2.0
        elif reg == 'qf':
            bandwith = 3.5
        elif reg == 'amazon':
            bandwith = 1.0
        print 'calculating density for ' + paramPA + ', please wait...'
        dens = r.density(parfPA_rc, bw=bandwith, kernel='vonmises')
        print 'density for ' + paramPA + ' OK'
        # get density values back to python
        densX = rpyn.ri2numpy(dens[1])
        densY = rpyn.ri2numpy(dens[2])
        gc.collect()
        # cartesian plot
        fig = reg + '_aspect_compass_VectAvg_PA_density_circular_bw.svg'
        print 'fig will be: ' + fig
        plt.ylabel('density')
        plt.xlabel('aspect_compass')
        plt.plot(densX, densY, label=paramPA)
        print 'plot for ' + paramPA + ' OK'
    plt.legend()
    plt.xticks((-270, -180, -90, 0, 90), ('E','S','W','N','E'))
    plt.savefig(fig)
    plt.clf()
    gc.collect()



#----------------------------------------------------
# paramAP - slope
# 8.6 - density plots (average-param = AP)
region = ['andes','amazon','qf']
for reg in region:
    grass.run_command('g.region', region=reg)
    param = 'slope'
    fig = reg + '_slope_AP_density.svg'
    print 'fig will be: ' + fig
    for res in resList[::-1]:
        grass.run_command('g.region', region=reg, res=res, flags='pa')
        paramAP = reg + '_avg_slope_' + res[-2:] + '_AP'
        parfAP = flatParam(paramAP)
        print 'flatParam for ' + paramAP + ' OK'
        densityAP = ss.kde.gaussian_kde(parfAP)
        gc.collect()
        x = np.arange(0., np.max(parfAP), .1)   
        yAP = densityAP(x)
        plt.ylabel('density')
        plt.xlabel(param)
        print 'density for ' + paramAP + ' OK'
        plt.plot(x, yAP, label=paramAP)
        print 'plot for ' + paramAP + ' OK'
        gc.collect()
    plt.legend()
    plt.savefig(fig)
    plt.clf()


#----------------------------------------------------
# paramAP  - aspect compass
# 8.7 - density plots (average-param = AP)
region = ['andes','amazon','qf']
for reg in region:
    for res in resList[::-1]:
        grass.run_command('g.region', region=reg, res=res, flags='pa')
        paramAP = reg + '_avg_aspect_compass_' + res[-2:] + '_AP'
        parfAP = flatParam(paramAP)
        print 'flatParam for ' + paramAP + ' OK'
        # create R object for circular density
        parfAP_r = IntVector(parfAP)
        parfAP_rc = r.circular(parfAP_r, units='degrees', template='geographics')
        # optimal bandwith determined based on a VonMises distribution 
        # bandwith = circular.bw_nrd_circular(parfAP_rc)
        # calculate circular density
        if reg == 'andes':
            bandwith = 2.0
        elif reg == 'qf':
            bandwith = 3.5
        elif reg == 'amazon':
            bandwith = 1.0
        print 'calculating density for ' + paramAP + ', please wait...'
        dens = r.density(parfAP_rc, bw=bandwith, kernel='vonmises')
        print 'density for ' + paramAP + ' OK'
        # get density values back to python
        densX = rpyn.ri2numpy(dens[1])
        densY = rpyn.ri2numpy(dens[2])
        gc.collect()
        # cartesian plot
        fig = reg + '_aspect_compass_avg_AP_density_circular_bw.svg'
        print 'fig will be: ' + fig
        plt.ylabel('density')
        plt.xlabel('aspect_compass')
        plt.plot(densX, densY, label=paramAP)
        print 'plot for ' + paramAP + ' OK'
    plt.legend()
    plt.xticks((-270, -180, -90, 0, 90), ('E','S','W','N','E'))
    plt.savefig(fig)
    plt.clf()
    gc.collect()





#----------------------------------------------------
#----------------------------------------------------
# 9 - correlation

# regions a bit smaller for random points, 
# to be sure nothing falls out of the computational region

# GRASS 6.4.3 (america_sul):~ > g.region n=1:20S s=2:58s w=67:28w e=64:57w res=0:0:03 save=amazon_rdn -pa
# projection: 3 (Latitude-Longitude)
# zone:       0
# datum:      wgs84
# ellipsoid:  wgs84
# north:      1:19:57S
# south:      2:58S
# west:       67:28W
# east:       64:57W
# nsres:      0:00:03
# ewres:      0:00:03
# rows:       1961
# cols:       3020
# cells:      5922220

# GRASS 6.4.3 (america_sul):~ > g.region n=19:22S s=20:42s w=44:58w e=43:02w res=0:0:03 save=qf_rdn -pa
# projection: 3 (Latitude-Longitude)
# zone:       0
# datum:      wgs84
# ellipsoid:  wgs84
# north:      19:22S
# south:      20:42S
# west:       44:58W
# east:       43:01:57W
# nsres:      0:00:03
# ewres:      0:00:03
# rows:       1600
# cols:       2321
# cells:      3713600

# GRASS 6.4.3 (america_sul):~ > g.region n=32:32S s=34:18s w=70:38w e=69:07w res=0:0:03 save=andes_rdn -pa
# projection: 3 (Latitude-Longitude)
# zone:       0
# datum:      wgs84
# ellipsoid:  wgs84
# north:      32:31:57S
# south:      34:18S
# west:       70:38W
# east:       69:06:57W
# nsres:      0:00:03
# ewres:      0:00:03
# rows:       2121
# cols:       1821
# cells:      3862341


#----------------------------------------------------
# 9.1 - create a set of random points for correlation (10.000 points)

# columns for sampling rasters
cols = 'cat INTEGER, orig DOUBLE, r10pa DOUBLE, r15pa DOUBLE, \
        r20pa DOUBLE, r25pa DOUBLE, r30pa DOUBLE, r35pa DOUBLE, \
        r40pa DOUBLE, r45pa DOUBLE, r50pa DOUBLE, r55pa DOUBLE, \
        r01pa DOUBLE, r10ap DOUBLE, r15ap DOUBLE, r20ap DOUBLE, \
        r25ap DOUBLE, r30ap DOUBLE, r35ap DOUBLE, r40ap DOUBLE, \
        r45ap DOUBLE, r50ap DOUBLE, r55ap DOUBLE, r01ap DOUBLE'

resList = ['0:0:10','0:0:15','0:0:20','0:0:25','0:0:30', \
        '0:0:35', '0:0:40','0:0:45','0:0:50','0:0:55','0:01']


region = ['amazon_rdn','qf_rdn','andes_rdn']
for regSmall in region:
    regLarge = regSmall[:-4]
    for param in ['aspect_compass','slope']:
        # make a region a bit smaller to create random points inside it
        grass.run_command('g.region', region=regSmall, flags='pa')
        # create points and add columns
        pts = regLarge + '_' + param + '_random'
        print '------------------------------------'
        print 'creating random points...'
        print '------------------------------------'
        grass.run_command('v.random', output=pts, n=10000, overwrite=ow)
        print '------------------------------------'
        print 'creating columns in table...'
        print '------------------------------------'
        grass.run_command('v.db.addtable', map=pts, layer=1, columns=cols)
        # sample raster of original derivatives
        # go back to region a bit larger to ensure all points are within
        grass.run_command('g.region', region=regLarge, flags='pa')
        mapOrig = regLarge + '_' + param 
        print '------------------------------------'
        print 'sampling original derivative raster ' + mapOrig
        print '------------------------------------'
        grass.run_command('v.what.rast', vector=pts, \
            raster=mapOrig, layer=1, column='orig') 
        # sample averaged rasters
        for kind in ['PA', 'AP']:
            for res in resList:
                # grass.run_command('g.region', region=regLarge, res=res, flags='pa')
                if kind == 'PA':
                    if param == 'aspect_compass':
                        map = regLarge + '_aspect_compass_VectAvg_' + res[-2:] + '_PA' 
                    else:
                        map = regLarge + '_slope_avg_' + res[-2:] + '_PA'
                elif kind == 'AP':
                    if param == 'aspect_compass':
                        map = regLarge + '_avg_aspect_compass_' + res[-2:] + '_AP'
                    else:
                        map = regLarge + '_avg_slope_' + res[-2:] + '_AP'
                col = 'r' + res[-2:] + kind.lower()
                print '------------------------------------'
                print 'sampling averaged (PA + AP) derivative raster ' + map 
                print '------------------------------------'
                grass.run_command('v.what.rast', vector=pts, \
                raster=map, layer=1, column=col)


#----------------------------------------------------
# 9.2 - perform some maintenance. remove lines that have NULL entries
colList = ['orig', 'r10pa', 'r15pa', 'r20pa', 'r25pa', 'r30pa', \
        'r35pa','r40pa', 'r45pa', 'r50pa', 'r55pa', 'r01pa', \
        'r10ap', 'r15ap', 'r20ap', 'r25ap', 'r30ap', 'r35ap', \
        'r40ap', 'r45ap', 'r50ap', 'r55ap', 'r01ap'] 

region = ['amazon','qf','andes']
for reg in region:
    for param in ['aspect_compass','slope']:
        vect = reg + '_' + param + '_random'
        print 'maintenance for vector ' + vect
        for col in colList:
            # get column from db as a dict
            col_dict = gvect.vector_db_select(map=vect, columns=col)['values']
            # get cats of NULL entries
            null_list = [int(i[1]) for i in col_dict.values() if i[0]=='']
            print 'removing NULL entries...'
            for n in null_list:
                grass.write_command("db.execute", \
                    stdin="DELETE FROM %s WHERE cat = %d" % (vect,n))
        grass.db_describe(vect)


#----------------------------------------------------
# 9.3 - correlation for slope (linear)

# files for results
os.chdir('/Volumes/HDD/Users/guano/Dropbox/artigos/derivadas_dem/stats')
fileOut = open('correlacao_PA_AP_original_slope.txt', 'w') 
fileOut.write('coeficientes de correlacao - SLOPE \n') 

colListPA = ['r10pa', 'r15pa', 'r20pa', 'r25pa', 'r30pa', \
        'r35pa', 'r40pa', 'r45pa', 'r50pa', 'r55pa', 'r01pa']

colListAP = ['r10ap', 'r15ap', 'r20ap', 'r25ap', 'r30ap', \
        'r35ap', 'r40ap', 'r45ap', 'r50ap', 'r55ap', 'r01ap'] 

resolutions = ['0:0:10','0:0:15','0:0:20','0:0:25','0:0:30','0:0:35', \
        '0:0:40','0:0:45','0:0:50','0:0:55','0:01']

region = ['amazon','qf','andes']
for reg in region:
    vect = reg + '_slope_random'
    # get values from 'orig' column
    col_dict = gvect.vector_db_select(map=vect, columns='orig')['values']
    origL = [float(i[0]) for i in col_dict.values()]
    # correlation for average-param (PA)
    corrcoef=[]
    for col in colListPA:
        # get column from db as a dict
        col_dict = gvect.vector_db_select(map=vect, columns=col)['values']
        # values as list
        paramL = [float(i[0]) for i in col_dict.values()]
        # correlation
        m, intercept, r_value, p_value, std_err = ss.linregress(paramL, origL)
        corrcoef.append(r_value**2)
        print corrcoef
        fileOut.write(reg + ' slope ' + col + '\n')
        fileOut.write('slope: %s \n' % (m))
        fileOut.write('intercept: %s \n' % (intercept))
        fileOut.write('R %s \n' % (r_value))
        fileOut.write('R2: %s \n' % (r_value**2))
        fileOut.write('\n')
        fileOut.write('\n')                
        gc.collect()
    # plot for PA
    fig = reg + '_slope_correlation_PA_orig.svg'
    print 'fig will be: ' + fig
    plt.plot(range(len(corrcoef)), corrcoef, label=vect + '_PA')
    plt.xticks(range(len(corrcoef)), resolutions)
    plt.legend()
    plt.savefig(fig)
    plt.clf()
    # correlation for average-param (PA)
    corrcoef=[]
    for col in colListAP:
        # get column from db as a dict
        col_dict = gvect.vector_db_select(map=vect, columns=col)['values']
        # values as list
        paramL = [float(i[0]) for i in col_dict.values()]
        # correlation
        m, intercept, r_value, p_value, std_err = ss.linregress(paramL, origL)
        corrcoef.append(r_value**2)
        print corrcoef
        fileOut.write(reg + ' slope ' + col + '\n')
        fileOut.write('slope: %s \n' % (m))
        fileOut.write('intercept: %s \n' % (intercept))
        fileOut.write('R %s \n' % (r_value))
        fileOut.write('R2: %s \n' % (r_value**2))
        fileOut.write('\n')
        fileOut.write('\n')                
        gc.collect()
    # plot for AP
    fig = reg + '_slope_correlation_AP_orig.svg'
    print 'fig will be: ' + fig
    plt.plot(range(len(corrcoef)), corrcoef, label=vect + '_AP')
    plt.xticks(range(len(corrcoef)), resolutions)
    plt.legend()
    plt.savefig(fig)
    plt.clf()


fileOut.close()



#----------------------------------------------------
# 9.4 - correlation for aspect (circular)

# files for results
os.chdir('/Volumes/HDD/Users/guano/Dropbox/artigos/derivadas_dem/stats')
fileOut = open('correlacao_PA_AP_original_aspect2.txt', 'w') 
fileOut.write('coeficientes de correlacao - ASPECT \n') 

colListPA = ['r10pa', 'r15pa', 'r20pa', 'r25pa', 'r30pa', \
        'r35pa', 'r40pa', 'r45pa', 'r50pa', 'r55pa', 'r01pa']

colListAP = ['r10ap', 'r15ap', 'r20ap', 'r25ap', 'r30ap', \
        'r35ap', 'r40ap', 'r45ap', 'r50ap', 'r55ap', 'r01ap'] 

resolutions = ['0:0:10','0:0:15','0:0:20','0:0:25','0:0:30','0:0:35', \
        '0:0:40','0:0:45','0:0:50','0:0:55','0:01']

region = ['amazon','qf','andes']
for reg in region:
    vect = reg + '_aspect_compass_random'
    # get values from 'orig' column
    col_dict = gvect.vector_db_select(map=vect, columns='orig')['values']
    origL = [float(i[0]) for i in col_dict.values()]
    # R object and circular
    origL_r = FloatVector(origL)
    origL_rc = r.circular(origL_r, units='degrees', template='geographics')
    # correlation for average-param (PA)
    corrcoef=[]
    for col in colListPA:
        # get column from db as a dict
        col_dict = gvect.vector_db_select(map=vect, columns=col)['values']
        # values as list
        paramL = [float(i[0]) for i in col_dict.values()]
        paramL_r = FloatVector(paramL)
        paramL_rc = r.circular(paramL_r, units='degrees', template='geographics')
        # circular correlation
        coef = circular.cor_circular(paramL_rc, origL_rc, test='FALSE')
        corrcoef.append(coef[0])
        print coef
        fileOut.write(reg + ' aspect_compass ' + col + '\n')
        fileOut.write('coef: %s \n' % (coef[0]))
        fileOut.write('\n')
        fileOut.write('\n')                
        gc.collect()
    # plot for PA
    fig = reg + '_aspect_compass_correlation_PA_orig.svg'
    print 'fig will be: ' + fig
    plt.plot(range(len(corrcoef)), corrcoef, label=vect + '_PA')
    plt.xticks(range(len(corrcoef)), resolutions)
    plt.legend()
    plt.savefig(fig)
    plt.clf()
    # correlation for average-param (PA)
    corrcoef=[]
    for col in colListAP:
        # get column from db as a dict
        col_dict = gvect.vector_db_select(map=vect, columns=col)['values']
        # values as list
        paramL = [float(i[0]) for i in col_dict.values()]
        paramL_r = FloatVector(paramL)
        paramL_rc = r.circular(paramL_r, units='degrees', template='geographics')
        # circular correlation
        coef = circular.cor_circular(paramL_rc, origL_rc, test='FALSE')
        corrcoef.append(coef[0])
        print coef
        fileOut.write(reg + ' aspect_compass ' + col + '\n')
        fileOut.write('coef: %s \n' % (coef[0]))
        fileOut.write('\n')
        fileOut.write('\n')                
        gc.collect()
    # plot for AP
    fig = reg + '_aspect_compass_correlation_AP_orig.svg'
    print 'fig will be: ' + fig
    plt.plot(range(len(corrcoef)), corrcoef, label=vect + '_AP')
    plt.xticks(range(len(corrcoef)), resolutions)
    plt.legend()
    plt.savefig(fig)
    plt.clf()


fileOut.close()




