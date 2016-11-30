#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
# Global DEMs (GDEMS) analysis and comparison
# comparing: 
# 30-sec resolution:
    # SRTM30_PLUS v9
    # don't use SRTM30 as it is v2.1 and it's the same used by SRTM30_PLUS
    # GLOBE
    # GMTED2010
    # ACE2
# lower res:
    # ETOPO1
    # ETOPO2
    # ETOPO5
    # Terdat (legacy)
    # TerrainBase (legacy)

# Carlos H. Grohmann - 2015
# guano (at) usp (dot) br
# http://carlosgrohmann.com


#----------------------------------------------------
# imports
import sys, os, csv
import random
import math as m
import numpy as np
import scipy as sp
import scipy.stats as ss
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statsmodels.tools import eval_measures as em
from statsmodels.base.model import GenericLikelihoodModel

import gc
gc.enable()

import grass.script as grass
import grass.script.array as garray
# import grass.script.vector as gvect


#----------------------------------------------------
# create float version of gdem (needed for histograms)
gdems_list = ['gdem_srtm30_plus_v9', 'gdem_globe', 'gdem_gtopo30', \
'gdem_gmted_land', 'gdem_ace2_int','gdem_etopo1_bedrock', \
'gdem_etopo1_ice', 'gdem_etopo2', 'gdem_etopo5', \
'gdem_terdat_round', 'gdem_terrainBase',]

for gdem in gdems_list:
    grass.run_command('g.region', rast=gdem, flags='pa')
    grass.run_command('g.region', n='90N', s='90S', w='180W', e='180E', flags='pa')
    grass.mapcalc("${out} = float(${rast1})",
        out = gdem + '_float',
        rast1 = gdem,
        overwrite = True)
    gdem = gdem + '_float'


#----------------------------------------------------
# helper func: return GRASS raster as flat list 
# (or numpy array), removing the nulls
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


#----------------------------------------------------
# helper func: round to nearest 5
def round5(x):
    rounded = int(round(x/5.0)*5.0)
    return rounded


#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
# this function runs the analysis per se.
# inputs:
# gdem: the gdem being analysed
# mapset: GRASS mapset of analysed gdem (in case its being run in another mapset)
# worldArea: text for plots and files to indicate which area is being analysed 
# xl1, xl2: x-axis limits for histograms
# yl1, yl2: y-axis limits for histograms
# use_stats: what kind of stats should be used (r.info or r.univar)
# bin_methd: boolean, 'fd' for Freedman–Diaconis rule or 'man' for manual method
# writeStats: boolean, write stats to an output file
# makeBoxPlot: boolean, produces a boxplot of the data
# makeHistogram: boolean, produces a histogram of the data (area-elevation)
# 
def gdemStats(gdem, mapset, worldArea, xl1, xl2, yl1, yl2, use_stats, bin_methd, writeStats, makeBoxPlot, makeHistogram):
    '''calculate stats for Global DEMS
        print stats to out file,
        make boxplots and histograms'''
    gdem = gdem + mapset
    print gdem
    grass.run_command('g.region', rast=gdem, flags='pa')
    gc.collect()
    # fast min/max
    # using F-D rule needs univar, so: 
    # bin_methd='fd' MUST go with use_stats='univar'
    if use_stats == 'rinfo':
        rinfo = grass.parse_command('r.info', map=gdem, flags='r')
        g_min = float(rinfo['min'])
        g_max = float(rinfo['max'])
    elif use_stats == 'univar':
    # summary stats, write to file
        print 'calculating stats...'
        univar=None
        univar = grass.parse_command('r.univar', map=gdem, flags= 'ge', percentile=100)
        print 'stats OK...'
        g_n = int(univar['cells'])
        g_null = int(univar['null_cells'])
        g_range = int(univar['range'])
        g_min = float(univar['min'])
        g_max = float(univar['max'])
        g_mean = float(univar['mean'])
        g_median = float(univar['median'])
        g_std = float(univar['stddev'])
        g_var = float(univar['variance'])
        g_cv = float(univar['coeff_var'])
        g_p25 = float(univar['first_quartile'])
        g_p75 = float(univar['third_quartile'])
    if bin_methd == 'fd':
        # Freedman–Diaconis rule to calculate optimal bin size
        # using non-null number of cells
        IQR = g_p75 - g_p25
        bin_width = 2. * IQR * (g_n-g_null)**(-1./3.)
        bin_width = round5(bin_width)
        nbins = round5((g_max - g_min)/bin_width)
        print 'bin width (Freedman–Diaconis rule): ', bin_width
        print 'number of bins: ', nbins
    elif bin_methd == 'man':
        # size of bins after manual (visual) inspection
        # 10min:150m, 5min:75m, 2min:50m, 1min:25m, 30sec:15m
        res_bins = {'0.166666666666667':150, '0.0833333333333':75, '0.0333333333333':50, '0.0166666666667':25, '0.00833333333333':15}
        for key in res_bins:
            gr = str(grass.region()['ewres'])
            if key.startswith(gr[0:5]):
                bin_width = res_bins[key]
        nbins = round5((g_max - g_min)/bin_width)
    if writeStats == True:
        fileOut.write(gdem + '-- ' + worldArea + '  \n')
        fileOut.write('min: max: mean: median: stddev: p25: p75: bin_width: nbins:\n')
        fileOut.write('%s %s %s %s %s %s %s %s %s\n' % \
            (g_min, g_max, g_mean, g_median, g_std, g_p25, g_p75, bin_width, nbins))
        fileOut.write('\n')
        fileOut.write('\n')
        gc.collect()
        print 'write stats OK'
    # plots:
    # make boxplot
    if makeBoxPlot == True:
        # boxplot
        # let's fake some data and change 
            # the boxplot values to real values calculated with numpy.
        # fake data:
        x1 = np.random.normal(0,5000,100)
        bp = plt.boxplot(x1, sym='', whis=1000)
        # # change boxplot values
        bp['medians'][0].set_ydata(np.array([g_median, g_median]))
        bp['whiskers'][0].set_ydata(np.array([g_p25, g_min]))
        bp['whiskers'][1].set_ydata(np.array([g_p75, g_max]))
        bp['caps'][0].set_ydata(np.array([g_min, g_min]))
        bp['caps'][1].set_ydata(np.array([g_max, g_max]))
        bp['boxes'][0].set_ydata(np.array([g_p25, g_p25, g_p75, g_p75, g_p25]))
        plt.ylim(yl1,yl2)
        plt.yticks(range(yl1,yl2+1000,2000))
        # # plt.show()
        figbox = 'boxplot_' + gdem + '_' + worldArea + '.svg'
        plt.savefig(figbox)
        plt.clf()
        print 'boxplot OK'
    # histogram
    if makeHistogram == True:
        # adapted from an example by Michael Barton (2008)
        print 'histogramming....'
        mapstats = grass.read_command('r.stats', input = gdem, \
            separator = ',', nsteps = nbins, flags = 'an')
        histlist = mapstats.splitlines()
        elevlist = []
        arealist = []
        # list of separators, for float maps
        listsep = ['0-', '1-', '2-', '3-', '4-', '5-', '6-', '7-', '8-', '9-']
        for pair in histlist:
            try:
                pairlist = pair.split(',')
                elevRange = pairlist[0]
                area = float(pairlist[1])/1000000
                for elem in listsep:
                    if elem in elevRange:
                        sep = elevRange.index(elem)+1
                        r1 = float(elevRange[0:sep])
                        r2 = float(elevRange[sep+1::])
                        elev = r1 + (bin_width/2) #center of bin
                elevlist.append(elev)
                arealist.append(area)
            except:
                pass
        plt.plot(elevlist,arealist)
        plt.xlabel('Elevation')
        plt.ylabel('Area km2')
        plt.title(gdem + '_' + worldArea + ' (' + str(bin_width) + 'm bins)')
        plt.xlim(xl1,xl2)
        fighist = 'histogram_' + gdem + '_' + worldArea + '_area_ok.svg'
        plt.savefig(fighist)
        plt.clf()
        print 'histogram OK'
    gdem=None
    gc.collect()


#----------------------------------------------------
#----------------------------------------------------
# Summary stats - WORLD

ow = True
grass.run_command('r.mask', flags='r', overwrite=ow)

# list of GDEMS 
gdems_list = ['gdem_srtm30_plus_v9', 'gdem_globe', 'gdem_gtopo30', \
'gdem_gmted_land', 'gdem_ace2_int','gdem_etopo1_bedrock', \
'gdem_etopo1_ice', 'gdem_etopo2', 'gdem_etopo5', \
'gdem_terdat_round', 'gdem_terrainBase',]

# files for results
workDir = '/Volumes/Macintosh HD2T/Dropbox/artigos/global_DEMs/stats/stats_plots_WORLD'
os.chdir(workDir)
fileOut = open('gdems_descriptive_stats_world.txt', 'w') 
fileOut.write('gdem_descriptive_stats_world \n') 

# stats and plots
for gdem in gdems_list[::-1]:
    gdem = gdem + '_float' 
    gdemStats(gdem, 'WORLD', -12000, 9000, -12000, 9000, use_stats='univar', bin_methd='man', writeStats=True, makeBoxPlot=True, makeHistogram=True)

# close stats file
fileOut.close()


#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
# Ultra Proeminent peaks
# correlation of real elevation vs. gdems
# 
# import ultras.shp
# (originally a kml from http://www.peaklist.org/ultras.html)

# v.in.gdal

# sample gdems rasters with vector
ultras = 'ultras' #vector

ow = True
grass.run_command('r.mask', flags='r', overwrite=ow)

# list of GDEMS 
# don't use all gdems, only those with 30'' and 01' resolution
gdems_list = ['gdem_srtm30_plus_v9', 'gdem_globe', 'gdem_gtopo30', \
'gdem_gmted_land', 'gdem_ace2_int', 'gdem_etopo1_ice',]

# files for results
workDir = '/Volumes/Macintosh HD2T/Dropbox/artigos/global_DEMs/stats/stats_plots_ULTRAS'
os.chdir(workDir)

# stats and plots
for gdem in gdems_list[::-1]: 
    grass.run_command('g.region', rast=gdem, flags='pa')
    col = gdem.split('_')[1]
    grass.run_command('v.db.addcol', map='ultras', columns=col+' integer') # create columns in vector file
    grass.run_command('v.what.rast', vector='ultras', raster=gdem, column=col) # sample

# check file
# grass.vector_db_select('ultras')['columns']
# ['cat', 'Name', 'long', 'lat', 'alt', 'proem', 'Country', 'ID', 'srtm30', 'etopo1', 'ace2', 'gmted', 'gtopo30', 'globe']

# get attribute data
attrs = grass.vector_db_select('ultras', columns = 'Name,alt,proem,Country,ID,srtm30,globe,gtopo30,gmted,ace2,etopo1')

# lists for vals
cat = []
name = []
alt = []
proem = []
country = []
pid = []
srtm30 = []
globe = []
gtopo30 = []
gmted = []
ace2 = []
etopo1 = []

# every attr in a list
for key in attrs['values']:
    # print key, attrs['values'][key][0]
    # print 
    cat.append(key)
    name.append(attrs['values'][key][0])
    alt.append(int(attrs['values'][key][1]))
    proem.append(int(attrs['values'][key][2]))
    country.append(attrs['values'][key][3])
    pid.append(attrs['values'][key][4])
    srtm30.append(int(attrs['values'][key][5]))
    globe.append(int(attrs['values'][key][6]))
    gtopo30.append(int(attrs['values'][key][7]))
    gmted.append(int(attrs['values'][key][8]))
    ace2.append(int(attrs['values'][key][9]))
    etopo1.append(int(attrs['values'][key][10]))

gdems = [srtm30,globe,gtopo30,gmted,ace2,etopo1]
gdems_str = ['srtm30','globe','gtopo30','gmted','ace2','etopo1']

for dem in gdems:
    plt.plot(alt,dem, 'o')
    plt.plot([-2000,9000],[-2000,9000]) # 45deg line
    plt.ylabel(dem)
    plt.xlabel('altitude (m)')
    plt.xlim(-2000,9000)
    plt.ylim(-2000,9000)
    # plt.show()  
    figpeak = 'ultras_' + gdems_str[gdems.index(dem)] + '_.svg'
    plt.savefig(figpeak)
    plt.clf()


# compare differences between altitude and gdems
def compDiff(array1, array2, gdem):
    ar1=np.array(array1)
    ar2=np.array(array2)
    dif = np.abs(ar1-ar2)
    # 
    mnd = min(dif)
    mxd = max(dif)
    avg = np.mean(dif)
    med = np.median(dif)
    std = np.std(dif)
    var = np.var(dif)
    p25 = np.percentile(dif,25)
    p75 = np.percentile(dif,75)
    # 
    print 'diff alt - ', gdem
    print 'min: max: mean: median: stddev: variance: p25: p75:'
    print mnd, mxd, avg, med, std, var, p25, p75
    # 
    IQR = p75 - p25
    bin_width = 2. * IQR * len(alt)**(-1./3.)
    bin_width = round5(bin_width)
    print 'bin width: ', bin_width
    nbins = round5((mxd - mnd)/bin_width)
    # 
    plt.hist(dif, bins=nbins)
    fig = 'hist_diff_' + gdem + '.svg'
    plt.savefig(fig)
    plt.clf()
    # "mode"
    h = np.histogram(dif, bins=nbins)
    mx = np.where(h[0]==max(h[0]))[0][0]
    print 'mode: ', h[1][mx]






