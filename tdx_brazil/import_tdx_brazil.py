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


# This file has functions used to import the data (DEMs) used in the paper


import sys, os, csv
import grass.script as grass
import itertools
import numpy as np

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------

# helper funcs
def import_tiles(dem_list, workDir):
    ''' import a list of DEM tiles '''
    for tile in dem_list:
        inName = workDir + tile[0]
        outName = tile[1]
        print inName, outName
        grass.run_command('r.in.gdal', input=inName, output=outName, overwrite=True)
        # grass.run_command('r.info', map=outName)

def clean_dems(dem_list):
    ''' clean up imported DEMs. delete values lower than zero, rename files '''
    for dem in dem_list:
        grass.run_command('g.region', raster=dem, flags='pa')
        clean = dem + '_clean'
        grass.mapcalc("${out} = if(${orig} <= 0, null(), ${orig})",
            out=clean,
            orig=dem,
            overwrite = True)
        grass.run_command('g.remove', flags='f', type='raster', name=dem)
        grass.run_command('g.rename', raster=(clean,dem))

def clean_coastal_dems(dem_list):
    ''' clean up coastal DEMs. apply mask, delete values lower than zero, rename files '''
    for dem in dem_list:
        dem_area = dem.split('_')[-1]
        poly = dem_area + '_poly'
        grass.run_command('g.region', raster=dem, flags='pa')
        grass.run_command('r.mask', vector=poly, overwrite=True)
        clean = dem + '_clean'
        grass.mapcalc("${out} = if(${orig} <= 0, null(), ${orig})",
            out=clean,
            orig=dem,
            overwrite = True)
        grass.run_command('g.remove', flags='f', type='raster', name=dem)
        grass.run_command('g.rename', raster=(clean,dem))
        grass.run_command('r.mask', flags='r')

def clean_dems_wan(dem_list,thresh):
    ''' clean up imported DEMs. delete values lower than zero, rename files '''
    for dem in dem_list:
        grass.run_command('g.region', raster=dem, flags='pa')
        clean = dem + '_clean'
        tdx_wam = dem.split('_')[0] + '_wam_' + dem.split('_')[1]
        grass.mapcalc("${out} = if(${wam} > ${val}, null(), ${orig})",
            out=clean,
            orig=dem,
            wam=tdx_wam,
            val=thresh,
            overwrite = True)
        grass.run_command('g.remove', flags='f', type='raster', name=dem)
        grass.run_command('g.rename', raster=(clean,dem))

def make_shade(dem_list):
    '''make shaded reliefs for visual inspection'''
    for dem in dem_list:
        grass.run_command('g.region', raster=dem, flags='pa')
        azim = 315
        inc = 20
        out_shade = dem + '_shade_' + str(azim) + '_' + str(inc)
        grass.run_command('r.relief', input=dem, output=out_shade, altitude=inc, azimuth=azim, overwrite=True)


#----------------------------------------------------
#----------------------------------------------------
#  Import TanDEM-X data
#----------------------------------------------------
#----------------------------------------------------

# set working directory
workDir = '/Volumes/MacintoshHD2/geodata/Global_DEMS/TanDEM-X/'
os.chdir(workDir)

# name the files - DEMs
Araca12 = 'GEOL_0538/TDM1_DEM__04_N00W064_V01_C_Araca/DEM/TDM1_DEM__04_N00W064_DEM.tif'
Barcelos12 = 'GEOL_0538/TDM1_DEM__04_S01W063_V01_C_Barcelos/DEM/TDM1_DEM__04_S01W063_DEM.tif'
Pantanal12 = 'GEOL_0538/TDM1_DEM__04_S19W057_V01_C_Pantanal/DEM/TDM1_DEM__04_S19W057_DEM.tif'
Itatiaia12 = 'GEOL_0538/TDM1_DEM__04_S23W045_V01_C_Itatiaia/DEM/TDM1_DEM__04_S23W045_DEM.tif'
RioClaro12 = 'GEOL_0538/TDM1_DEM__04_S23W048_V01_C_RioClaro/DEM/TDM1_DEM__04_S23W048_DEM.tif'
Paraty12 = 'GEOL_0538/TDM1_DEM__04_S24W045_V01_C_Paraty/DEM/TDM1_DEM__04_S24W045_DEM.tif'
RCSB12 = 'GEOL_0538/TDM1_DEM__04_S24W046_V01_C_RCSB/DEM/TDM1_DEM__04_S24W046_DEM.tif'
Iporanga12 = 'GEOL_0538/TDM1_DEM__04_S25W049_V01_C_Iporanga/DEM/TDM1_DEM__04_S25W049_DEM.tif'
Floripa12 = 'GEOL_0538/TDM1_DEM__04_S28W049_V01_C_Floripa/DEM/TDM1_DEM__04_S28W049_DEM.tif'
Garopaba12 = 'GEOL_0538/TDM1_DEM__04_S29W049_V01_C_Garopaba/DEM/TDM1_DEM__04_S29W049_DEM.tif'

Araca30 = 'GEOL_0538/TDM1_DEM__10_N00W064_V01_C_Araca/DEM/TDM1_DEM__10_N00W064_DEM.tif'
Barcelos30 = 'GEOL_0538/TDM1_DEM__10_S01W063_V01_C_Barcelos/DEM/TDM1_DEM__10_S01W063_DEM.tif'
Pantanal30 = 'GEOL_0538/TDM1_DEM__10_S19W057_V01_C_Pantanal/DEM/TDM1_DEM__10_S19W057_DEM.tif'
Itatiaia30 = 'GEOL_0538/TDM1_DEM__10_S23W045_V01_C_Itatiaia/DEM/TDM1_DEM__10_S23W045_DEM.tif'
RioClaro30 = 'GEOL_0538/TDM1_DEM__10_S23W048_V01_C_RioClaro/DEM/TDM1_DEM__10_S23W048_DEM.tif'
Paraty30 = 'GEOL_0538/TDM1_DEM__10_S24W045_V01_C_Paraty/DEM/TDM1_DEM__10_S24W045_DEM.tif'
RCSB30 = 'GEOL_0538/TDM1_DEM__10_S24W046_V01_C_RCSB/DEM/TDM1_DEM__10_S24W046_DEM.tif'
Iporanga30 = 'GEOL_0538/TDM1_DEM__10_S25W049_V01_C_Iporanga/DEM/TDM1_DEM__10_S25W049_DEM.tif'
Floripa30 = 'GEOL_0538/TDM1_DEM__10_S28W049_V01_C_Floripa/DEM/TDM1_DEM__10_S28W049_DEM.tif'
Garopaba30 = 'GEOL_0538/TDM1_DEM__10_S29W049_V01_C_Garopaba/DEM/TDM1_DEM__10_S29W049_DEM.tif'


tdx_tiles = [[Araca12,      'tdx12_Araca'], \
            [Barcelos12,    'tdx12_Barcelos'], \
            [Pantanal12,    'tdx12_Pantanal'], \
            [Itatiaia12,    'tdx12_Itatiaia'], \
            [RioClaro12,    'tdx12_RioClaro'], \
            [Paraty12,      'tdx12_Paraty'], \
            [RCSB12,        'tdx12_RCSB'], \
            [Iporanga12,    'tdx12_Iporanga'], \
            [Floripa12,     'tdx12_Floripa'], \
            [Garopaba12,    'tdx12_Garopaba'], \
            [Araca30,       'tdx30_Araca'], \
            [Barcelos30,    'tdx30_Barcelos'], \
            [Pantanal30,    'tdx30_Pantanal'], \
            [Itatiaia30,    'tdx30_Itatiaia'], \
            [RioClaro30,    'tdx30_RioClaro'], \
            [Paraty30,      'tdx30_Paraty'], \
            [RCSB30,        'tdx30_RCSB'], \
            [Iporanga30,    'tdx30_Iporanga'], \
            [Floripa30,     'tdx30_Floripa'], \
            [Garopaba30,    'tdx30_Garopaba']]


# import files - DEMs
import_tiles(tdx_tiles,workDir)
tdx = [tile[1] for tile in tdx_tiles]
# clean up - delete values lower than zero, rename files
clean_dems(tdx)
# make shaded reliefs for visual inspection
make_shade(tdx)

tdx = ['tdx12_Araca', 'tdx12_Barcelos', 'tdx12_Pantanal', 'tdx12_Itatiaia', 'tdx12_RioClaro', \
'tdx12_Paraty', 'tdx12_RCSB', 'tdx12_Iporanga', 'tdx12_Floripa', 'tdx12_Garopaba', 'tdx30_Araca', \
'tdx30_Barcelos', 'tdx30_Pantanal', 'tdx30_Itatiaia', 'tdx30_RioClaro', 'tdx30_Paraty', 'tdx30_RCSB', \
'tdx30_Iporanga', 'tdx30_Floripa', 'tdx30_Garopaba']



#----------------------------------------------------
#----------------------------------------------------
#  Water mask - WAM
# __________________________________________________________
# set working directory
workDir = '/Volumes/MacintoshHD2/geodata/Global_DEMS/TanDEM-X/'
os.chdir(workDir)

# name the files - WAMs
Araca12 = 'GEOL_0538/TDM1_DEM__04_N00W064_V01_C_Araca/AUXFILES/TDM1_DEM__04_N00W064_WAM.tif'
Barcelos12 = 'GEOL_0538/TDM1_DEM__04_S01W063_V01_C_Barcelos/AUXFILES/TDM1_DEM__04_S01W063_WAM.tif'
Pantanal12 = 'GEOL_0538/TDM1_DEM__04_S19W057_V01_C_Pantanal/AUXFILES/TDM1_DEM__04_S19W057_WAM.tif'
Itatiaia12 = 'GEOL_0538/TDM1_DEM__04_S23W045_V01_C_Itatiaia/AUXFILES/TDM1_DEM__04_S23W045_WAM.tif'
RioClaro12 = 'GEOL_0538/TDM1_DEM__04_S23W048_V01_C_RioClaro/AUXFILES/TDM1_DEM__04_S23W048_WAM.tif'
Paraty12 = 'GEOL_0538/TDM1_DEM__04_S24W045_V01_C_Paraty/AUXFILES/TDM1_DEM__04_S24W045_WAM.tif'
RCSB12 = 'GEOL_0538/TDM1_DEM__04_S24W046_V01_C_RCSB/AUXFILES/TDM1_DEM__04_S24W046_WAM.tif'
Iporanga12 = 'GEOL_0538/TDM1_DEM__04_S25W049_V01_C_Iporanga/AUXFILES/TDM1_DEM__04_S25W049_WAM.tif'
Floripa12 = 'GEOL_0538/TDM1_DEM__04_S28W049_V01_C_Floripa/AUXFILES/TDM1_DEM__04_S28W049_WAM.tif'
Garopaba12 = 'GEOL_0538/TDM1_DEM__04_S29W049_V01_C_Garopaba/AUXFILES/TDM1_DEM__04_S29W049_WAM.tif'

Araca30 = 'GEOL_0538/TDM1_DEM__10_N00W064_V01_C_Araca/AUXFILES/TDM1_DEM__10_N00W064_WAM.tif'
Barcelos30 = 'GEOL_0538/TDM1_DEM__10_S01W063_V01_C_Barcelos/AUXFILES/TDM1_DEM__10_S01W063_WAM.tif'
Pantanal30 = 'GEOL_0538/TDM1_DEM__10_S19W057_V01_C_Pantanal/AUXFILES/TDM1_DEM__10_S19W057_WAM.tif'
Itatiaia30 = 'GEOL_0538/TDM1_DEM__10_S23W045_V01_C_Itatiaia/AUXFILES/TDM1_DEM__10_S23W045_WAM.tif'
RioClaro30 = 'GEOL_0538/TDM1_DEM__10_S23W048_V01_C_RioClaro/AUXFILES/TDM1_DEM__10_S23W048_WAM.tif'
Paraty30 = 'GEOL_0538/TDM1_DEM__10_S24W045_V01_C_Paraty/AUXFILES/TDM1_DEM__10_S24W045_WAM.tif'
RCSB30 = 'GEOL_0538/TDM1_DEM__10_S24W046_V01_C_RCSB/AUXFILES/TDM1_DEM__10_S24W046_WAM.tif'
Iporanga30 = 'GEOL_0538/TDM1_DEM__10_S25W049_V01_C_Iporanga/AUXFILES/TDM1_DEM__10_S25W049_WAM.tif'
Floripa30 = 'GEOL_0538/TDM1_DEM__10_S28W049_V01_C_Floripa/AUXFILES/TDM1_DEM__10_S28W049_WAM.tif'
Garopaba30 = 'GEOL_0538/TDM1_DEM__10_S29W049_V01_C_Garopaba/AUXFILES/TDM1_DEM__10_S29W049_WAM.tif'


tdx_wam_tiles = [[Araca12,  'tdx12_wam_Araca'], \
            [Barcelos12,    'tdx12_wam_Barcelos'], \
            [Pantanal12,    'tdx12_wam_Pantanal'], \
            [Itatiaia12,    'tdx12_wam_Itatiaia'], \
            [RioClaro12,    'tdx12_wam_RioClaro'], \
            [Paraty12,      'tdx12_wam_Paraty'], \
            [RCSB12,        'tdx12_wam_RCSB'], \
            [Iporanga12,    'tdx12_wam_Iporanga'], \
            [Floripa12,     'tdx12_wam_Floripa'], \
            [Garopaba12,    'tdx12_wam_Garopaba'], \
            [Araca30,       'tdx30_wam_Araca'], \
            [Barcelos30,    'tdx30_wam_Barcelos'], \
            [Pantanal30,    'tdx30_wam_Pantanal'], \
            [Itatiaia30,    'tdx30_wam_Itatiaia'], \
            [RioClaro30,    'tdx30_wam_RioClaro'], \
            [Paraty30,      'tdx30_wam_Paraty'], \
            [RCSB30,        'tdx30_wam_RCSB'], \
            [Iporanga30,    'tdx30_wam_Iporanga'], \
            [Floripa30,     'tdx30_wam_Floripa'], \
            [Garopaba30,    'tdx30_wam_Garopaba']]


# import files - WAMs
import_tiles(tdx_wam_tiles,workDir)
tdx_wam = [tile[1] for tile in tdx_wam_tiles]


tdx_wam = ['tdx12_wam_Araca', 'tdx12_wam_Barcelos', 'tdx12_wam_Pantanal', 'tdx12_wam_Itatiaia', 'tdx12_wam_RioClaro', \
'tdx12_wam_Paraty', 'tdx12_wam_RCSB', 'tdx12_wam_Iporanga', 'tdx12_wam_Floripa', 'tdx12_wam_Garopaba', 'tdx30_wam_Araca', \
'tdx30_wam_Barcelos', 'tdx30_wam_Pantanal', 'tdx30_wam_Itatiaia', 'tdx30_wam_RioClaro', 'tdx30_wam_Paraty', 'tdx30_wam_RCSB', \
'tdx30_wam_Iporanga', 'tdx30_wam_Floripa', 'tdx30_wam_Garopaba']

# clean up - delete values lower than zero, rename files
clean_dems_wan(tdx,65)
# make shaded reliefs for visual inspection
# make_shade(tdx)


# set categories and colors (reclassified maps)
rule = '''0 thru 2 = NULL\n
3 thru 4 = 1 1 x AMP_thresh_1\n
5 thru 6 = 2 2 x AMP_thresh_1\n
7 thru 8 = 3 3 x AMP_thresh_1\n
9 thru 16 = 4 1 x AMP_thresh_2\n
17 thru 25 = 5 2 x AMP_thresh_2\n
26 thru 32 = 6 3 x AMP_thresh_2\n
32 thru 64 = 7 1 x COH_thresh_1\n
65 thru 96 = 8 2 x COH_thresh_1\n
97 thru 127 = 9 3 x COH_thresh_1\n
'''
for tdx in tdx_wam:
    grass.write_command('r.reclass', input=tdx, output=tdx+'_cats', rules='-', stdin=rule, overwrite=True)


#----------------------------------------------------
#----------------------------------------------------
#  Import SRTM 01 sec data
#----------------------------------------------------
#----------------------------------------------------

# set working directory
workDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/gdems/srtm_30m/'
os.chdir(workDir)

# name the files
Araca30 =       'N00W064.hgt'
Barcelos30 =    'S01W063.hgt'
Pantanal30 =    'S19W057.hgt'
Itatiaia30 =    'S23W045.hgt'
RioClaro30 =    'S23W048.hgt'
Paraty30 =      'S24W045.hgt'
RCSB30 =        'S24W046.hgt'
Iporanga30 =    'S25W049.hgt'
Floripa30 =     'S28W049.hgt'
Garopaba30 =    'S29W049.hgt'


srtm_tiles = [[Araca30,     'srtm30_Araca'], \
            [Barcelos30,    'srtm30_Barcelos'], \
            [Pantanal30,    'srtm30_Pantanal'], \
            [Itatiaia30,    'srtm30_Itatiaia'], \
            [RioClaro30,    'srtm30_RioClaro'], \
            [Paraty30,      'srtm30_Paraty'], \
            [RCSB30,        'srtm30_RCSB'], \
            [Iporanga30,    'srtm30_Iporanga'], \
            [Floripa30,     'srtm30_Floripa'], \
            [Garopaba30,    'srtm30_Garopaba']]

# import files
import_tiles(srtm_tiles,workDir)
# get list of raster tiles
srtm = [tile[1] for tile in srtm_tiles]
# clean up - delete values lower than zero, rename files
clean_dems(srtm)
# make shaded reliefs for visual inspection
# make_shade(srtm)

srtm = ['srtm30_Araca', 'srtm30_Barcelos', 'srtm30_Pantanal', 'srtm30_Itatiaia', 'srtm30_RioClaro', \
'srtm30_Paraty', 'srtm30_RCSB', 'srtm30_Iporanga', 'srtm30_Floripa', 'srtm30_Garopaba']


#----------------------------------------------------
#----------------------------------------------------
#  Import ASTER GDEM data
#----------------------------------------------------
#----------------------------------------------------

# set working directory
workDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/gdems/aster_gdem_v2/'
os.chdir(workDir)

# name the files
Araca30 =       'ASTGTM2_N00W064/ASTGTM2_N00W064_dem.tif'
Barcelos30 =    'ASTGTM2_S01W063/ASTGTM2_S01W063_dem.tif'
Pantanal30 =    'ASTGTM2_S19W057/ASTGTM2_S19W057_dem.tif'
Itatiaia30 =    'ASTGTM2_S23W045/ASTGTM2_S23W045_dem.tif'
RioClaro30 =    'ASTGTM2_S23W048/ASTGTM2_S23W048_dem.tif'
Paraty30 =      'ASTGTM2_S24W045/ASTGTM2_S24W045_dem.tif'
RCSB30 =        'ASTGTM2_S24W046/ASTGTM2_S24W046_dem.tif'
Iporanga30 =    'ASTGTM2_S25W049/ASTGTM2_S25W049_dem.tif'
Floripa30 =     'ASTGTM2_S28W049/ASTGTM2_S28W049_dem.tif'
Garopaba30 =    'ASTGTM2_S29W049/ASTGTM2_S29W049_dem.tif'


aster_tiles = [[Araca30,    'aster30_Araca'], \
            [Barcelos30,    'aster30_Barcelos'], \
            [Pantanal30,    'aster30_Pantanal'], \
            [Itatiaia30,    'aster30_Itatiaia'], \
            [RioClaro30,    'aster30_RioClaro'], \
            [Paraty30,      'aster30_Paraty'], \
            [RCSB30,        'aster30_RCSB'], \
            [Iporanga30,    'aster30_Iporanga'], \
            [Floripa30,     'aster30_Floripa'], \
            [Garopaba30,    'aster30_Garopaba']]

# import files
import_tiles(aster_tiles,workDir)
# get list of raster tiles
aster = [tile[1] for tile in aster_tiles]
# clean up - delete values lower than zero, rename files
clean_dems(aster)
# make shaded reliefs for visual inspection
# make_shade(aster)

aster = ['aster30_Araca', 'aster30_Barcelos', 'aster30_Pantanal', 'aster30_Itatiaia', 'aster30_RioClaro', \
'aster30_Paraty', 'aster30_RCSB', 'aster30_Iporanga', 'aster30_Floripa', 'aster30_Garopaba']


#----------------------------------------------------
#----------------------------------------------------
#  Import ALOS W3D30 data - AVERAGE
#----------------------------------------------------
#----------------------------------------------------

# set working directory
workDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/gdems/alos_aw3d30/'
os.chdir(workDir)

# name the files
Araca30 =       'N000W064/AVERAGE/N000W064_AVE_DSM.tif'
Barcelos30 =    'S001W063/AVERAGE/S001W063_AVE_DSM.tif'
Pantanal30 =    'S019W057/AVERAGE/S019W057_AVE_DSM.tif'
Itatiaia30 =    'S023W045/AVERAGE/S023W045_AVE_DSM.tif'
RioClaro30 =    'S023W048/AVERAGE/S023W048_AVE_DSM.tif'
Paraty30 =      'S024W045/AVERAGE/S024W045_AVE_DSM.tif'
RCSB30 =        'S024W046/AVERAGE/S024W046_AVE_DSM.tif'
Iporanga30 =    'S025W049/AVERAGE/S025W049_AVE_DSM.tif'
Floripa30 =     'S028W049/AVERAGE/S028W049_AVE_DSM.tif'
Garopaba30 =    'S029W049/AVERAGE/S029W049_AVE_DSM.tif'


alos_tiles_avg = [[Araca30,     'aw3d30_avg_Araca'], \
            [Barcelos30,    'aw3d30_avg_Barcelos'], \
            [Pantanal30,    'aw3d30_avg_Pantanal'], \
            [Itatiaia30,    'aw3d30_avg_Itatiaia'], \
            [RioClaro30,    'aw3d30_avg_RioClaro'], \
            [Paraty30,      'aw3d30_avg_Paraty'], \
            [RCSB30,        'aw3d30_avg_RCSB'], \
            [Iporanga30,    'aw3d30_avg_Iporanga'], \
            [Floripa30,     'aw3d30_avg_Floripa'], \
            [Garopaba30,    'aw3d30_avg_Garopaba']]

# import files
import_tiles(alos_tiles_avg,workDir)
# get list of raster tiles
alos_avg = [tile[1] for tile in alos_tiles_avg]
# clean up - delete values lower than zero, rename files
clean_dems(alos_avg)
# make shaded reliefs for visual inspection
# make_shade(alos_avg)

alos_avg = ['aw3d30_avg_Araca', 'aw3d30_avg_Barcelos', 'aw3d30_avg_Pantanal', \
'aw3d30_avg_Itatiaia', 'aw3d30_avg_RioClaro', 'aw3d30_avg_Paraty', 'aw3d30_avg_RCSB', \
'aw3d30_avg_Iporanga', 'aw3d30_avg_Floripa', 'aw3d30_avg_Garopaba']


#----------------------------------------------------
#----------------------------------------------------
#  Patch areas that have 2 tiles
#----------------------------------------------------
#----------------------------------------------------

# Garopaba
# TDX 12m
grass.run_command('g.region', n=-27, s=-29, w=-49, e=-48, res='0:0:00.4', flags='pa')
grass.run_command('r.patch', input=['tdx12_Floripa','tdx12_Garopaba'], output='tdx12_SantaCatarina', overwrite=True)
# TDX 30m
grass.run_command('g.region', n=-27, s=-29, w=-49, e=-48, res='0:0:01', flags='pa')
grass.run_command('r.patch', input=['tdx30_Floripa','tdx30_Garopaba'], output='tdx30_SantaCatarina', overwrite=True)
# SRTM30m
grass.run_command('r.patch', input=['srtm30_Floripa','srtm30_Garopaba'], output='srtm30_SantaCatarina', overwrite=True)
# ASTER 30m
grass.run_command('r.patch', input=['aster30_Floripa','aster30_Garopaba'], output='aster30_SantaCatarina', overwrite=True)
# ALOS 30m - AVG
grass.run_command('r.patch', input=['aw3d30_avg_Floripa','aw3d30_avg_Garopaba'], output='aw3d30_avg_SantaCatarina', overwrite=True)

# Itatiaia - Serra do Mar
# TDX 12m
grass.run_command('g.region', n=-22, s=-23.5, w=-45, e=-44, res='0:0:00.4', flags='pa')
grass.run_command('r.patch', input=['tdx12_Itatiaia','tdx12_Paraty'], output='tdx12_SerraMar', overwrite=True)
# TDX 30m
grass.run_command('g.region', n=-22, s=-23.5, w=-45, e=-44, res='0:0:01', flags='pa')
grass.run_command('r.patch', input=['tdx30_Itatiaia','tdx30_Paraty'], output='tdx30_SerraMar', overwrite=True)
# SRTM30m
grass.run_command('r.patch', input=['srtm30_Itatiaia','srtm30_Paraty'], output='srtm30_SerraMar', overwrite=True)
# ASTER 30m
grass.run_command('r.patch', input=['aster30_Itatiaia','aster30_Paraty'], output='aster30_SerraMar', overwrite=True)
# ALOS 30m - AVG
grass.run_command('r.patch', input=['aw3d30_avg_Itatiaia','aw3d30_avg_Paraty'], output='aw3d30_avg_SerraMar', overwrite=True)

# new list of patched tiles
patched = ['tdx12_SantaCatarina','tdx30_SantaCatarina','srtm30_SantaCatarina',\
'aster30_SantaCatarina','aw3d30_SantaCatarina', 'tdx12_SerraMar', \
'tdx30_SerraMar', 'srtm30_SerraMar', 'aster30_SerraMar','aw3d30_SerraMar']

# make shades of patched tiles
# make_shade(patched)


#----------------------------------------------------
#----------------------------------------------------
#  Create vector mask of coastal areas 
#----------------------------------------------------
#----------------------------------------------------

# Santa Catarina
# SRTM30m
dem = 'srtm30_SantaCatarina'
poly = 'SantaCatarina_poly'
clean = poly + '_clean'
grass.run_command('g.region', n=-27, s=-29, w=-49, e=-48, res='0:0:03', flags='pa')
grass.run_command('r.mask', raster=dem, maskcats='1 thru 1500', overwrite=True)
grass.run_command('r.to.vect', input='MASK', output=poly, type='area', overwrite=True)
grass.run_command('v.clean', input=poly, output=clean, type='area', tool='rmarea', thres=1000000, overwrite=True)
grass.run_command('g.remove', flags='f', type='vector', name=poly)
grass.run_command('g.rename', vector=(clean,poly))

# Itatiaia - Serra do Mar
# SRTM30m
dem = 'srtm30_SerraMar'
poly = 'SerraMar_poly'
clean = poly + '_clean'
grass.run_command('g.region', n=-22, s=-23.5, w=-45, e=-44, res='0:0:03', flags='pa')
grass.run_command('r.mask', raster=dem, maskcats='1 thru 3000', overwrite=True)
grass.run_command('r.to.vect', input='MASK', output=poly, type='area', overwrite=True)
grass.run_command('v.clean', input=poly, output=clean, type='area', tool='rmarea', thres=10000000, overwrite=True)
grass.run_command('g.remove', flags='f', type='vector', name=poly)
grass.run_command('g.rename', vector=(clean,poly))


#----------------------------------------------------
#----------------------------------------------------
#  Clean DEMs of coastal areas
#----------------------------------------------------
#----------------------------------------------------
coast = ['tdx12_SantaCatarina','tdx30_SantaCatarina','tdx12_SerraMar','tdx30_SerraMar','srtm30_SantaCatarina','srtm30_SerraMar','aster30_SantaCatarina'\
,'aster30_SerraMar','aw3d30_SantaCatarina','aw3d30_SerraMar']

clean_coastal_dems(coast)
# make_shade(coast)



#----------------------------------------------------
#----------------------------------------------------
#  Fix elevation of DEMs that are referenced to EGM96:
#  SRTM, ASTER GDEM, ALOS AW3D30
#----------------------------------------------------
#----------------------------------------------------

# 15-minute geoid height grid 
# http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
# C implementation of f477 (to calculate point geoid undulations at any given WGS 84 latitude and longitude by spherical harmonic synthesis):
# https://sourceforge.net/projects/egm96-f477-c/
# compile with: gcc -o f477 f477.c

# create text file of lat/long for south america (interval 0.1 degrees):
workdir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/gdems/egm96/egm96_c/'
os.chdir(workdir)

# make a XY file with points at lat/lonmg every 0.1 degrees
def make_grid_ll(north,south,east,west,interval):
    '''make X,Y values for a grid'''
    lx = []
    ly = []
    longs = np.arange(west, east-interval, interval)
    lats = np.arange(south, north+interval, interval)
    for y,x in itertools.product(lats,longs):
        ly.append(y)
        lx.append(x)
    return lx, ly
# south america
lx, ly = make_grid_ll(north=14, south=-57, east=-30, west=-85, interval=0.1)

# export grid to file
ll4pd = [('lat', ly),('long', lx)]
ll_pd = pd.DataFrame.from_items(ll4pd)
ll_pd.to_csv('input.dat', sep=' ',columns=('lat','long'), header=False, index=False)

# calculate geoid height for south america
# (run f477 manually)

# f477 results in a file that GRASS can't read. Open from csv with pandas and export a new one
egm = pd.read_csv('OUTF477.DAT', sep='\s+', header=None, names=('lat','long','und'))
egm.to_csv('egm96.xyz', sep=' ',columns=('lat','long','und'), header=False, index=False)
# import into GRASS as vector points
grass.run_command('v.in.ascii', input=workdir+'egm96.xyz', output='egm96', separator=' ', x=2, y=1, z=3, flags='ztb', overwrite=True)


# lists of dems
srtm = ['srtm30_Araca', 'srtm30_Barcelos', 'srtm30_Pantanal', 'srtm30_SerraMar', 'srtm30_RioClaro', \
'srtm30_Iporanga', 'srtm30_SantaCatarina']

aster = ['aster30_Araca', 'aster30_Barcelos', 'aster30_Pantanal', 'aster30_SerraMar', 'aster30_RioClaro', \
'aster30_Iporanga', 'aster30_SantaCatarina']

alos = ['aw3d30_Araca', 'aw3d30_Barcelos', 'aw3d30_Pantanal', 'aw3d30_SerraMar', \
'aw3d30_RioClaro', 'aw3d30_Iporanga', 'aw3d30_SantaCatarina']

# for each area, interpolate egm96 with bilinear interpolation
for dem in aster:
    area_name = dem.split('_')[1]
    grass.run_command('g.region', raster=dem, flags='a')
    grass.run_command('v.surf.bspline', input='egm96', raster_output='egm96_30m_'+area_name, method='bilinear', overwrite=True)

# add egm96 to SRTM, ASTER, ALOS
for dem in itertools.chain.from_iterable([srtm, aster, alos]):
    grass.run_command('g.region', raster=dem, flags='pa')
    area_name = dem.split('_')[-1]
    grass.mapcalc("${out} = ${rast1} + ${rast2}", \
                out = dem.split('_')[0]+'_wgs84_'+dem.split('_')[-1],
                rast1 = dem,
                rast2 = 'egm96_30m_'+area_name, overwrite=True)





#----------------------------------------------------
#----------------------------------------------------
#  Final lists of tiles (WGS84)
#----------------------------------------------------
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
#----------------------------------------------------
#  make screenshots of shaded relief maps
#----------------------------------------------------
#----------------------------------------------------
# set working directory
workDir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/stats_plots/'
os.chdir(workDir)

# make shaded relief maps for all tiles
for area_list in areas:
    make_shade(area_list)

# fix  colors
# set colortable of TDX12 to 'haxby'
# for dem in itertools.chain.from_iterable(areas):
tdx12_list = ['tdx12_Araca', 'tdx12_Barcelos', 'tdx12_Pantanal', 'tdx12_SerraMar', 'tdx12_RioClaro', \
'tdx12_Iporanga', 'tdx12_SantaCatarina']

for dem in tdx12_list:
    grass.run_command('g.region', raster=dem, flags='pa')
    grass.run_command('r.colors', map=dem, color='haxby', flags='e')

# copy tdx12 colortable
for area in areas:
    area_name = area[0].split('_')[-1]
    for dem in area: 
        tdx12 = 'tdx12_' + area_name
        if not dem.startswith('tdx12'):
            grass.run_command('g.region', raster=dem, flags='pa')
            grass.run_command('r.colors', map=dem, raster=tdx12)


# dict with range for colorscales
range_dict = {'Araca':(0,1500),'Barcelos':(0,100),'Pantanal':(80,180),'Iporanga':(0,1250),\
            'RioClaro':(400,1000),'SantaCatarina':(0,1200),'SerraMar':(0,2600)}

# export display as pngs
for area in areas:
    area_name = area[0].split('_')[-1]
    grass.run_command('g.region', raster=area[0], flags='pa', res='0:0:01')
    for dem in area:
        out_shade = dem + '_shade_315_20'
        grass.run_command('d.mon', start='cairo', output=dem+'.png', resolution = '2', \
            height=600, width=600, overwrite=True)
        grass.run_command('d.shade', shade=out_shade, color=dem)
        grass.run_command('d.grid', size='0.25', text_color='black', color='black', \
            border_color='black', fontsize=18, flags='c')
        grass.run_command('d.legend', raster=dem, flags='fsb', at='4,25,2,4', font='Helvetica', \
            fontsize=13, bgcolor='240:240:240', range=range_dict[area_name], title='(m)')
        grass.run_command('d.mon', stop='cairo')
        print('PNG export OK. map: ' + dem)



#----------------------------------------------------
#----------------------------------------------------
#  Import GeoPNGs exported manually from QGIS with Google Sat Layer
#  (had to fix some images in GIMP)
#----------------------------------------------------
#----------------------------------------------------

# set working directory
workdir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/pngs_gearth/'
os.chdir(workdir)

# PNG files
flist = os.listdir(workdir)
flist = [f for f in flist if f.endswith('png')]

# gdalwarp
from subprocess import call
for img in flist:
    img_out = img[:-4]+'.tif'
    call(['/usr/local/bin/gdalwarp', '-s_srs', 'EPSG:3857', '-t_srs', 'EPSG:4326', img, img_out])

# TIFF files 
flist = os.listdir(workdir)
flist = [f for f in flist if f.endswith('tif')]

# import TIFFs
for img in flist:
    img_out = img[:-4]
    grass.run_command('r.in.gdal', input=workdir+img, output=img_out, overwrite=True)

# make RGB GRASS rasters
for img in flist:
    img_imp = img[:-4]
    r = img_imp+'.red'
    g = img_imp+'.green'
    b = img_imp+'.blue'
    grass.run_command('g.region', raster=r, flags='pa')
    grass.run_command('r.composite', red=r, green=g, blue=b, levels=32, output=img_imp, overwrite=True)

# patch maps
# araca
patchlist = ['araca4', 'araca3', 'araca2', 'araca1']
grass.run_command('g.region', raster='srtm30_Araca', flags='pa', res='0:00:01')
grass.run_command('r.patch', input=patchlist, output='araca_ge', overwrite=True)

# barcelos
patchlist = ['barcelos4', 'barcelos3', 'barcelos2', 'barcelos1']
grass.run_command('g.region', raster='srtm30_Barcelos', flags='pa', res='0:00:01')
grass.run_command('r.patch', input=patchlist, output='barcelos_ge', overwrite=True)

# pantanal
patchlist = ['pantanal4', 'pantanal3', 'pantanal2', 'pantanal1']
grass.run_command('g.region', raster='srtm30_Pantanal', flags='pa', res='0:00:01')
grass.run_command('r.patch', input=patchlist, output='pantanal_ge', overwrite=True)

# rio claro
patchlist = ['rioclaro4', 'rioclaro3', 'rioclaro2', 'rioclaro1']
grass.run_command('g.region', raster='srtm30_RioClaro', flags='pa', res='0:00:01')
grass.run_command('r.patch', input=patchlist, output='rioclaro_ge', overwrite=True)

# iporanga
patchlist = ['iporanga1', 'iporanga2', 'iporanga3', 'iporanga4']
grass.run_command('g.region', raster='srtm30_Iporanga', flags='pa', res='0:00:01')
grass.run_command('r.patch', input=patchlist, output='iporanga_ge', overwrite=True)

# santa catarina
patchlist = ['floripa5', 'floripa4', 'floripa3', 'floripa2', 'floripa1', 'floripa_bk']
grass.run_command('g.region', raster='srtm30_SantaCatarina', flags='pa', res='0:00:01')
grass.run_command('r.patch', input=patchlist, output='santacatarina_ge', flags='z', overwrite=True)

# serra do mar
patchlist = ['paraty3', 'paraty2', 'paraty1', 'itatiaia4', 'itatiaia3', 'itatiaia2', 'itatiaia1', 'rift' ]
grass.run_command('g.region', raster='srtm30_SerraMar', flags='pa', res='0:00:01')
grass.run_command('r.patch', input=patchlist, output='serramar_ge', flags='z', overwrite=True)

# remove imported files
suffixes = ['.alpha','.red','.green','.blue']
for img in flist:
    img_imp = img[:-4]
    # grass.run_command('g.remove', flags='f', type='raster', name=img_imp)
    for suffix in suffixes:
        grass.run_command('g.remove', flags='f', type='raster', name=img_imp+suffix)


# export display as pngs
workdir = '/Volumes/MacintoshHD2/Dropbox/USP/projetosPesquisa/TanDEM-X/analysis_brazil/pngs/'
os.chdir(workdir)

maps_gearth = ['araca_ge','barcelos_ge','pantanal_ge','rioclaro_ge','iporanga_ge','santacatarina_ge','serramar_ge']

for map_ge in maps_gearth:
    grass.run_command('g.region', raster=map_ge, flags='pa', res='0:0:01')
    grass.run_command('d.mon', start='cairo', output=map_ge+'.png', resolution = '2', \
        height=500, width=500, overwrite=True)
    grass.run_command('d.rast', map=map_ge)
    grass.run_command('d.grid', size='0.25', text_color='white', border_color='white', flags='c')
    grass.run_command('d.mon', stop='cairo')

