#!/usr/bin/env python
#-*- coding:utf-8 -*-
#
#
##############################################################################
# Script used in the paper:
# GAROPABA 3D
# by
# Carlos H. Grohmann et al - 2019/2020
# guano (at) usp (dot) br
# Institute of Energy and Environment - University of Sao Paulo
# 
# this is a helper script holding the azim() function, which will 
# ingest an ESRI shapefile with linear features and return the
# azimuth and length of each line
# 
# the fucntion does not check for projection, so UTM or other planar
# projection is recommended
##############################################################################

# import modules
import osgeo.ogr as ogr
import sys, math

# =============================================================================
def azim(infile):
    ''' this function reads an ESRI shapefile with linear features and return 
        the azimuth and length of each line.
        requires python GDAL/OGR libraries 
    '''

    # get the ogr driver
    ogrdriver = ogr.GetDriverByName("ESRI Shapefile")

    # open the data source
    dataSource = ogrdriver.Open(infile, 0)
    if dataSource is None:
       print 'File not found: ' + infile
       sys.exit(1)

    # get the data layer
    layer = dataSource.GetLayer()

    # count features in layer
    numFeatures = layer.GetFeatureCount()
    print('Linear features in file:', numFeatures)
    print('')

    #lists for vertices coordinates
    coord_X0 = []
    coord_X1 = []
    coord_Y0 = []
    coord_Y1 = []

    # loop through the features 
    f = layer.GetNextFeature()
    while f is not None:
        geometry = f.GetGeometryRef()
        x0 = geometry.GetX(0)
        x1 = geometry.GetX(1)
        y0 = geometry.GetY(0)
        y1 = geometry.GetY(1)
        coord_X0.append(x0)
        coord_X1.append(x1)
        coord_Y0.append(y0)
        coord_Y1.append(y1)
        f = layer.GetNextFeature()

    # need if looping again
    layer.ResetReading() 

    # close the data source
    dataSource.Destroy()

    #empty lists for results
    azimuth = []
    length = []

    # start calculations...
    # the loop is a bit cumbersome but it checks for special cases (N-S or E-W)
    for i in range(numFeatures):
        # case of a N-S line
        if coord_X0[i] == coord_X1[i]:
            if coord_Y0[i] > coord_Y1[i]: # line goes from north to south
                az = 180.
            else: # line goes from south to north
                az = 0.
            hyp = math.fabs(coord_Y0[i] - coord_Y1[i])
        # case of an E-W line
        elif coord_Y0[i] == coord_Y1[i]:
            if coord_X0[i] < coord_X1[i]:  # line goes from west to east
                az = 90.
            else: # line goes from east to west
                az = 270.
            hyp = math.fabs(coord_Y0[i] - coord_Y1[i])
        # other cases, first point of line W of second point
        elif coord_X0[i] < coord_X1[i]:
            m = (coord_Y1[i]-coord_Y0[i]) / (coord_X1[i]-coord_X0[i])
    #       print "X0 < X1: m = %f" % m
            dx = math.fabs(coord_X1[i]-coord_X0[i])
            dy = math.fabs(coord_Y1[i]-coord_Y0[i])
    #       print "dx = %f, dy = %f" % (dx,dy)
            hyp = math.hypot(dx,dy)
            arc = math.atan(dx/dy)
            azim = math.degrees(arc) # (arc*180)/PI
            if m < 0.:
                azim = 180. - azim
        # other cases, first point of line E of second point
        else:   # elif coord_X0[i] > coord_X1[i]
            m = (coord_Y1[i]-coord_Y0[i]) / (coord_X1[i]-coord_X0[i])
    #       print "X0 > X1: m = %f" % m
            dx = math.fabs(coord_X1[i]-coord_X0[i])
            dy = math.fabs(coord_Y1[i]-coord_Y0[i])
    #       print "dx = %f, dy = %f" % (dx,dy)
            hyp = math.hypot(dx,dy)
            arc = math.atan(dx/dy)
            azim = math.degrees(arc) # (arc*180)/PI
            if m < 0.:
                azim = 360. - azim
            else:
                azim = 180. + azim
        azimuth.append(azim)
        length.append(hyp)

    print('len azimuth:', len(azimuth))
    # return results
    return (azimuth,length)



