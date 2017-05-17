#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from osgeo import gdal
from raster2polygon import raster2polygon
import shapefile

filename = "b20170517_densityMap_Korsblomst.tif"

ds = gdal.Open(filename, gdal.GA_ReadOnly)
#ds = gdal.Open("foulum_vest_4_index_temperature [°c].tif", gdal.GA_ReadOnly)
(X, deltaX, rotation, Y, rotation, deltaY) = ds.GetGeoTransform()
srs_wkt = ds.GetProjection()
Nx = ds.RasterXSize
Ny = ds.RasterYSize

rasterimage = ds.ReadAsArray()
polygons = raster2polygon(rasterimage, truncate=0.5,useGrid=False)


filename=filename.replace('.tif','')

#Write shapefile
sfwriter = shapefile.Writer(shapeType=shapefile.POLYGON)
sfwriter.field('DN','N','11')

for group in polygons.keys():
    listofpolygons=[p for p in polygons[group]]
    
    if len(listofpolygons)>0:
    
        #Convert to georefered numbers
        for ixpol,polygon in enumerate(listofpolygons):
            for ixpoint, point in enumerate(polygon):
                listofpolygons[ixpol][ixpoint] = [point[0]*deltaX+X,point[1]*deltaY+Y]
                    
        sfwriter.poly(shapeType=5, parts=listofpolygons)
        sfwriter.record(str(group))

sfwriter.save("%s.shp" % filename)


#Test
#sf = shapefile.Reader("%s.shp" % filename)
sf = shapefile.Reader("sample.shp")

#pol=sf.shapeRecords()


# create the PRJ file
with open("%s.prj" % filename, "w") as prj:
    prj.write(srs_wkt)
