#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 13:33:23 2017

@author: mads
"""

import pandas as pd
import matplotlib.path as mplPath
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import shapefile
import utm


def inpolygon(point,polygon):
    bbPath = mplPath.Path(np.array(polygon))
    return bbPath.contains_point(point)


def getPolygon(coordinate):
    sf = shapefile.Reader("/media/mads/Eksternt drev/MarkPolygoner/Marker2016oShape/Marker_download_offentlig_20160809.shp")
    polygons = sf.shapeRecords()
    
    coordinateUTM = utm.from_latlon(coordinate[0],coordinate[1])
    

    polygons=[polygon for polygon in polygons if (polygon.shape.bbox[0]<coordinateUTM[0] and polygon.shape.bbox[2]>coordinateUTM[0])]
    polygons=[polygon for polygon in polygons if (polygon.shape.bbox[1]<coordinateUTM[1] and polygon.shape.bbox[3]>coordinateUTM[1])]

    
    for polygon in polygons:
        polyPoints=polygon.shape.points
        if inpolygon(coordinateUTM,polyPoints):
            return polygon


def getPixelsInPolygon(polygon,scale=1):
    print('scale er i forhold til meter. 1=1meter, 0.001 = 1mm, 10=10meter etc.')
    #The point that is considered (0,0) in image
    
    #polygon.shape.points=[utm.to_latlon(point[0],point[1],32,'u') for point in polygon.shape.points]
    
    minlat= min([p[0] for p in polygon.shape.points])
    maxlat= max([p[0] for p in polygon.shape.points])
    minlong=min([p[1] for p in polygon.shape.points])
    maxlong=max([p[1] for p in polygon.shape.points])
    
                                  
    offsetCoordinate=(minlat, minlong)
    
    ncolumns=np.ceil((maxlong-minlong)/scale).astype(np.int)
    nrows=np.ceil((maxlat-minlat)/scale).astype(np.int)
    
    fieldRaster=np.zeros((ncolumns,nrows))
    
    
    for row in range(fieldRaster.shape[0]):
        for column in range(fieldRaster.shape[1]):
            coordinate=(row+offsetCoordinate[0],column+offsetCoordinate[1])
            if inpolygon(coordinate,polygon.shape.points):
                fieldRaster[row,column]=1
    
    
    
    pixelwidthInMeters=scale
    
    
    return fieldRaster, pixelwidthInMeters, offsetCoordinate


def distanceRaster(fieldraster,pixelwidthInMeters,offsetCoordinate,listOfDetections, species):
    
    detectionArea=0.25
    
    
    fieldDensities={}
    for species in species:
        fieldDensities[species]=np.nan*fieldraster
                      
    #For all detections determine the density for a pixel in the raster
    rows,columns=np.where(fieldraster>0)
    for row, col in zip(rows,columns):
        #Find the detections in this pixel
        r=row+offsetCoordinate[0]
        c=col+offsetCoordinate[1]
        
        #Index of all detected plants that falls in the index
        inpixel = [1 if (x[0]>r and x[0]<r+pixelwidthInMeters and x[1]>c and x[1]<c+pixelwidthInMeters) else 0 for x in plantsDF['UTM']]

        
        if np.any(inpixel):
            #Create dataframe with plants in pixel
            inpixelDF=plantsDF.ix[np.where(inpixel)]
            #Loop over all unique species that is in this pixel
            
            #Count the number of occurences and set the pixel value to the number of occurances.
    
    
            print('All the plants in this pixel should ')
    
    
    
    
    rows,columns=np.where(fieldRaster>0)
    for row, column in zip(rows,columns):
        pass
    
    return fieldraster



def createMap(plantDF,listOfWeedSpecies,scaleinmeters=1):
    
    #Convert plantDF to UTM
    plantDF['UTM'] = plantDF.apply(lambda row: utm.from_latlon(row['longitude'], row['longitude']), axis=1)

    
    polygon = getPolygon((plantDF.iloc[0].longitude,plantDF.iloc[0].latitude))
    
    fieldRaster, pixelwidthInMeters, offsetCoordinate = getPixelsInPolygon(polygon,scale=1)
    
    distanceRaster(fieldraster=fieldRaster,pixelwidthInMeters=pixelwidthInMeters,offsetCoordinate=offsetCoordinate, species=listOfWeedSpecies, listOfDetections=plantDF)
    
    pass
#    print('convert polygons to latitude-longitude')
#    polygon.shape.points=[utm.to_latlon(point[0],point[1],32,'u') for point in polygon.shape.points]
    #lons=[p[0] for p in polygon.shape.points]
    #lats=[p[1] for p in polygon.shape.points]

    
#    m = Basemap(projection='merc',llcrnrlat=min(lats),urcrnrlat=max(lats),llcrnrlon=min(lons),urcrnrlon=max(lons),lat_ts=20,resolution='h')



    verts = np.array(polygon.shape.points)
    coll = PolyCollection(np.expand_dims(verts,0),edgecolors='red',facecolor=None)

#    x,y = m(lons,lats)


    fig, ax = plt.subplots()
    ax.add_collection(coll)
    ax.autoscale_view()
    
    
    pass





csvPath = '/media/mads/79131cf3-d458-45c7-bb29-830bc120a265/Phd/Software/DeepNetworks/SingleImageClassificationUsingCaffe/Results_2017_04_03_242.csv'

plantsDF = pd.read_csv(csvPath,sep=';')
listOfWeedSpecies = list(set(plantsDF['predictSpecies']))
createMap(plantsDF,listOfWeedSpecies,scaleinmeters=1)


