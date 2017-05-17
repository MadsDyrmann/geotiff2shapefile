#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import cv2
import numpy as np

import itertools

#Todo: make it possible to specify grid and interval
def raster2polygon(rasterimage, truncate=1, intervals=[], minimumlength=4, useGrid=False):

    
    if truncate and len(intervals)>0:
        print('Both "truncate" and "intervals" are specified. Going to ignore specified intervals')
    
    
    #Truncate to nearest number divisible with truncate 
    if truncate:
        rasterimage=np.round(1.0*rasterimage/truncate)*truncate
    elif len(intervals)>0:
        intervals=[0]+intervals #make sure 0 is represented
        for i in range(len(intervals)-1):
            rasterimage[np.where(np.bitwise_and(rasterimage>=intervals[i] , rasterimage<intervals[i+1]))]=intervals[i]
        rasterimage[np.where(rasterimage>intervals[-1])]=intervals[-1]
    
                            
        
        
        
    if useGrid:
        rows= range(rasterimage.shape[0])
        cols= range(rasterimage.shape[1])
        contours = list(itertools.product(rows, cols))
        contours = [[[x[1],x[0]],[x[1]+1,x[0]],[x[1]+1,x[0]+1],[x[1],x[0]+1]] for x in contours]
        
        uniquePixels = np.unique(rasterimage)
        polygons={pixelvalue:[] for pixelvalue in uniquePixels[~ np.isnan(uniquePixels)]}
        
        for cnt in contours:
            pixelvalue = rasterimage[cnt[0][1],cnt[0][0]]
            if ~np.isnan(pixelvalue):
                polygons[pixelvalue].append(cnt)
          
        
    else: #make shapes based on equal intensities        
        polygons={}
        #Loop throuhg all pixel-values
        for pixelvalue in np.unique(rasterimage):
            imtmp=(rasterimage==pixelvalue).astype(np.uint8)
    
            contours, hierarchy = cv2.findContours(imtmp, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)        
            
            #remove small contours and make contours end in the start point
            contours=[c.squeeze().tolist() for c in contours if len(c)>minimumlength]
            contours=[c+[c[0]] for c in contours]
            
            if len(contours)>0:
                polygons[pixelvalue]=contours

    return polygons