#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 12:43:03 2018

@author: ericdougherty
"""

import rasterio

file_list = ['Mean_Greenness_2009.tif',
             'Mean_Wetness_2009.tif',
             'Road_Density.tif']

with rasterio.open(file_list[0]) as src0:
    meta = src0.meta
    
meta.update(count=len(file_list))

with rasterio.open('Final_Predictors_2009.tif', 'w', **meta) as dst:
    for id, layer in enumerate(file_list):
        with rasterio.open(layer) as src1:
            dst.write_band(id + 1, src1.read(1))
          
####################################################

file_list = ['Mean_Greenness_2010.tif',
             'Mean_Wetness_2010.tif',
             'Road_Density.tif']

with rasterio.open(file_list[0]) as src0:
    meta = src0.meta
    
meta.update(count=len(file_list))

with rasterio.open('Final_Predictors_2010.tif', 'w', **meta) as dst:
    for id, layer in enumerate(file_list):
        with rasterio.open(layer) as src1:
            dst.write_band(id+1, src1.read(1))
