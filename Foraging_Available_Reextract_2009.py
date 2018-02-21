#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 09:04:47 2018

@author: ericdougherty
"""

def weighted_means(z, name_list, dist_vars):    
    temp_name = name_list[0][z]
    all_pts = gpd.GeoDataFrame.from_file('Movement_Shapefiles/' + ''.join([temp_name, '_Python.shp']))
    all_pts.crs = {'init' :'epsg:32733'}

    from scipy.stats import gamma
    shape = float(dist_vars[z]['gamma.shape'])
    rate = float(dist_vars[z]['gamma.rate'])
    rad = gamma.ppf(0.975, a=shape, scale=1/rate)

    big_buffers = all_pts.geometry.buffer(rad)
    #small_buffers = all_pts.geometry.buffer(30)
    #buffers_diff = big_buffers.difference(small_buffers)

    avail_green = []
    avail_wet = []
    avail_roads = []

    all_buff = big_buffers.geometry.values
    sample_iter = range(0, len(all_buff))
    from shapely.geometry import mapping
    for i in sample_iter:
        geoms = [mapping(big_buffers.geometry.values[i])]
        from rasterio.mask import mask       
        with rasterio.open("ENP_Predictors/Final_Predictors_2009.tif") as src:
            out_image, out_transform = mask(src, geoms, crop=True)

            no_data = -3.39999995e+38
            Green_band = out_image.data[0]
            Wet_band = out_image.data[1]
            Road_band = out_image.data[2]
            row, col = np.where(Green_band != no_data) 
            green = np.extract(Green_band != no_data, Green_band)
            wet = np.extract(Wet_band != no_data, Wet_band)
            roads = np.extract(Road_band != no_data, Road_band)

            from rasterio import Affine # or from affine import Affine
            T1 = out_transform * Affine.translation(0.5, 0.5) # reference the pixel centre
            rc2xy = lambda r, c: (c, r) * T1

            d = gpd.GeoDataFrame({'col':col,'row':row,
                                  'green':green,
                                  'wet':wet,
                                  'roads':roads})
            # coordinate transformation
            d['x'] = d.apply(lambda row: rc2xy(row.row,row.col)[0], axis=1)
            d['y'] = d.apply(lambda row: rc2xy(row.row,row.col)[1], axis=1)
            # geometry
            from shapely.geometry import Point
            d['geometry'] = d.apply(lambda row: Point(row['x'], row['y']), axis=1)

            from scipy.stats import gamma
            shape=float(dist_vars[z]['gamma.shape'])
            rate=float(dist_vars[z]['gamma.rate'])

            pt_iter = range(0, len(d))
            temp_weights = []
            temp_green_vals = []
            temp_wet_vals = []
            temp_roads_vals = []
            for j in pt_iter:
                temp_dist = d.loc[j].geometry.distance(all_pts['geometry'][i])
                weight = gamma.pdf(temp_dist, a=shape, scale=1/rate)
                temp_weights.append(weight)

                temp_green = d.loc[j].green
                temp_wet = d.loc[j].wet
                temp_roads = d.loc[j].roads
                temp_green_vals.append(temp_green)
                temp_wet_vals.append(temp_wet)
                temp_roads_vals.append(temp_roads)

            weighted_green = sum(temp_green_vals[g] * temp_weights[g] for g in range(len(temp_green_vals))) / sum(temp_weights)
            weighted_wet = sum(temp_wet_vals[g] * temp_weights[g] for g in range(len(temp_wet_vals))) / sum(temp_weights)
            weighted_roads = sum(temp_roads_vals[g] * temp_weights[g] for g in range(len(temp_roads_vals))) / sum(temp_weights)
            avail_green.append(weighted_green.mean())
            avail_wet.append(weighted_wet.mean())
            avail_roads.append(weighted_roads.mean())
            
    import pandas

    gdf = gpd.GeoDataFrame(geometry=all_pts['geometry'])
    x_test = gdf.geometry.apply(lambda p: p.x)
    y_test = gdf.geometry.apply(lambda p: p.y)
    out_df = pandas.DataFrame(data={"x": x_test, "y": y_test,
                            "green_avail": avail_green,
                            "wet_avail": avail_wet,
                            "roads_avail": avail_roads})
    return out_df;
            
import csv
import geopandas as gpd
import numpy as np
import rasterio

dist_vars = []
with open('Foraging_StepDist_Parameters.csv') as dist:
    reader = csv.DictReader(dist)
    for row in reader:
        dist_vars.append(row)

name_list = [('AG059_2009', 'AG061_2009', 'AG062_2009', 'AG063_2009', 
              'AG068_2009', 'AG063_2010', 'AG068_2010', 'AG252_2010', 
              'AG253_2010', 'AG255_2010', 'AG256_2010')]
    
#for q in len(name_list):
from joblib import Parallel, delayed
out_df = Parallel(n_jobs=5)(delayed(weighted_means)(z = q, name_list = name_list, dist_vars = dist_vars) for q in range(5))
for w in range(5):
    temp_df = out_df[w]
    temp_df.to_csv(''.join(['./', name_list[0][w], '_Foraging_Available.csv']), sep=',',index=False)
