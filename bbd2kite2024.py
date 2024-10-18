#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:01:28 2024

@author: hog
"""

import geopandas
import pyogrio
import numpy as np
from timeit import default_timer as timer




geofile  = '/media/hog/docCrucial1T/data/BBD2022/l2b_schleswig-holstein_clipped.gpkg'
#geolayer = 'ASCE_015_01'
geolayer = 'DESC_066_02'

tl2_file = '/media/hog/docCrucial1T/data/data_sc/BBD_TL2/BBD_S1_PSI_TL2_SH_HH_SGD.gdb'
tl2_layer= 'TL2__DELVY_20181025_V001__ASCE__Defo_Params'

tl3_shp = '/media/hog/docCrucial1T/lab/BBD/BBDTL3/2020/BBD-TL3-BadSegeberg/BBD-TL3B-ASCE-117-01-BadSegebergUTM_sp.shp'

tl5_file  = '/media/hog/docCrucial1T/tools/nextcloud/kandidat/Roenne_overview/aoi_msc_gpk/tl5_l2b_aoi_msc_gpkg.gpkg'
tl5_layer = 'tl5_a_044_01_mscaoi'
# gdf = geopandas.read_file(
#     geofile,
#     layer = geolayer
# )

class DataStruct(dict):
    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, value):
        self[attr] = value
        
def read_geofile(geofile, layer=None, engine='fiona',
                 target_crs='EPSG:25832'):
    
    data = DataStruct()
    
    #specific to BBD TL5 internal GPKG files
    column_list_TL2 = ['PS_ID', 'Mean_Velo', 
               'Var_Mean_V', 'temp_coh', 
               'Los_Up', 'Los_North', 'Los_East']
    column_list_TL3 = ['PS_ID', 'Mean_Velo', 
               'Var_Mean', 'temp_coh', 
               'Los_Up', 'Los_North', 'Los_East']
    column_list_TL5 = ['PS_ID', 'mean_velocity', 
               'var_mean_velocity', 'temp_coh', 
               'LOS_Up', 'LOS_North', 'LOS_East']
    
    column_list = column_list_TL3
    
    cached_engine = geopandas.options.io_engine
    
    geopandas.options.io_engine = engine # PYthon OGR IO
    
    # include_fields for fiona 
    # columns for pyogrio 
    if geopandas.options.io_engine == 'pyogrio':
        gdf = geopandas.read_file(
            geofile,
            layer = layer,
            columns = column_list,
    )
    else:
        gdf = geopandas.read_file(
            geofile,
            layer = layer,
            include_fields = column_list,
    )
    
    if gdf.crs != target_crs:
        gdf = gdf.to_crs(target_crs)
    if gdf.crs.utm_zone is None:
        print('WARNING: CRS is NOT ETRS89/UTM z32N.')
    elif len(gdf.crs.utm_zone) != 3:
        print('WARNING: UTM zone naming convention probably violated.')
    else:    
        data.bbox = list(gdf.unary_union.envelope.bounds)
        los_u = gdf[column_list[4]].to_numpy()
        los_e = gdf[column_list[6]].to_numpy()
        los_n = gdf[column_list[5]].to_numpy()
        data.phi   = np.arctan2(los_n, los_e)
        data.theta = np.arcsin(los_u)
        
        data.easts = gdf['geometry'].x.to_numpy()
        data.norths= gdf['geometry'].y.to_numpy()
        
        data.ps_mean_v   = gdf[column_list[1]].to_numpy()
        data.ps_mean_var = gdf[column_list[2]].to_numpy()
        
        zone   = gdf.crs.utm_zone[0:2]
        letter = gdf.crs.utm_zone[2]
        
    geopandas.options.io_engine = cached_engine
    
    return data, zone, letter



    
start = timer()  
#df = read_geofile(geofile, geolayer, 'pyogrio')
#df = read_geofile(tl2_file, tl2_layer, 'pyogrio')
df = read_geofile(tl3_shp, engine='pyogrio')



end = timer()
print(end - start)