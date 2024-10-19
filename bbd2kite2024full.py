#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:01:28 2024

@author: hog
"""
import logging
import os.path as op
import re
import shapefile
import utm
from scipy import stats
import geopandas
import pyogrio
import numpy as np
from kite.scene import Scene, SceneConfig

from timeit import default_timer as timer
from bbd2kite2024 import tl5_layer





log = logging.getLogger("bbd2kite24")

d2r = np.pi / 180.0
r2d = 180.0 / np.pi



# geofile  = '/media/hog/docCrucial1T/data/BBD2022/l2b_schleswig-holstein_clipped.gpkg'
# #geolayer = 'ASCE_015_01'
# geolayer = 'DESC_066_02'

# tl2_file = '/media/hog/docCrucial1T/data/data_sc/BBD_TL2/BBD_S1_PSI_TL2_SH_HH_SGD.gdb'
# tl2_layer= 'TL2__DELVY_20181025_V001__ASCE__Defo_Params'

# tl3_shp = '/media/hog/docCrucial1T/lab/BBD/BBDTL3/2020/BBD-TL3-BadSegeberg/BBD-TL3B-ASCE-117-01-BadSegebergUTM_sp.shp'



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
    
    log.info("Attempt loading data from %s", geofile)
    print('reading layer: ', layer)
    
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
    column_list_TL5b = ['PS_ID', 'nw_mv', 
               'var_mean_velocity', 'temp_coh', 
               'LOS_Up', 'LOS_North', 'LOS_East']
    
    column_list = column_list_TL5b
    
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
        
        # grond needs this the other way round -> *-1, -180 resp.
        data.phi   = np.arctan2(los_n, los_e)
        data.theta = np.arcsin(los_u)
        
        if gdf.iloc[0]['geometry'].geom_type == 'MultiPoint':
            dummy = gdf.explode(index_parts=True)
            gdf = dummy
            
        data.easts = gdf['geometry'].x.to_numpy()
        data.norths= gdf['geometry'].y.to_numpy()
        
        data.ps_mean_v   = gdf[column_list[1]].to_numpy()
        data.ps_mean_var = gdf[column_list[2]].to_numpy()
        
        zone   = gdf.crs.utm_zone[0:2]
        letter = gdf.crs.utm_zone[2]
        
    geopandas.options.io_engine = cached_engine
    
    return data, int(zone), letter


def bin_ps_data(data, bins=(200, 200)):
    log.debug("Binning mean velocity data...")
    bin_vels, edg_E, edg_N, _ = stats.binned_statistic_2d(
        data.easts, data.norths, data.ps_mean_v, statistic="mean", bins=bins
    )

    log.debug("Binning LOS angles...")
    bin_phi, _, _, _ = stats.binned_statistic_2d(
        data.easts, data.norths, data.phi, statistic="mean", bins=bins
    )
    bin_theta, _, _, _ = stats.binned_statistic_2d(
        data.easts, data.norths, data.theta, statistic="mean", bins=bins
    )

    log.debug("Binning mean velocity variance...")
    bin_mean_var, _, _, _ = stats.binned_statistic_2d(
        data.easts, data.norths, data.ps_mean_var, statistic="mean", bins=bins
    )

    
    data.bin_mean_var = np.float32(bin_mean_var).T
    data.bin_ps_mean_v = np.float32(bin_vels).T
    


    data.bin_phi = np.float32(bin_phi).T
    data.bin_theta = np.float32(bin_theta).T

    data.bin_edg_N = edg_N.T
    data.bin_edg_E = edg_E.T

    return data

def bbd2kite(filename, layer=None, px_size=(500, 500), 
             import_var=False, convert_m=True):
    """Convert BGR BodenBewegungsdienst PS velocity data to a Kite Scene

    Loads the mean PS velocities (from e.g. ``ps_plot(..., -1)``) from a
    BGR BodenBewegungsdienst, and grids the data into mean velocity bins.
    The LOS velocities will be converted to a Kite Scene
    (:class:`~kite.Scene`).

    :param filename: Name of the BGR BBD as ESRI shapefile.
    :type filename: str
    :param px_size: Size of pixels in North and East in meters.
        Default (500, 500).
    :type px_size: tuple
    :param convert_m: Convert displacement to meters, default True.
    :type convert_m: bool
    :param import_var: Import the mean velocity variance, this information
        is used by the Kite scene to define the covariance.
    :param import_var: bool
    """
    geofile_data = read_geofile(geofile=filename, 
                                layer=layer, 
                                engine='pyogrio')
    data   = geofile_data[0]
    zone   = geofile_data[1]
    letter = geofile_data[2]

    print('bbd2kite layer: ', layer)

    if convert_m:
        data.ps_mean_v /= 1e3
        data.ps_mean_var /= 1e3

    # lengthN = od.distance_accurate50m(
    #     data.bbox[1], data.bbox[0],
    #     data.bbox[3], data.bbox[0])
    # lengthE = od.distance_accurate50m(
    #     data.bbox[1], data.bbox[0],
    #     data.bbox[1], data.bbox[2])

    lengthE = data.bbox[2] - data.bbox[0]
    lengthN = data.bbox[3] - data.bbox[1]
    #print('Dimensions bbox [N, E]', lengthN, lengthE)
    #print(data.bbox)
   

    bins = (round(lengthE / px_size[0]), round(lengthN / px_size[1]))

    data = bin_ps_data(data, bins=bins) #added "data = " to stay correct

    log.debug("Setting up the Kite Scene")
    config = SceneConfig()
    #zone, letter = read_projection(filename)

    # ll = lower left
    llLat, llLon = utm.to_latlon(data.bbox[0], data.bbox[1], zone, letter)
    config.frame.llLat = llLat
    config.frame.llLon = llLon
    #print(llLat, llLon)

    config.frame.dE = (data.bin_edg_E[1] - data.bin_edg_E[0])
    config.frame.dN = (data.bin_edg_N[1] - data.bin_edg_N[0])
    config.frame.spacing = "meter"

    scene_name = op.basename(op.abspath(filename))
    #config.meta.scene_title = "%s (BBD import) %s " % layer % scene_name
    config.meta.scene_title = str(layer) + ' BBD import ' + str(scene_name)
    config.meta.scene_id = layer
    config.meta.satellite_name = 'Sentinel I'
    # config.meta.time_master = data.tmin.timestamp()
    # config.meta.time_slave = data.tmax.timestamp()
    

    scene = Scene(
        theta=data.bin_theta,
        phi=data.bin_phi,
        displacement=data.bin_ps_mean_v,
        config=config,
    )
    if import_var:
        scene.displacement_px_var = data.bin_mean_var
        
    # FIX LOS-TURN
    scene.phi = np.pi * (scene.phiDeg - 180.)/180.
    scene.phiDeg = scene.phiDeg - 180.
    scene.theta *= -1.
    scene.thetaDeg *= -1.
    scene.los.unitE *= -1.
    scene.los.unitN *= -1
    scene.los.unitU *= -1.
    return scene

#EXE##########################################################################    
start = timer()  
#df = read_geofile(geofile, geolayer, 'pyogrio')
#df = read_geofile(tl2_file, tl2_layer, 'pyogrio')
#df = read_geofile(tl3_shp, engine='pyogrio')

tl5_file  = '/media/hog/docCrucial1T/tools/nextcloud/kandidat/Roenne_overview/aoi_msc_gpk/tl5_l2b_aoi_msc_gpkg.gpkg'
tl5_file  = '/home/hog/data/mscthesis/aoi_msc_gpk/tl5_l2b_aoi_msc_gpkg.gpkg'
tl5_layer = 'tl5_d_139_01_mscaoi_sent1ABmv'


scene = bbd2kite(
    filename = tl5_file,
    layer = tl5_layer,
    px_size=(50,50),
    import_var=True,
    convert_m=not False,
)

end = timer()
# print('Time to prepare "spool" [s]', end - start)
out_file = '/home/hog/Downloads/'+tl5_layer
scene.save(out_file)
scene.spool()
