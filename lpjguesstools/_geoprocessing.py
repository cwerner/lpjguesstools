# -*- coding: utf-8 -*-
"""lpjguesstools._geoprocessing: calculate slope, aspect, etc."""

import logging
import math
import numpy as np

import gdal
import gdalconst

from osgeo import osr, gdal, gdalconst

log = logging.getLogger(__name__)


def assert_gdal_dataset(d, fname='object'):
    """check if we have a valid gdal.dataset instance"""
    if not isinstance(d, gdal.Dataset):
        print("ERROR: %s not of type 'gdal.Dataset'" % fname)
        sys.exit(-1)

def assert_ogr_datasource(d, fname='object'):
    """check if we have a valid ogr.datasource instance"""
    if not isinstance(d, ogr.DataSource):
        print("ERROR: %s not of type 'ogr.DataSource'" % fname)
        sys.exit(-1)


def project2UTM(src, GDALResampleAlg, verbose=False, force=False):
    """Reproject latlon to UTM"""
    
    assert_gdal_dataset(src, fname='project2UTM.src')

    ## get the spatial reference
    src_wkt = src.GetProjectionRef()
    src_proj = osr.SpatialReference()
    src_proj.ImportFromWkt(src_wkt)

    ## check if unit is already 'metre'
    if src_proj.GetAttrValue('UNIT') == 'metre' and not force:
        print("Unit is already 'metre', returning unmodified!")
        return(src)

    ## create the UTM information
    src_gtrans = src.GetGeoTransform()
    src_ix = src.RasterXSize
    src_iy = src.RasterYSize
    src_cx = src_gtrans[0] + src_ix * src_gtrans[1] / 2.0
    src_cy = src_gtrans[3] + src_iy * src_gtrans[5] / 2.0
    ## TODO: check if longitude is larger than 180
    utm_zone = math.floor((src_cx + 186) / 6)
    dst_epsg = int(32000 + (6 + (1 - np.sign(src_cy)) / 2) * 100 + utm_zone)
    if verbose:
        print("Source projection:   {}".format(src_proj))

    ## Create the destination GDAL object
    dst_proj = osr.SpatialReference()
    dst_proj.ImportFromEPSG(dst_epsg)
    if verbose:
        print("Destination EPSG code:  %i" % dst_epsg)
        print("Destination Projection: {}".format(dst_proj))
    ## transform
    dst = gdal.AutoCreateWarpedVRT(src, src_proj.ExportToWkt(), dst_proj.ExportToWkt(), GDALResampleAlg)
    return(dst)


def rasterizeMask(mask, dst):
    """Rasterize mask shapefile (water bodies)"""

    assert_gdal_dataset(dst, fname='rasterizeMask.dst')
    assert_ogr_datasource(mask, fname='rasterizeMask.mask')

    lmask = mask.GetLayer()

    mem_drv = gdal.GetDriverByName('MEM')
    dest_ds = mem_drv.Create('', dst.RasterXSize, dst.RasterYSize, 1, gdal.GDT_Float64)
    dest_ds.SetGeoTransform(dst.GetGeoTransform())
    dest_ds.SetProjection(dst.GetProjection())

    ## Rasterize
    ## TODO: Solve the Projection warning here
    gdal.RasterizeLayer(dest_ds, [1], lmask)
    return dest_ds.ReadAsArray(1)


def process_dem(fname, verbose=False, force=False, shp_mask=None, GDALResampleAlg=gdal.GRA_Bilinear):
    """Calculate slope, aspect and tpi from source DEM"""
    
    #src = gdal.Open('/vsizip/' + os.path.join(tile, bfilename), gdalconst.GA_ReadOnly)
    print fname
    src = gdal.Open(fname, gdalconst.GA_ReadOnly)
    
    ## test if 'src' is a valid GDAL
    assert_gdal_dataset(src, fname='processDEM.src')

    ## reproject to UTM
    utm = project2UTM(src, GDALResampleAlg, verbose=verbose, force=force)

    ## Create slope and aspect in UTM coordinates
    ## TODO: maybe make computeEdges customizable
    utm_slp = gdal.DEMProcessing('', utm, 'slope', format='MEM', computeEdges=True)
    utm_asp = gdal.DEMProcessing('', utm, 'aspect', format='MEM', computeEdges=True)

    ## Reproject slope/aspect back to input projection
    mem_drv = gdal.GetDriverByName('MEM')

    dst_slp = mem_drv.Create('', src.RasterXSize,  src.RasterYSize, 1, gdal.GDT_Float64)
    dst_slp.SetGeoTransform(src.GetGeoTransform())
    dst_slp.SetProjection(src.GetProjection())  
    gdal.Warp(dst_slp, utm_slp)

    dst_asp = mem_drv.Create('', src.RasterXSize,  src.RasterYSize, 1, gdal.GDT_Float64)
    dst_asp.SetGeoTransform(src.GetGeoTransform())
    dst_asp.SetProjection(src.GetProjection())  
    gdal.Warp(dst_asp, utm_asp)

    ## put input/slope/aspect into one 3-Band GDAL object
    dst = mem_drv.Create('', src.RasterXSize,  src.RasterYSize, 4, gdal.GDT_Float64)
    dst.SetGeoTransform(src.GetGeoTransform())
    dst.SetProjection(src.GetProjection())
    
    dst.GetRasterBand(1).WriteArray(src.GetRasterBand(1).ReadAsArray())
    dst.GetRasterBand(2).WriteArray(dst_slp.GetRasterBand(1).ReadAsArray())
    dst.GetRasterBand(3).WriteArray(dst_asp.GetRasterBand(1).ReadAsArray())

    if shp_mask != None:
        rmask = rasterizeMask(mask, src)
    else:
        rmask = np.zeros_like(src.GetRasterBand(1).ReadAsArray())
    dst.GetRasterBand(4).WriteArray(rmask)

    return dst
