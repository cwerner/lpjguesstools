# -*- coding: utf-8 -*-
"""lpjguesstools._geoprocessing: calculate slope, aspect, etc."""

import logging
import math
import numpy as np

import gdal
import gdalconst

from osgeo import osr, gdal, gdalconst

log = logging.getLogger(__name__)

def assign_boundary_cond(dem):
    """Pad grid boundaries for proper slope calc at adges."""

    # This creates a grid 2 rows and 2 columns larger than the input

    ny, nx = dem.shape
    dem_padded = np.zeros((ny + 2, nx + 2))
    dem_padded[1:-1,1:-1] = dem  # Insert old grid in center

    # Assign boundary conditions - sides
    dem_padded[0, 1:-1] = dem[0, :]
    dem_padded[-1, 1:-1] = dem[-1, :]
    dem_padded[1:-1, 0] = dem[:, 0]
    dem_padded[1:-1, -1] = dem[:,-1]

    # Assign boundary conditions - corners
    dem_padded[0, 0] = dem[0, 0]
    dem_padded[0, -1] = dem[0, -1]
    dem_padded[-1, 0] = dem[-1, 0]
    dem_padded[-1, -1] = dem[-1, 0]

    return dem_padded


def calc_slope_components(dem, dx):
    """Calculate finite slopes."""
    # sx,sy = calcFiniteDiffs(elevGrid,dx)
    # calculates finite differences in X and Y direction using the 
    # 2nd order/centered difference method.
    # Applies a boundary condition such that the size and location 
    # of the grids in is the same as that out.

    # Assign boundary conditions
    dem_padded = assign_bound_cond(dem)

    #Compute finite differences
    Sx = (dem_padded[1:-1, 2:] - dem_padded[1:-1, :-2])/(2*dx)
    Sy = (dem_padded[2:,1:-1] - dem_padded[:-2, 1:-1])/(2*dx)
    return (Sx, Sy)


def calc_aspect(Sx, Sy):
    """Calculate aspect given X and Y slope components (unit: deg)."""
    return np.rad2deg(np.arctan2(Sy, Sx))


def calc_slope(Sx, Sy):
    """Calculate slope given X and Y slope components (unit: deg)."""
    return np.rad2deg(np.sqrt(Sx**2 + Sy**2))


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


def process_dem(fname, verbose=False, force=False, shp_mask=None, GDALResampleAlg=gdal.GRA_Bilinear, dump=None):
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

    ## put input/slope/aspect into one 4-Band GDAL object
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

    # for DEBUG, also dump the file
    if dump != None:
        driver = gdal.GetDriverByName('GTiff')
        fdst = driver.Create(dump, src.RasterXSize,  src.RasterYSize, 4, gdal.GDT_Float64)
        fdst.SetGeoTransform(src.GetGeoTransform())
        fdst.SetProjection(src.GetProjection())
        fdst.GetRasterBand(1).WriteArray(src.GetRasterBand(1).ReadAsArray())
        fdst.GetRasterBand(2).WriteArray(dst_slp.GetRasterBand(1).ReadAsArray())
        fdst.GetRasterBand(3).WriteArray(dst_asp.GetRasterBand(1).ReadAsArray())
        fdst.GetRasterBand(4).WriteArray(rmask)

        
    return dst

def read_processed_geotiff(fname):
    src = gdal.Open(fname, gdalconst.GA_ReadOnly)
    return src 