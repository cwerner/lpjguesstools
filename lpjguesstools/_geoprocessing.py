# -*- coding: utf-8 -*-
"""lpjguesstools._geoprocessing: calculate slope, aspect, etc."""

import logging
import math
import numpy as np
import xarray as xr
import gdal
import gdalconst

# read lat/lon GTiff
import fiona
import rasterio
from rasterio.warp import calculate_default_transform
from rasterio.enums import Resampling
from rasterio.mask import mask
import scipy

from osgeo import osr, gdal, gdalconst

log = logging.getLogger(__name__)

NODATA = -9999.0
defaultAttrsDA = {
        '_FillValue':    NODATA,
        'missing_value': NODATA
        }



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
    dem_padded = assign_boundary_cond(dem)

    #Compute finite differences
    Sx = (dem_padded[1:-1, 2:] - dem_padded[1:-1, :-2])/(2*dx)
    Sy = (dem_padded[2:,1:-1] - dem_padded[:-2, 1:-1])/(2*dx)
    return (Sx, Sy)


def calculate_utm_crs(lon, lat):
    """Calculate the UTM crs string from lon and lat coordinate."""
    code = 32700-int(np.round((45.0+lat)/90,0))*100+int(np.round((183.0+lon)/6,0))
    return 'EPSG:%d' % code    


def calc_aspect(Sx, Sy):
    """Calculate aspect given X and Y slope components (unit: deg)."""
    aspect = np.rad2deg( np.arctan2(Sy, -Sx) )
    aspect = np.mod((450.0 - aspect), 360.)
    aspect[aspect==360] = 0
    
    # reclass
    asp_cl = np.zeros_like(aspect)
    asp_cl[(aspect >= 315) | (aspect <  45)] = 1    # North
    asp_cl[(aspect >= 45)  & (aspect < 135)] = 2    # East
    asp_cl[(aspect >= 135) & (aspect < 225)] = 3    # South
    asp_cl[(aspect >= 225) & (aspect < 315)] = 4    # West
    return (aspect, asp_cl)

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

def create_dem_dataset(dem, dem_mask, slope, aspect, info=None):
    """Create a datasets from dem, dem_mask, slope and aspect."""
    
    # if a rasterio transfrom info is passed
    info = None
    if info != None:
        dx, _, leftc, _, dy, lowerc, _, _, _ = info['transform']
        cellsx = info['width']
        cellsy = info['height']
        lons = np.linspace(leftc, leftc+((cellsx+1)*dx), cellsx)
        lats = np.linspace(lowerc, lowerc+((cellsy+1)*-dy), cellsy)
        
        COORDS = dict(lat=lats[::-1], lon=lons)
        DIMS = ['lat', 'lon']
    else:
        log.warn('No spatial information provided. y-axis likely flipped.')
        COORDS={}
        DIMS=['dim_0', 'dim_1']
    
    ds = xr.Dataset()
    ds['dem'] = xr.DataArray(dem, coords=COORDS, dims=DIMS).fillna(NODATA)
    ds['mask'] = xr.DataArray(dem_mask.astype('bool'), coords=COORDS, dims=DIMS)
    ds['slope'] = (xr.DataArray(slope, coords=COORDS, dims=DIMS) * ds['mask']).fillna(NODATA)
    ds['aspect'] = (xr.DataArray(aspect, coords=COORDS, dims=DIMS) * ds['mask']).fillna(NODATA)
    
    for v in ['dem', 'mask', 'slope', 'aspect']:
        ds[v].attrs.update(defaultAttrsDA)
    return ds


def compute_mask_slope_aspect(fname, fname_shp=None):
    """Take a GTiff file name and return a xarray datasets of dem, slope, 
    aspect and water mask layers."""

    # open source GTiff file (in WGS84)
    with rasterio.open(fname) as src:    
        msrc_kwargs = src.meta.copy()
        msrc_kwargs.update(count=4)
        msrc_kwargs.update(dtype='Float64')
        
        # read the dem band and shp polygons
        dem = src.read(1)

        if fname_shp != None:
            with fiona.open(fname_shp) as shp:
                geoms = [feature["geometry"] for feature in shp]

                # mask dem
                dem_mask, out_transform = rasterio.mask.mask(src, geoms, crop=False, invert=True)
                dem_mask = dem_mask.squeeze()
        else:
            dem_mask = np.zeros_like(dem)
            
        # create a in-mem copy of input dem (4 bands: dem, mask, slope, aspect)
        with rasterio.io.MemoryFile() as memfile_geo1:
            with memfile_geo1.open(**msrc_kwargs) as ds_geo1:
                print ds_geo1.shape
                ds_geo1.write(dem.astype('Float64'), 1)
                ds_geo1.write(dem_mask.astype('Float64'), 2)
                ds_geo1.write(np.zeros_like(dem_mask, 'Float64'), 3)
                ds_geo1.write(np.zeros_like(dem_mask, 'Float64'), 4)
                
                # utm projection
                # get center coordinate of tile
                lon, lat = ds_geo1.transform * (ds_geo1.width * 0.5, ds_geo1.height * 0.5)
                dst_crs = calculate_utm_crs(lon, lat)
                               
                # calc transform for UTM dst
                affine, width, height = calculate_default_transform(
                    ds_geo1.crs, dst_crs, ds_geo1.width, ds_geo1.height, *ds_geo1.bounds)
                
                # modify meta-data for dst after transform
                kwargs = ds_geo1.meta.copy()
                kwargs.update({
                    'crs': dst_crs,
                    'transform': affine,
                    'affine': affine,
                    'width': width,
                    'height': height
                })
                
                # reproject to another in-mem file
                with rasterio.io.MemoryFile() as memfile_utm:
                    with memfile_utm.open(**kwargs) as ds_utm:
                        for i in range(1, ds_geo1.count + 1):
                            dst_array = np.empty((height, width), dtype='Float64')
                                
                            rasterio.warp.reproject(
                                source=ds_geo1.read(i),
                                src_crs=ds_geo1.crs,
                                src_transform=ds_geo1.transform,
                                destination=dst_array, #_utm,
                                dst_transform=affine,
                                dst_crs=dst_crs,
                                resampling=Resampling.bilinear,
                                num_threads=2)

                            ds_utm.write(dst_array, i)
                        
                        # buffer dem at mask edge
                        dem  = ds_utm.read(1)
                        dem_mask = ds_utm.read(2)
                        
                        # gapfill data
                        indices = scipy.ndimage.distance_transform_edt(np.invert(dem_mask.astype('bool')), 
                            return_distances=False, 
                            return_indices=True)
                        dem_filled = dem[tuple(indices)]
                        
                        # calculate slope & aspect
                        dx, dy = affine[0], affine[4]
                        if dx != -dy:
                            log.error("Cell sizes not square. Abort.")
                        
                        # alternative slope/ aspect code
                        #slope = slope_a(dem_filled, cell_size=dx, degrees=True, keepdims=False)
                        #aspect = aspect_a(dem_filled)
                        
                        Sx, Sy = calc_slope_components(dem_filled, dx)
                        slope = calc_slope(Sx, Sy)
                        aspect, aspect_classed = calc_aspect(Sx, Sy)
                        
                        # write slope, aspect to ds_utm
                        ds_utm.write(slope.astype('Float64'), 3)
                        ds_utm.write(aspect.astype('Float64'), 4)
                        
                        # transform back to LatLon
                        with rasterio.io.MemoryFile() as memfile_geo2:
                            with memfile_geo2.open(**ds_geo1.meta.copy()) as ds_geo2:    

                                # take info from in-memory geo file
                                dst_crs = ds_geo1.crs
                                dst_transform = ds_geo1.transform
                                dst_height = ds_geo1.height
                                dst_width = ds_geo1.width

                                out_kwargs = ds_utm.profile.copy()
                                out_kwargs.update({
                                    'crs': dst_crs,
                                    'transform': dst_transform,
                                    'width': dst_width,
                                    'height': dst_height
                                })    

                                rasterio.warp.reproject(
                                    source=rasterio.band(ds_utm, list(range(1, ds_utm.count + 1))),
                                    destination=rasterio.band(ds_geo2, list(range(1, ds_utm.count + 1))),
                                    src_transform=ds_utm.transform,
                                    src_crs=ds_utm.crs,
                                    #src_nodata=ds_geods_utmsrc_nodata,
                                    dst_transform=out_kwargs['transform'],
                                    dst_crs=out_kwargs['crs'],
                                    #dst_nodata=dst_nodata,
                                    resampling=Resampling.bilinear,
                                    num_threads=2)

                                dem_mask = ds_geo2.read(2).astype(bool)
                                dem = np.ma.masked_array(ds_geo2.read(1), mask=~dem_mask)
                                slope = np.ma.masked_array(ds_geo2.read(3), mask=~dem_mask)
                                aspect = np.ma.masked_array(ds_geo2.read(4), mask=~dem_mask)

    # create dataset    
    ds = create_dem_dataset(dem, dem_mask, slope, aspect, info=msrc_kwargs)
    
    return ds


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