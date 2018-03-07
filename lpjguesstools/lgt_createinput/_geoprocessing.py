# -*- coding: utf-8 -*-
"""lpjguesstools._geoprocessing: calculate slope, aspect, etc."""

import fiona
import logging
import numpy as np
import os
import rasterio
from rasterio.warp import calculate_default_transform
from rasterio.enums import Resampling
from rasterio.mask import mask
import scipy
import xarray as xr

from ._tpi import calculate_tpi

log = logging.getLogger(__name__)

# import constants
from . import NODATA
from . import ENCODING

def enlarge_array(a):
    """Pad grid boundaries for proper slope calc at adges."""

    ny, nx = a.shape
    b = np.zeros((ny + 2, nx + 2))
    b[1:-1,1:-1] = a  # Insert old grid in center

    # Assign boundary conditions - sides
    b[0, 1:-1] = a[0, :]
    b[-1, 1:-1] = a[-1, :]
    b[1:-1, 0] = a[:, 0]
    b[1:-1, -1] = a[:,-1]

    # Assign boundary conditions - corners
    b[0, 0] = a[0, 0]
    b[0, -1] = a[0, -1]
    b[-1, 0] = a[-1, 0]
    b[-1, -1] = a[-1, 0]
    return b


def calc_slope_components(dem, dx):
    """Calculate finite slopes."""
    # sx,sy = calcFiniteDiffs(elevGrid,dx)
    # calculates finite differences in X and Y direction using the 
    # 2nd order/centered difference method.
    # Applies a boundary condition such that the size and location 
    # of the grids in is the same as that out.

    # Assign boundary conditions
    dem_padded = enlarge_array(dem)

    #Compute finite differences
    Sx = (dem_padded[1:-1, 2:] - dem_padded[1:-1, :-2])/(2*dx)
    Sy = (dem_padded[2:,1:-1] - dem_padded[:-2, 1:-1])/(2*dx)
    return (Sx, Sy)


def calculate_utm_crs(lon, lat):
    """Calculate the UTM crs string from lon and lat coordinate."""
    code = 32700-int(np.round((45.0+lat)/90,0))*100+int(np.round((183.0+lon)/6,0))
    return 'EPSG:%d' % code    


def apply_mask(a, m):
    """Apply a mask from another masked_array."""
    return np.ma.masked_where(np.ma.getmask(m), a)


def calc_aspect(Sx, Sy):
    """Calculate aspect given X and Y slope components (unit: deg)."""
    aspect = np.rad2deg( np.arctan2(Sy, -Sx) )
    aspect = np.mod((450.0 - aspect), 360.)
    aspect[aspect==360] = 0
    return aspect


def calc_slope(Sx, Sy):
    """Calculate slope given X and Y slope components (unit: deg)."""
    return np.rad2deg(np.sqrt(Sx**2 + Sy**2))


def derive_coordinates(info):
    """Calculate tile lat lon information from GTiff info."""
    dx, _, leftc, _, dy, upperc, _, _, _ = info['transform']
    cellsx = info['width']
    cellsy = info['height']
    lowerc = upperc - (cellsy*abs(dy))
    lons = np.linspace(leftc, leftc+((cellsx+1)*dx), cellsx)
    lats = np.linspace(lowerc, lowerc+((cellsy+1)*abs(dy)), cellsy)
    # flipped lats
    return dict(lon=lons, lat=lats[::-1])


def create_tile(dem, dem_mask, slope, aspect, landform, info=None, source=None):
    """Create a tile dataset from dem, dem_mask, slope and aspect."""
    
    # if a rasterio transfrom info is passed
    if info != None:
        COORDS = derive_coordinates(info)
        DIMS = ['lat', 'lon']
    else:
        log.warn('No spatial information provided. Y-axis likely flipped.')
        COORDS={}
        DIMS=['dim_0', 'dim_1']
    
    # default mask
    m = np.ma.masked_where(dem_mask == 0, dem_mask)

    # special encoding (force output as Int16)
    ENCODING_INT = dict(ENCODING)
    ENCODING_INT.update({'dtype': np.int16})

    ds = xr.Dataset()
    ds['elevation'] = xr.DataArray(apply_mask(dem,m), coords=COORDS, dims=DIMS, 
                                   encoding=ENCODING_INT)
    ds['mask'] = xr.DataArray(dem_mask.astype('bool'), coords=COORDS, dims=DIMS, 
                                   encoding=ENCODING_INT)
    ds['slope'] = xr.DataArray(apply_mask(slope,m), coords=COORDS, dims=DIMS,
                                   encoding=ENCODING_INT)
    ds['aspect'] = xr.DataArray(apply_mask(aspect,m), coords=COORDS, dims=DIMS,
                                   encoding=ENCODING_INT)
    ds['landform'] = xr.DataArray(apply_mask(landform,m), coords=COORDS, dims=DIMS,
                                   encoding=ENCODING_INT)
    
    # add scale_factor to slope encoding
    ds['slope'].tile.update_encoding(dict(scale_factor=0.1))
    
    if source != None:
        ds.tile.set('source', source)
    return ds


def analyze_filename_dem(fname):
    """Analyze passed filename for zip components"""
    if fname[-4:] == '.zip':
        # default hgt in zip (SRTM1) - specific naming convention for SRTM1 1arc files
        bname = os.path.basename(fname).replace('.zip', '').split('.')[0] + '.hgt'
        fname = 'zip://%s!%s' % (fname, bname)
        source_name = bname
    else:
        if fname[-4:] not in ['.tif', '.hgt']:
            log.error('DEM file has unknown file suffix.')
            exit()
        source_name = os.path.basename(fname)
    return (fname, source_name)


def analyze_filename_shp(fname):
    """Analyze passed filename for zip components"""
    if fname[-4:] == '.zip':
        # default hgt in zip (SRTM1) - specific naming convention for SRTM1 1arc files
        bname = os.path.basename(fname).replace('.zip', '').split('.')[0] + '.shp'
        fname = 'zip://%s' % fname #, bname)
        source_name = bname
    else:
        if fname[-4:] not in ['.shp']:
            log.error('Shapefile file has unknown file suffix.')
            exit()
        source_name = os.path.basename(fname)
    return (fname, source_name)


def compute_spatial_dataset(fname_dem, fname_shp=None):
    """Take a GTiff file name and return a xarray datasets of dem, slope, 
    aspect and water mask layers."""
    
    fname_dem, source_name_dem = analyze_filename_dem(fname_dem)

    log.info('Opening file %s ...' % fname_dem)

    # open source GTiff file (in WGS84)
    with rasterio.open(fname_dem) as src:    
        msrc_kwargs = src.meta.copy()
        msrc_kwargs.update(count=5)
        msrc_kwargs.update(dtype='float64')
        msrc_kwargs.update(driver='GTiff')
        
        # read dem (as maskedarray) and create land mask (with gtiff nodata if present)
        dem = src.read(1, masked=True)
        # invert the bool array (0=missing, 1=valid)
        dem_mask = ~dem.mask #~np.ma.getmaskarray(dem)
        
        if fname_shp != None:
            fname_shp, source_name_shp = analyze_filename_shp(fname_shp)
            log.info("Masking water bodies")
            with fiona.open(fname_shp) as shp:
                geoms = [feature["geometry"] for feature in shp]
                dmask, _ = rasterio.mask.mask(src, geoms, nodata=NODATA, crop=False, invert=True)
                dmask = np.where(dmask == NODATA, 0, 1)
                # union of the dem mask and the waterfile mask
                dem_mask = dem_mask * dmask.squeeze()
        else:
            log.warn("No water mask shapefile found: %s" % fname_shp)
            
        # create a in-mem copy of input dem (4 bands: dem, mask, slope, aspect)
        with rasterio.io.MemoryFile() as memfile_geo1:
            with memfile_geo1.open(**msrc_kwargs) as ds_geo1:
                ds_geo1.write(dem.astype('float64'), 1)                 # dem
                ds_geo1.write(dem_mask.astype('float64'), 2)            # dem_mask
                ds_geo1.write(np.zeros_like(dem_mask, 'float64'), 3)    # slope
                ds_geo1.write(np.zeros_like(dem_mask, 'float64'), 4)    # aspect
                ds_geo1.write(np.zeros_like(dem_mask, 'float64'), 5)    # tpi300
                
                # derive utm projection (get center coordinate of tile)
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
                            dst_array = np.empty((height, width), dtype='float64')
                                
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

                        if dem_mask.sum() == 0:
                            dem_filled = dem.copy()                        
                        else:
                            log.debug('We have NoData gaps in DEM... filling')
                            # gapfill data
                            indices = scipy.ndimage.distance_transform_edt(np.invert(dem_mask.astype('bool')), 
                                return_distances=False, 
                                return_indices=True)
                            dem_filled = dem[tuple(indices)]

                        
                        # calculate slope & aspect
                        dx, dy = affine[0], affine[4]
                        if dx != -dy:
                            log.error("Cell sizes not square. Abort.")
                            exit()
                        
                        Sx, Sy = calc_slope_components(dem_filled, dx)
                        slope = calc_slope(Sx, Sy)
                        aspect = calc_aspect(Sx, Sy)
                        
                        # calculate tpi (now in utm)
                        landform = calculate_tpi(dem_filled, slope, 300, 
                                                 res=dx, TYPE='SIMPLE')

                        # write slope, aspect to ds_utm
                        ds_utm.write(slope.astype('float64'), 3)
                        ds_utm.write(aspect.astype('float64'), 4)
                        ds_utm.write(landform.astype('float64'), 5)
                        
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
                                    resampling=Resampling.nearest,
                                    num_threads=2)

                                dem_mask = ds_geo2.read(2).astype(bool)
                                dem_mask[:,-1] = dem_mask[:,-2]
                                dem = np.ma.masked_array(ds_geo2.read(1), mask=~dem_mask)
                                slope = np.ma.masked_array(ds_geo2.read(3), mask=~dem_mask)
                                aspect = np.ma.masked_array(ds_geo2.read(4), mask=~dem_mask)
                                landform = np.ma.masked_array(ds_geo2.read(5), mask=~dem_mask)
                                

    # create tile dataset    
    ds = create_tile(dem, dem_mask, slope, aspect, landform, 
                     info=msrc_kwargs, source=source_name_dem)
    
    return ds

# xarray-based methods

def classify_aspect(ds, TYPE='SIMPLE'):
    """Classify dataarray from continuous aspect to 1,2,3,4. or 1, 2"""
    
    aspect = ds['aspect'].to_masked_array()
    asp_cl = ds['aspect'].to_masked_array()
    
    # silence numpy warning in the comaprisons nan in masked_array
    import warnings
    warnings.filterwarnings("ignore",category=RuntimeWarning)
    if TYPE in ['WEISS', 'SIMPLE']:
        asp_cl[(aspect >= 315) | (aspect <  45)] = 1    # North
        asp_cl[(aspect >= 45)  & (aspect < 135)] = 2    # East
        asp_cl[(aspect >= 135) & (aspect < 225)] = 3    # South
        asp_cl[(aspect >= 225) & (aspect < 315)] = 4    # West
    #elif TYPE == 'SIMPLE':
    #    asp_cl[(aspect >= 270) | (aspect <  90)] = 1    # North
    #    asp_cl[(aspect  < 270) & (aspect >= 90)] = 3    # South
    else:
        log.error('Currently only classifiation schemes WEISS, SIMPLE supported.')

    # special encoding (force output as Int16)
    ENCODING_INT = dict(ENCODING)
    ENCODING_INT.update({'dtype': np.int16})    

    asp_cl = np.ma.masked_where(ds['mask'] == 0, asp_cl)
    da_asp_cl = xr.full_like(ds['aspect'], np.nan)
    ds['aspect_class'] = da_asp_cl
    ds['aspect_class'][:] = asp_cl
    ds['aspect_class'].tile.update_encoding(ENCODING_INT)
    return ds


def calculate_asp_slope(ds):
    ds['asp_slope'] = ds['slope'] * np.abs( np.cos(np.radians(ds['aspect'])) )
    
    # special encoding (force output as Int16)
    ENCODING_INT = dict(ENCODING)
    ENCODING_INT.update({'dtype': np.int16})
    ENCODING_INT.update({'scale_factor': 0.1})
    ds['asp_slope'].tile.update_encoding(ENCODING_INT)
    
    return ds


def classify_landform(ds, elevation_levels=[], TYPE='SIMPLE'):
    """Subdivide landform classes by aspect class."""        
    SHAPE = ds['mask'].shape
    lf_cl = np.ma.masked_array(np.ones_like(ds['mask'].values), mask=ds['mask'].values)
    
    # depending on classifiaction scheme we need different slope classes that 
    # have an aspect component
    if TYPE == 'SIMPLE':
        aspect_lf = [3]
    elif TYPE == 'WEISS':
        aspect_lf = [2,3,5]
    else:
        log.error('Currently only classifiation schemes WEISS, SIMPLE supported.')
    ds.tile.set('classification', TYPE.lower())
    
    aspect_lfs = (ds['aspect_class'].to_masked_array() > 0) & \
                  (np.in1d(ds['landform'].to_masked_array(), aspect_lf).reshape(SHAPE))
    
    lf_cl = np.ma.where(aspect_lfs, ds['landform'] * 10 + ds['aspect_class'],
                                    ds['landform'] * 10).filled(NODATA)
    lf_cl = np.ma.masked_where(ds['mask'] == 0, lf_cl)
    
    # if we have elevation levels subdivide the landform classes
    ele = ds['elevation'].to_masked_array()
    if len(elevation_levels) > 0:
        # add global elevation step attribute (second element, first is lower boundary)
        ds.tile.set('elevation_step', elevation_levels[1])

        for i, (lb, ub) in enumerate(zip(elevation_levels[:-1], elevation_levels[1:])):
            lf_cl = np.ma.where(((ele >= lb) & (ele < ub)), lf_cl + (i+1) * 100, lf_cl)   

    # special encoding (force output as Int16)
    ENCODING_INT = dict(ENCODING)
    ENCODING_INT.update({'dtype': np.int16})    

    lf_cl = np.ma.masked_where(ds['mask'] == 0, lf_cl)
    da_lf_cl = xr.full_like(ds['landform'], np.nan)
    ds['landform_class'] = da_lf_cl
    ds['landform_class'][:] = lf_cl
    ds['landform_class'].tile.update_encoding(ENCODING_INT)
    return ds
