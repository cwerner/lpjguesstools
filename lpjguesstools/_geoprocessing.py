# -*- coding: utf-8 -*-
"""lpjguesstools._geoprocessing: calculate slope, aspect, etc."""

import fiona
import logging
import math
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform
from rasterio.enums import Resampling
from rasterio.mask import mask
import os
import scipy
import xarray as xr

from _tpi import calculate_tpi

log = logging.getLogger(__name__)

# import constants
from . import NODATA
from . import defaultAttrsDA

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
    return aspect


def calc_slope(Sx, Sy):
    """Calculate slope given X and Y slope components (unit: deg)."""
    return np.rad2deg(np.sqrt(Sx**2 + Sy**2))


def create_dem_dataset(dem, dem_mask, slope, aspect, landform, info=None, source=None):
    """Create a datasets from dem, dem_mask, slope and aspect."""
    
    # if a rasterio transfrom info is passed
    if info != None:
        dx, _, leftc, _, dy, upperc, _, _, _ = info['transform']
        cellsx = info['width']
        cellsy = info['height']
        lowerc = upperc - (cellsy*abs(dy))
        lons = np.linspace(leftc, leftc+((cellsx+1)*dx), cellsx)
        lats = np.linspace(lowerc, lowerc+((cellsy+1)*abs(dy)), cellsy)
        
        COORDS = dict(lat=lats[::-1], lon=lons)
        DIMS = ['lat', 'lon']
    else:
        log.warn('No spatial information provided. y-axis likely flipped.')
        COORDS={}
        DIMS=['dim_0', 'dim_1']
    
    ds = xr.Dataset()
    ds['elevation'] = xr.DataArray(dem, coords=COORDS, dims=DIMS).fillna(NODATA)
    ds['mask'] = xr.DataArray(dem_mask.astype('bool'), coords=COORDS, dims=DIMS)
    ds['slope'] = (xr.DataArray(slope, coords=COORDS, dims=DIMS) * ds['mask']).fillna(NODATA)
    ds['aspect'] = (xr.DataArray(aspect, coords=COORDS, dims=DIMS) * ds['mask']).fillna(NODATA)
    ds['landform'] = (xr.DataArray(landform, coords=COORDS, dims=DIMS) * ds['mask']).fillna(NODATA)
    
    for v in ['elevation', 'slope', 'aspect', 'landform']:
        ds[v].attrs.update(defaultAttrsDA)
    
    if source != None:
        set_global_attr(ds, 'source', source)
    return ds


def analyze_filename(fname):
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


def compute_spatial_dataset(fname, fname_shp=None):
    """Take a GTiff file name and return a xarray datasets of dem, slope, 
    aspect and water mask layers."""
    
    fname, source_name = analyze_filename(fname)

    log.info('Opening file %s ...' % fname)

    # open source GTiff file (in WGS84)
    with rasterio.open(fname) as src:    
        msrc_kwargs = src.meta.copy()
        msrc_kwargs.update(count=5)
        msrc_kwargs.update(dtype='float64')
        msrc_kwargs.update(driver='GTiff')
        
        # read dem (as maskedarray) and create land mask (with gtiff nodata if present)
        dem = src.read(1, masked=True)
        dem_mask = ~np.ma.getmaskarray(dem)
        
        if fname_shp != None:
            log.info("Masking water bodies")
            with fiona.open(fname_shp) as shp:
                geoms = [feature["geometry"] for feature in shp]
                dmask, _ = rasterio.mask.mask(src, geoms, nodata=0, crop=False, invert=True)
                dmask = np.where(dmask > 0, 1, 0)
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
                            log.info('We have ZERO gaps... Proceed.')
                            dem_filled = dem.copy()                        
                        else:
                            log.info('We have gaps... Filling.')
                            # gapfill data
                            indices = scipy.ndimage.distance_transform_edt(np.invert(dem_mask.astype('bool')), 
                                return_distances=False, 
                                return_indices=True)
                            dem_filled = dem[tuple(indices)]

                        
                        # calculate slope & aspect
                        dx, dy = affine[0], affine[4]
                        if dx != -dy:
                            log.error("Cell sizes not square. Abort.")
                        
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
                                dem = np.ma.masked_array(ds_geo2.read(1), mask=~dem_mask)
                                slope = np.ma.masked_array(ds_geo2.read(3), mask=~dem_mask)
                                aspect = np.ma.masked_array(ds_geo2.read(4), mask=~dem_mask)
                                landform = np.ma.masked_array(ds_geo2.read(5), mask=~dem_mask)
                                
                                #log.info("Dumping GTiff 1arc file for debugging.")
                                #print msrc_kwargs
                                #with rasterio.open('test.tiff', 'w', **msrc_kwargs) as dst:
                                #    for i in range(1,6):
                                #        dst.write(ds_geo2.read(i), i)

    # create dataset    
    ds = create_dem_dataset(dem, dem_mask, slope, aspect, landform, 
                            info=msrc_kwargs, source=source_name)
    
    return ds

# xarray-based methods

def get_center_coord(ds):
    """Return the (lon, lat) of dataset (center)"""
    lat_c = min(ds.lat.values) + (max(ds.lat.values) - min(ds.lat.values)) * 0.5
    lon_c = min(ds.lon.values) + (max(ds.lon.values) - min(ds.lon.values)) * 0.5
    return (lon_c, lat_c)
        

def classify_aspect(ds, TYPE='SIMPLE'):
    """Classify dataarray from continuous aspect to 1,2,3,4. or 1, 2"""
    
    # TODO: Check if we need only northern and southern aspect classes for
    #       implementation     
    aspect = ds['aspect'].to_masked_array()
    asp_cl = ds['aspect'].to_masked_array()
    
    if TYPE == 'WEISS':
        asp_cl[(aspect >= 315) | (aspect <  45)] = 1    # North
        asp_cl[(aspect >= 45)  & (aspect < 135)] = 2    # East
        asp_cl[(aspect >= 135) & (aspect < 225)] = 3    # South
        asp_cl[(aspect >= 225) & (aspect < 315)] = 4    # West
    elif TYPE == 'SIMPLE':
        asp_cl[(aspect >= 270) | (aspect <  90)] = 1    # North
        asp_cl[(aspect  < 270) & (aspect >= 90)] = 3    # South
    else:
        log.error('Currently only classifiation schemes WEISS, SIMPLE supported.')

    asp_cl = np.ma.masked_where(ds['mask'] == 0, asp_cl).filled(NODATA)
    ds['aspect_class'] = xr.full_like(ds['aspect'], NODATA)
    ds['aspect_class'][:] = asp_cl
    ds['aspect_class'].attrs.update(defaultAttrsDA)
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
    set_global_attr(ds, 'lgt.classification', TYPE.lower())
    
    aspect_lfs = (ds['aspect_class'].to_masked_array() > 0) & \
                  (np.in1d(ds['landform'].to_masked_array(), aspect_lf).reshape(SHAPE))
    
    lf_cl = np.ma.where(aspect_lfs, ds['landform'] * 10 + ds['aspect_class'],
                                    ds['landform'] * 10).filled(NODATA)
    lf_cl = np.ma.masked_where(ds['mask'] == 0, lf_cl)
    
    # if we have elevation levels subdivide the landform classes
    ele = ds['elevation'].to_masked_array()
    if len(elevation_levels) > 0:
        # add global elevation step attribute (second element, first is lower boundary)
        set_global_attr(ds, 'lgt.elevation_step', "%s" % elevation_levels[1])

        for i, (lb, ub) in enumerate(zip(elevation_levels[:-1], elevation_levels[1:])):
            lf_cl = np.ma.where(((ele >= lb) & (ele < ub)), lf_cl + (i+1) * 100, lf_cl)   
    
    ds['landform_class'] = xr.full_like(ds['landform'], NODATA)
    ds['landform_class'][:] = lf_cl.filled(NODATA)
    ds['landform_class'].attrs.update(defaultAttrsDA)
    return ds


def is_empty_tile(ds):
    """Check if this tile has no data (sum(mask)==0)."""
    if ds['mask'].sum() == 0:
        return True
    return False


def split_srtm1_dataset(ds):
    """Split a 1arc SRTM1 dataset into 4 0.5x0.5 tiles."""
    lats_ix = np.arange(len(ds['lat'].values))
    lons_ix = np.arange(len(ds['lon'].values))
    lats = [x.tolist() for x in np.array_split(lats_ix, 2)]
    lons = [x.tolist() for x in np.array_split(lons_ix, 2)]

    # if we have an uneven length of split arrays (srtm1 data with 3601 px)
    if len(lats[0]) != len(lats[1]):
        lats[1] = [lats[0][-1]] + lats[1]

    if len(lons[0]) != len(lons[1]):
        lons[1] = [lons[0][-1]] + lons[1]

    # split into 4 tiles [0.5x0.5 deg]
    ds1 = ds[dict(lat=lats[0], lon=lons[0])]
    ds2 = ds[dict(lat=lats[0], lon=lons[1])]
    ds3 = ds[dict(lat=lats[1], lon=lons[1])]
    ds4 = ds[dict(lat=lats[1], lon=lons[0])]
    
    return_tiles = []
    for i, ds_ in enumerate([ds1, ds2, ds3, ds4]):
        if is_empty_tile(ds_):
            log.info("Empty tile at quadrant %d detected. Ignored." % (i+1))
            return_tiles.append(None)
        else:
            return_tiles.append(ds_)
    
    return return_tiles


def get_global_attr(ds, attr_name):
    """Get the global dataset attribute."""
    if ds.attrs.has_key(attr_name):
        return ds.attrs[attr_name]
    else:
        return None


def set_global_attr(ds, attr_name, value, overwrite=False):
    """Set the global dataset attribute."""
    if ds.attrs.has_key(attr_name) and overwrite==False:
        log.error("Trying to set attr %s to %s (it already exists)." % (
                  attr_name, str(value)))
    else:
        ds.attrs[attr_name] = value
