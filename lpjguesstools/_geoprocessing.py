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
import scipy
import xarray as xr

from _tpi import calculate_tpi

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
    return aspect

def calc_slope(Sx, Sy):
    """Calculate slope given X and Y slope components (unit: deg)."""
    return np.rad2deg(np.sqrt(Sx**2 + Sy**2))


def create_dem_dataset(dem, dem_mask, slope, aspect, landform, info=None):
    """Create a datasets from dem, dem_mask, slope and aspect."""
    
    # if a rasterio transfrom info is passed
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
    ds['landform'] = (xr.DataArray(landform, coords=COORDS, dims=DIMS) * ds['mask']).fillna(NODATA)
    
    for v in ['dem', 'mask', 'slope', 'aspect', 'landform']:
        ds[v].attrs.update(defaultAttrsDA)
    return ds


def compute_spatial_dataset(fname, fname_shp=None):
    """Take a GTiff file name and return a xarray datasets of dem, slope, 
    aspect and water mask layers."""

    # open source GTiff file (in WGS84)
    with rasterio.open(fname) as src:    
        msrc_kwargs = src.meta.copy()
        msrc_kwargs.update(count=5)
        msrc_kwargs.update(dtype='Float64')
        
        # create water mask from shapefile and dem
        dem = src.read(1)
        if fname_shp != None:
            with fiona.open(fname_shp) as shp:
                geoms = [feature["geometry"] for feature in shp]
                dem_mask, out_transform = rasterio.mask.mask(src, geoms, crop=False, invert=True)
                dem_mask = dem_mask.squeeze()
        else:
            dem_mask = np.zeros_like(dem)
            
        # create a in-mem copy of input dem (4 bands: dem, mask, slope, aspect)
        with rasterio.io.MemoryFile() as memfile_geo1:
            with memfile_geo1.open(**msrc_kwargs) as ds_geo1:
                print ds_geo1.shape
                ds_geo1.write(dem.astype('Float64'), 1)                 # dem
                ds_geo1.write(dem_mask.astype('Float64'), 2)            # dem_mask
                ds_geo1.write(np.zeros_like(dem_mask, 'Float64'), 3)    # slope
                ds_geo1.write(np.zeros_like(dem_mask, 'Float64'), 4)    # aspect
                ds_geo1.write(np.zeros_like(dem_mask, 'Float64'), 5)    # tpi300
                
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
                        
                        Sx, Sy = calc_slope_components(dem_filled, dx)
                        slope = calc_slope(Sx, Sy)
                        aspect = calc_aspect(Sx, Sy)
                        
                        # calculate tpi (now in utm)
                        landform = calculate_tpi(dem, slope, 300, res=dx)

                        # write slope, aspect to ds_utm
                        ds_utm.write(slope.astype('Float64'), 3)
                        ds_utm.write(aspect.astype('Float64'), 4)
                        ds_utm.write(landform.astype('Float64'), 5)
                        
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

    # create dataset    
    ds = create_dem_dataset(dem, dem_mask, slope, aspect, landform, info=msrc_kwargs)
    
    return ds

# xarray-based methods

def get_center_coord(ds):
    """Return the (lon, lat) of dataset (center)"""
    lat_c = min(ds.lat.values) + (max(ds.lat.values) - min(ds.lat.values)) * 0.5
    lon_c = min(ds.lon.values) + (max(ds.lon.values) - min(ds.lon.values)) * 0.5
    return (lon_c, lat_c)
        

def classify_aspect(ds):
    """Classify dataarray from continuous aspect to 1,2,3,4."""        
    aspect = ds['aspect'].to_masked_array()
    asp_cl = ds['aspect'].to_masked_array()
    asp_cl[(aspect >= 315) | (aspect <  45)] = 1    # North
    asp_cl[(aspect >= 45)  & (aspect < 135)] = 2    # East
    asp_cl[(aspect >= 135) & (aspect < 225)] = 3    # South
    asp_cl[(aspect >= 225) & (aspect < 315)] = 4    # West
    asp_cl = np.ma.masked_where(ds['mask'] == 0, asp_cl).filled(NODATA)
    ds['aspect_class'] = xr.full_like(ds['aspect'], NODATA)
    ds['aspect_class'][:] = asp_cl
    ds['aspect_class'].attrs.update(defaultAttrsDA)
    return ds


def classify_landform(ds, elevation_levels=[]):
    """Subdivide landform classes by aspect class."""        
    SHAPE = ds['mask'].shape
    lf_cl = np.ma.masked_array(np.ones_like(ds['mask'].values), mask=ds['mask'].values)
    
    aspect_lfs = (ds['aspect_class'].to_masked_array() > 0) & (np.in1d(ds['landform'].to_masked_array(), [2,3,5]).reshape(SHAPE))
    
    lf_cl = np.ma.where(aspect_lfs, ds['landform'] * 10 + ds['aspect_class'],
                                    ds['landform'] * 10).filled(NODATA)
    lf_cl = np.ma.masked_where(ds['mask'] == 0, lf_cl)
    
    # if we have elevation levels subdivide the landform classes
    ele = ds['dem'].to_masked_array()
    if len(elevation_levels) > 0:
        # add global elevation step attribute (second element, first is lower boundary)
        ds.attrs['landform_elevation_step'] = "%s" % elevation_levels[1]

        for i, (lb, ub) in enumerate(zip(elevation_levels[:-1], elevation_levels[1:])):
            lf_cl = np.ma.where(((ele >= lb) & (ele < ub)), lf_cl + (i+1) * 100, lf_cl)   
    
    ds['landform_class'] = xr.full_like(ds['landform'], NODATA)
    ds['landform_class'][:] = lf_cl.filled(NODATA)
    ds['landform_class'].attrs.update(defaultAttrsDA)
    
    
    return ds