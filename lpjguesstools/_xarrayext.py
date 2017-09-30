# -*- coding: utf-8 -*-
"""lpjguesstools._xarrayext: Extensions of xarray objects."""


# Extensions to xarray dataset and dataarray types
#
# functions are accessed via the tile namespace:
# example> ds.lgt.set('step_width', 200)

import logging
import xarray as xr

log = logging.getLogger(__name__)


# lgt xarray extensions

@xr.register_dataset_accessor('tile')
class TileAccessor(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self.lgtattrs = set()
        
    def _get_lgt_attrs(self):
        # parse any existing lgt attributes
        lgtattrs = [x.replace('lgt.','') for x in self._obj.attrs.keys() if 'lgt.' in x]
        self.lgtattrs = set(lgtattrs)        

    def set(self, attr_name, attr_value, overwrite=False):
        """Set a lgt attribute"""
        aname = 'lgt.' + attr_name
        if attr_name in self.lgtattrs:
            if not overwrite:
                log.warn('LGT attribute "%s" exists but overwrite is False.' % aname)
                return    
        self.lgtattrs.union(set(attr_name))
        self._obj.attrs[aname] = attr_value

    def get(self, attr_name):
        """Get a lgt attribute"""
        if len(self.lgtattrs) == 0:
            self._get_lgt_attrs()
        if attr_name in self.lgtattrs:
            return self._obj.attrs['lgt.' + attr_name]
        else:
            log.warn('LGT attribute "%s" does not exist in dataset.' % attr_name)
            return None
    
    def copy_attrs(self, src, overwrite=False):
        """Copy lgt attributes from src to dst dataset."""
        # check src.tile.lgtattrs
        if len(src.tile.lgtattrs) == 0:
            src.tile._get_lgt_attrs()
        for attr_name in src.tile.lgtattrs:
            self.set(attr_name, src.tile.get(attr_name), overwrite=overwrite)

            
@xr.register_dataarray_accessor('tile')
class TileAccessor(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def update_attrs(self, *args, **kwargs):
        """Update the attributes in a xarray DataArray"""
        def _update_attrs(obj, *args, **kwargs):
            obj.attrs.update(*args, **kwargs)
            return obj
        self._obj.pipe(_update_attrs, *args, **kwargs)

    def update_encoding(self, *args, **kwargs):
        """Update the encoding in a xarray DataArray"""
        def _update_encoding(obj, *args, **kwargs):
            obj.encoding.update(*args, **kwargs)
            return obj
        self._obj.pipe(_update_encoding, *args, **kwargs)


# geo xarray extensions

def center(self, decimals=2):
    """Return the geographic center point of this tile."""
    if self._center is None:
        # we can use a cache on our accessor objects, because accessors
        # themselves are cached on instances that access them.
        lon = self._obj.lon
        lat = self._obj.lat
        self._center = (round(float(lon.mean()), decimals), \
                        round(float(lat.mean()), decimals))
    return self._center

def contains(self, coord):
    """Check if coordinate is within extent of dataset"""
    
    # allow coord tuple and 2-coord extent
    if (type(coord) in [list, tuple]) and (len(coord) == 4):
        coords = [(coord[0], coord[1]), (coord[2], coord[3])]
    elif (type(coord) in [list, tuple]) and (len(coord) == 2):
        coords = [coord]
    else:
        log.warn('Invalid extent: %s' % str(coords))
        return False
    
    lon, lat = self._obj.lon, self._obj.lat
    checks = []
    for c in coords:
        inside = False
        if (c[0] >= lon.min() and c[0] <= lon.max()):
            if (c[1] >= lat.min() and c[1] <= lat.max()):
                inside = True
        checks.append(inside)
    if all(checks):
        return True
    return False

def clip(self, extent):
    """Return clipped copy of dataset."""
    if type(extent) == type(self._obj):
        lon1, lat1, lon2, lat2 = x.extent
    elif self.contains(x):
        lon1, lat1, lon2, lat2 = x
    else:
        log.warn("Clip failed.")
        return None
    return self._obj.sel(lon=(self._obj.lon >= lon1 | self._obj.lon <= lon2),
                         lat=(self._obj.lat >= lat2 | self._obj.lat <= lat2))

@property
def extent(self):
    lon, lat = self._obj.lon, self._obj.lat
    extent = [min(lon), min(lat), max(lon), max(lat)]
    if self._extent != None:
        self._extent = extent
    return extent

@property
def is_3d(self):
    if len(self._obj.dims) == 3:
        return True
    return False


@xr.register_dataset_accessor('geo')
class GeoAccessor(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

    center = center
    clip = clip
    contains = contains
    extent = extent    
    is_3d = is_3d


@xr.register_dataarray_accessor('geo')
class GeoAccessor(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

    center = center
    clip = clip
    contains = contains
    extent = extent    
    is_3d = is_3d