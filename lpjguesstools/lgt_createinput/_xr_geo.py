# -*- coding: utf-8 -*-
"""lpjguesstools._xr_geo: Geo extensions of xarray objects."""

# extensions to xarray dataset and dataarray types
# functions and properties are accessed via the geo namespace:

import logging
import xarray as xr

log = logging.getLogger(__name__)


# geo xarray extensions

def center(self, decimals=2):
    """Return the geographic center point of this dataset/ dataarray."""
    if self._center is None:
        # we can use a cache on our accessor objects, because accessors
        # themselves are cached on instances that access them.
        lon = self._obj.lon
        lat = self._obj.lat
        self._center = (round(float(lon.mean()), decimals), \
                        round(float(lat.mean()), decimals))
    return self._center

def contains(self, coord):
    """Check if coordinate is within extent of dataset/ dataarray."""
    
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
    """Return clipped copy of dataset/ dataarray."""
    if type(extent) == type(self._obj):
        lon1, lat1, lon2, lat2 = x.extent
    elif self.contains(x):
        lon1, lat1, lon2, lat2 = x
    else:
        log.warn("Clip failed.")
        return None
    return self._obj.sel(lon=(self._obj.lon >= lon1 | self._obj.lon <= lon2),
                         lat=(self._obj.lat >= lat2 | self._obj.lat <= lat2))

# properties
@property
def extent(self):
    lon, lat = self._obj.lon, self._obj.lat
    extent = [min(lon.values), min(lat.values), max(lon.values), max(lat.values)]
    if self._extent is None:
        self._extent = extent
    return self._extent

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
        self._extent = None

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
        self._extent = None

    center = center
    clip = clip
    contains = contains
    extent = extent    
    is_3d = is_3d
