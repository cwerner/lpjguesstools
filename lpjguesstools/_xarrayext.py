# -*- coding: utf-8 -*-
"""lpjguesstools._xarrayext: Extensions of xarray objects."""


# Extensions to xarray dataset and dataarray types
#
# functions are accessed via the tile namespace:
# example> ds.lgt.set('step_width', 200)

import xarray as xr


# lgt xarray extensions

@xr.register_dataset_accessor('tile')
class TileAccessor(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def set(self, attr_name, attr_value, overwrite=False):
        """Set a lgt attribute"""
        aname = 'lgt.' + attr_name
        if self._obj.attrs.has_key(aname):
            if not overwrite:
                log.warn('Attribute %s exists but overwrite is False.' % aname)
        self._obj.attrs[aname] = attr_value

    def get(self, attr_name):
        """Get a lgt attribute"""
        aname = 'lgt.' + attr_name
        if self._obj.attrs.has_key(aname):
            return self._obj.attrs[aname]
        else:
            log.warn('Attribute %s does not exist.' % aname)
            return None


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

@xr.register_dataset_accessor('geo')
class GeoAccessor(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

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


@xr.register_dataarray_accessor('geo')
class GeoAccessor(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        self._center = None

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

    @property
    def is_3d(self):
        if len(self._obj.dims) == 3:
            return True
        return False
    
