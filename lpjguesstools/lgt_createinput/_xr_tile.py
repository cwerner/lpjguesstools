# -*- coding: utf-8 -*-
"""lpjguesstools._xr_tile: Tile extensions of xarray objects."""

# extensions to xarray dataset and dataarray types
# functions and properties are accessed via the tile namespace:

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
