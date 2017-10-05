# -*- coding: utf-8 -*-
"""lpjguesstools._srtm1: SRTM1 DEM specific routines."""


import numpy as np


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
            return_tiles.append(None)
        else:
            return_tiles.append(ds_)
    
    return return_tiles
