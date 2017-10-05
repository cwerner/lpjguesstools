# -*- coding: utf-8 -*-
"""lpjguesstools._tpi: topographic position index computations according to Weiss 2001."""

# (T)opographic (P)osition (I)ndex calculations
#
# according to: Weiss, A.: Topographic position and landforms analysis, pp. 200–200. 2001.

import logging
import math
import numpy as np
import os
from scipy.ndimage.filters import generic_filter

log = logging.getLogger(__name__)

# import constants
from . import NODATA

def create_kernel(radius=2, invert=False):
    """Define a kernel"""
    
    if invert:
        value = 0
        k = np.ones((2*radius+1, 2*radius+1))                
    else:
        value = 1
        k = np.zeros((2*radius+1, 2*radius+1))
    
    y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2
    k[mask] = value
    return k


def calculate_tpi(dem, slope, scalefactor, res=30, return_unclassed=False, TYPE='SIMPLE'):
    """Classify DEM to tpi300 array according to Weiss 2001 """

    # Parameters:
    # - scalefactor: outerradius in map units (300 := 300m)
    # - res: resolution of one pixel in map units (default SRTM1: 30m)
    # - return_unclassed: return the continuous tpi values
    # - dx: cell size
    
    # inner and outer tpi300 kernels
    k_smooth = create_kernel(radius=2)
    
    radius_outer = int(math.ceil(scalefactor / float(res)))
    if radius_outer > 5:
        radius_inner  = radius_outer - 5
    else:
        log.error("Scalefactor, resolution error in tpi calc")
    
    k_outer = create_kernel(radius=radius_outer)
    k_inner = create_kernel(radius=radius_inner, invert=True)

    x = y = (k_outer.shape[0] - k_inner.shape[0]) / 2
    k_outer[x:x+k_inner.shape[0], y:y+k_inner.shape[1]] = k_inner

    # compute tpi
    tpi = dem - generic_filter(dem, np.mean, footprint=k_outer, mode="reflect") + 0.5
    tpi = generic_filter(tpi, np.mean, footprint=k_smooth, mode="reflect").astype(int)

    if TYPE == 'WEISS':
        # values from poster
        mz10, mz05, pz05, pz10 = np.percentile(tpi, [100-84.13, 100-69.15, 69.15, 84.13])
        # TODO: check if this should be a decision tree (we have unclassified cells) 
        tpi_classes = np.ones( tpi.shape ) * NODATA
        tpi_classes[(tpi > pz10)]                                   = 1 # ridge
        tpi_classes[((tpi > pz05)  & (tpi <= pz10))]                = 2 # upper slope
        tpi_classes[((tpi > mz05)  & (tpi <  pz05) & (slope >  5))] = 3 # middle slope
        tpi_classes[((tpi >= mz05) & (tpi <= pz05) & (slope <= 5))] = 4 # flats slope
        tpi_classes[((tpi >= mz10) & (tpi <  mz05))]                = 5 # lower slopes
        tpi_classes[(tpi < mz10)]                                   = 6 # valleys

    # simplified:
    if TYPE == 'SIMPLE':
        # according to Tagil & Jenness (2008) Science Alert doi:10.3923/jas.2008.910.921
        mz10, pz10 = np.percentile(tpi, [100-84.13, 84.13])
        tpi_classes = np.ones( tpi.shape ) * NODATA
        tpi_classes[(tpi >= mz10)]                               = 1 # hilltop
        tpi_classes[(tpi >= mz10) & (tpi < pz10) & (slope >= 6)] = 3 # mid slope
        tpi_classes[(tpi >  mz10) & (tpi < pz10) & (slope < 6)]  = 4 # flat surface
        tpi_classes[(tpi <= mz10)]                               = 6 # valley
    
    if return_unclassed:
        return tpi
    return tpi_classes


def classify_tpi300x2000(dem, slope):
    """Combine tpi300 and tpi2000 classification according to Weiss 2001"""

    tpi300  = calculate_tpi(dem, slope, 300, return_unclassed=True)
    tpi2000 = calculate_tpi(dem, slope, 2000, return_unclassed=True)

    tpi300_sd   = np.std( tpi300 )
    tpi300_mean = np.mean( tpi300 )

    tpi2000_sd   = np.std( tpi2000 )
    tpi2000_mean = np.mean( tpi2000 )

    tp3sd  = (((tpi300 - tpi300_mean)/tpi300_sd)*100 + 0.5).astype(int) 
    tp20sd = (((tpi2000 - tpi2000_mean)/tpi2000_sd)*100 + 0.5).astype(int) 

    lf3x20 = np.zeros( tpi2000.shape )
    lf3x20[( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd > -100) & (tp20sd < 100) & (slope <= 5))] = 5
    lf3x20[( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd > -100) & (tp20sd < 100) & (slope >  5))] = 6
    lf3x20[( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd >= 100))]                                 = 7
    lf3x20[( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd <= -100))]                                = 4
    lf3x20[( (tp3sd <= -100) & (tp20sd > -100) & (tp20sd < 100))]                                = 2
    lf3x20[( (tp3sd >= 100)  & (tp20sd > -100) & (tp20sd < 100))]                                = 9
    lf3x20[( (tp3sd <= -100) & (tp20sd >= 100))]                                                 = 3
    lf3x20[( (tp3sd <= -100) & (tp20sd <= -100))]                                                = 1
    lf3x20[( (tp3sd >= 100)  & (tp20sd >= 100))]                                                 = 10
    lf3x20[( (tp3sd >= 100)  & (tp20sd <= -100))]                                                = 8

    return lf3x20

# lookup table for Weiss landform classification
TPI300x200_LUT = {1: 'canyons, deeply incised streams', \
                  2: 'midslope drainages, shallow valleys', \
                  3: 'upland drainages, headwaters', \
                  4: 'U-shaped valleys', \
                  5: 'plains', \
                  6: 'open slopes', \
                  7: 'upper slopes, mesas', \
                  8: 'local ridges, hills in valleys', \
                  9: 'midslope ridges, small hills in plains', \
                  10: 'mountain tops, high ridges'}