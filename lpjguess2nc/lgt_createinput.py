"""FILE lgt_createinput.py.

This script creates condensed LPJ netcdf files
for landforms and soil properties

 landforms.nc:
 - lfcnt (landid) number of landforms in cell
 - frac (landid, lfid/ standid) area fraction this landform represents
 - slope (landid, lfid/ standid)
 - elevation (landid, lfid/ standid) avg. elevation in this landform
 - soildepth (landid, lfid/ standid) [implemented later const in model for now]

 sites.nc:
 - soildepth
 - clay
 - silt
 - sand
 - totc
 - elevation (reference elevation for grid, 0.5deg)

Christian Werner, SENCKENBERG Biodiversity and Climate Research Centre (BiK-F)
email: christian.werner@senkenberg.de
2017/02/07
"""

from collections import Counter, OrderedDict
import copy
import datetime
import gdal
import gdalconst
import glob 
import math
import numpy as np
import os
from osgeo import osr, gdal, gdalconst
import pandas as pd
import rasterio as rio
from scipy.ndimage.filters import generic_filter as gf
import string
import xarray as xr


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
    """rasterize mask shapefile (water bodies)"""

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
    return(dest_ds)


def processDEM(src, verbose=False, force=False, mask=None, GDALResampleAlg=gdal.GRA_Bilinear):
    """Calculate slope, aspect and tpi from source DEM"""
    
    ## test if 'src' is a valid GDAL
    if not isinstance(src, gdal.Dataset):
        print("ERROR (project2UTM): 'src' not of type 'gdal.Dataset'!")
        sys.exit(-999)

    ## check if a mask layer is desired
    if mask == None:
        nLayer = 3
    else:
        nLayer = 4

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

    ## put input/slope/aspect into one 3-Band GDAL object
    dst = mem_drv.Create('', src.RasterXSize,  src.RasterYSize, 3, gdal.GDT_Float64)
    dst.SetGeoTransform(src.GetGeoTransform())
    dst.SetProjection(src.GetProjection())
    
    dst.GetRasterBand(1).WriteArray(src.GetRasterBand(1).ReadAsArray())
    dst.GetRasterBand(2).WriteArray(dst_slp.GetRasterBand(1).ReadAsArray())
    dst.GetRasterBand(3).WriteArray(dst_asp.GetRasterBand(1).ReadAsArray())

    if not mask == None:
        rmask = rasterizeMask(mask, src)
        dst.GetRasterBand(4).WriteArray(rmask.GetRasterBand(1).ReadAsArray())

    return(dst)



    
def compute_tpi(geofile, scale=300):
    """compute topograition index based on DEM""" 
    
    # scale in meters: tested 300, 2000 based on Weiss Poster

    return tpi



def call_custom_tpi():

    #tiles = glob.glob('SRTM1/*_bil.zip')
    
    # Rules:
    # - slope and aspect calc should be conducted on water-clipped grid
    # - tpi calculations should be carried out on original data (since they)
    #   involve a potentially larger kernel (otherwise we have NODATA increase)
    #
    # TODO:
    # - include a new step to clip srtm1_filled tifs with their shp files
    #   this was formerly done offline using GRASS 
    #
    # -----8< --------------------------------------------------
    # !/bin/bash
    # tifdir=data/srtm1_filled
    # shpdir=data/srtm1_water
    # outdir=data/srtm1_masked
    # 
    # GISDBASE=/data/Grass
    # 
    # LOCATION_NAME=Chile
    # MAPSET=PERMANENT
    # 
    # for file in $tifdir/*.tif; do
    #   file=$(basename $file)
    #   lat=$(echo $file | grep -Po '[ns][0-9]{2}')
    #   lon=$(echo $file | grep -Po '[ew][0-9]{3}')
    #   echo $lon $lat
    # 
    #   grass -c ${tifdir}/${lat}_${lon}_1arc_v3.tif -e $GISDBASE/$LOCATION_NAME
    #   grass ${GISDBASE}/${LOCATION_NAME}/${MAPSET}/ --exec r.import --quiet input=${tifdir}/${lat}_${lon}_1arc_v3.tif output=dem
    #   grass ${GISDBASE}/${LOCATION_NAME}/${MAPSET}/ --exec v.import --quiet -o input=${shpdir}/${lon}${lat}s.shp output=water
    #   grass ${GISDBASE}/${LOCATION_NAME}/${MAPSET}/ --exec r.mask --quiet -i vector=water
    #   grass ${GISDBASE}/${LOCATION_NAME}/${MAPSET}/ --exec r.out.gdal -mc --quiet --overwrite input=dem output=${outdir}/${lat}_${lon}_1arc_v3.tif type=Int16 nodata=-32768
    #   rm -r $GISDBASE/$LOCATION_NAME
    # done
    # -----8< --------------------------------------------------
    
    # 
    tiles = glob.glob('srtm1_filled/*.tif')

    for tile in sorted( tiles ):
        print 'processing: ', tile, '(', datetime.datetime.now(), ')'
        #tile = "SRTM1/s17_w071_1arc_v3_bil.zip"

        #zfilename = os.path.basename(tile)
        #bfilename = zfilename[:-8] + '.bil'

        #src = gdal.Open('/vsizip/' + os.path.join(tile, bfilename), gdalconst.GA_ReadOnly)

        src = gdal.Open(tile, gdalconst.GA_ReadOnly)
        rawdata = processDEM(src)

        DEM       = src.GetRasterBand(1).ReadAsArray()
        SLOPE     = src.GetRasterBand(2).ReadAsArray()
        ASPECT    = src.GetRasterBand(3).ReadAsArray()
        WATERMASK = src.GetRasterBand(4).ReadAsArray()
        
        NODATA = src.GetRasterBand(1).GetNoDataValue()
        if NODATA is not None:
            DEM    = np.ma.masked_equal(DEM, NODATA)
            SLOPE  = np.ma.masked_equal(SLOPE, NODATA)
            ASPECT = np.ma.masked_equal(ASPECT, NODATA)
            DEMM   = np.ma.masked_where(WATERMASK == 1, DEM) 


        # tpi 300 first / small scale
        print '   - tpi300 calculation'
        # kernel

        # radii are for for 30m DEM

        # smoothing circle (5x5 window)
        radius = 2
        k1 = np.zeros((2*radius+1, 2*radius+1))
        y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
        mask = x**2 + y**2 <= radius**2
        k1[mask] = 1

        # apply kernels to DATA
        radius = 10
        kernel = np.zeros((2*radius+1, 2*radius+1))
        y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
        mask = x**2 + y**2 <= radius**2
        kernel[mask] = 1

        # cut out inner area again
        sradius = 5
        skernel = np.ones((2*sradius+1, 2*sradius+1))
        y,x = np.ogrid[-sradius:sradius+1, -sradius:sradius+1]
        smask = x**2 + y**2 <= sradius**2
        skernel[smask] = 0



        x = y = (kernel.shape[0] - skernel.shape[0]) / 2
        kernel[x:x+skernel.shape[0], y:y+skernel.shape[1]] = skernel

        tpi300 = DEM - gf(DEM, np.mean, footprint=kernel, mode="reflect") + 0.5
        tpi300 = gf(tpi300, np.mean, footprint=k1, mode="reflect").astype(int)


        # stats
        tpi300_sd   = np.std( tpi300 )
        tpi300_mean = np.mean( tpi300 )

        pz05_val = np.percentile(tpi300, 69.15)
        pz10_val = np.percentile(tpi300, 84.13)
        mz05_val = np.percentile(tpi300, 100 - 69.15)
        mz10_val = np.percentile(tpi300, 100 - 84.13)

        # reclassification
        tpi300_classes = np.zeros( tpi300.shape )
        tpi300_classes[np.where( tpi300 > pz10_val)]                                           = 1 # ridge
        tpi300_classes[np.where( (tpi300 > pz05_val)  & (tpi300 <= pz10_val)) ]                = 2 # upper slope
        tpi300_classes[np.where( (tpi300 > mz05_val)  & (tpi300 < pz05_val) & (SLOPE > 6) ) ]  = 3 # middle slope
        tpi300_classes[np.where( (tpi300 >= mz05_val) & (tpi300 <= pz05_val) & (SLOPE <= 6)) ] = 4 # flats slope
        tpi300_classes[np.where( (tpi300 >= mz10_val) & (tpi300 <  mz05_val)) ]                = 5 # lower slopes
        tpi300_classes[np.where( tpi300 < mz10_val)]                                           = 6 # valleys


        print '   - tpi2000 calculation'

        # apply kernels to DATA
        radius = 67

        kernel = np.zeros((2*radius+1, 2*radius+1))
        y,x = np.ogrid[-radius:radius+1, -radius:radius+1]
        mask = x**2 + y**2 <= radius**2
        kernel[mask] = 1

        # cut out inner area again
        sradius = 62
        skernel = np.ones((2*sradius+1, 2*sradius+1))
        y,x = np.ogrid[-sradius:sradius+1, -sradius:sradius+1]
        smask = x**2 + y**2 <= sradius**2
        skernel[smask] = 0

        x = y = (kernel.shape[0] - skernel.shape[0]) / 2
        kernel[x:x+skernel.shape[0], y:y+skernel.shape[1]] = skernel

        tpi2000 = DEM - gf(DEM, np.mean, footprint=kernel, mode="reflect") + 0.5
        tpi2000 = gf(tpi2000, np.mean, footprint=k1, mode="reflect").astype(int)

        tpi2000_sd   = np.std( tpi2000 )
        tpi2000_mean = np.mean( tpi2000 )


        pz05_val = np.percentile(tpi2000, 69.15)
        pz10_val = np.percentile(tpi2000, 84.13)
        mz05_val = np.percentile(tpi2000, 100 - 69.15)
        mz10_val = np.percentile(tpi2000, 100 - 84.13)


        # reclassification
        tpi2000_classes = np.zeros( tpi2000.shape )
        tpi2000_classes[np.where( tpi2000 > pz10_val)]                                            = 1 # ridge
        tpi2000_classes[np.where( (tpi2000 > pz05_val)  & (tpi2000 <= pz10_val)) ]                = 2 # upper slope
        tpi2000_classes[np.where( (tpi2000 > mz05_val)  & (tpi2000 < pz05_val) & (SLOPE > 5) ) ]  = 3 # middle slope
        tpi2000_classes[np.where( (tpi2000 >= mz05_val) & (tpi2000 <= pz05_val) & (SLOPE <= 5)) ] = 4 # flats slope
        tpi2000_classes[np.where( (tpi2000 >= mz10_val) & (tpi2000 <  mz05_val)) ]                = 5 # lower slopes
        tpi2000_classes[np.where( tpi2000 < mz10_val)]                                            = 6 # valleys


        # join classsifications
        # -----------------------------------------------------------------

        print '   - landform classification calculation'

        tp3sd  = (((tpi300 - tpi300_mean)/tpi300_sd)*100 + 0.5).astype(int) 
        tp20sd = (((tpi2000 - tpi2000_mean)/tpi2000_sd)*100 + 0.5).astype(int) 

        lf3x20 = np.zeros( tpi2000.shape )
        lf3x20[np.where( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd > -100) & (tp20sd < 100) & (SLOPE <= 5))] = 5
        lf3x20[np.where( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd > -100) & (tp20sd < 100) & (SLOPE >  5))] = 6
        lf3x20[np.where( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd >= 100))]                                 = 7
        lf3x20[np.where( (tp3sd > -100)  & (tp3sd < 100) & (tp20sd <= -100))]                                = 4
        lf3x20[np.where( (tp3sd <= -100) & (tp20sd > -100) & (tp20sd < 100))]                                = 2
        lf3x20[np.where( (tp3sd >= 100)  & (tp20sd > -100) & (tp20sd < 100))]                                = 9
        lf3x20[np.where( (tp3sd <= -100) & (tp20sd >= 100))]                                                 = 3
        lf3x20[np.where( (tp3sd <= -100) & (tp20sd <= -100))]                                                = 1
        lf3x20[np.where( (tp3sd >= 100)  & (tp20sd >= 100))]                                                 = 10
        lf3x20[np.where( (tp3sd >= 100)  & (tp20sd <= -100))]                                                = 8

        # return arrays
        # TODO: convert to an xarray dataset
        
        return (DEM, SLOPE, ASPECT, tpi300, tpi300_classes, tpi2000, tpi2000_classes, lf3x20)


# ----------------- PART2 -------------------------------
# (formerly known as computeLandformStats5.py)
#
# Computes LandForm stats as TXT files for all pixels
#
# Input files:
# - srtm1_masked [gap-filled and water-clipped elevation data]
# - tpi300NEW, slopeNEW, aspectNEW
#
# Files produced:
# - netcdfs_tpi300 [0.5deg netcdf files with aspect, slope, elevation, landform]
# - elevations_full_tpi300.txt
# - landforms_full_tpi300.txt
# - slopes_full_tpi300.txt
#
# NOTES:
# - txt data is filtered for 1% area (smaller units -1)
# - netcdfs are original, for spatial mapping of filtered data onto netcdfs
#   run script 'create_lf_avg_netcdfs.py'
#
# Christian Werner
# christian.werner@senckenerg.de
# 2017/03



LUT = {1: 'canyons, deeply incised streams', \
       2: 'midslope drainages, shallow valleys', \
       3: 'upland drainages, headwaters', \
       4: 'U-shaped valleys', \
       5: 'plains', \
       6: 'open slopes', \
       7: 'upper slopes, mesas', \
       8: 'local ridges, hills in valleys', \
       9: 'midslope ridges, small hills in plains', \
       10: 'mountain tops, high ridges'}


    # 1: ridge
    # 2: upper slope
    # 3: midslope
    # 4: flats
    # 5: lower slope
    # 6: valleys


def rename_lat_lon(lat, lon):
    latS = 's'
    lonS = 'w'
    
    if lat >= 0: latS = 'n'
    if lon >= 0: lonS = 'e'
    
    outS = '%s%2.2f%s%2.2f_0.5arc_%s.nc' % (latS, math.fabs(lat), lonS, math.fabs(lon), TYPE)
    return outS

def call_computeLandformStats5():
    # call old computeLandformStats5.py code
    

    CUTOFF = 1.0        # min percent landmass for LF to be considered
    TYPE = 'tpi300'     # tpi2000

    #for tile in sorted( glob.glob("landforms/*.tif") ):

    f_el = open('elevations_full_%s.txt' % TYPE, 'w')
    f_lf = open('landforms_full_%s.txt' % TYPE, 'w')
    f_sl = open('slopes_full_%s.txt' % TYPE, 'w')

    # landforms (TPI300x2000)
    #header = '\t'.join(['lat', 'lon', 'lf_cnt'] + ['LF%d' % x for x in range(1,11)]) + '\n'

    # elevation boundaries
    ele_breaks = [-1000] + range(400, 4801, 400) + [10000]
    ele_cnt = range(1, len(ele_breaks))

    lf_set = [10,21,22,23,24,31,32,33,34,40,51,52,53,54,60]

    lf_full_set = []
    for e in ele_cnt:
        lf_full_set += [x+(100*e) for x in lf_set]

    # new landforms
    header = '\t'.join(['lat', 'lon', 'lf_cnt'] + ['LF%d' % x for x in lf_full_set]) + '\n'

    f_el.write(header)
    f_lf.write(header)
    f_sl.write(header)

    tiles = glob.glob('/Volumes/Drobo/projects/EarthShape/subpixel/srtm1_masked/*.tif')

    for tile in sorted(tiles):

        # get lower left coord of 1 deg tile, all in the west/south hemisphere
        b = [-int(x[1:]) for x in os.path.basename(tile).split('_')[0:2]]
        blat, blon = b

        print blat, blon
        
        #if blat > -54:
        #    continue

        # calc slope file
        #os.system('gdaldem slope -s 111122 %s %s' % (tile, string.replace(tile[:-4], 'srtm1', 'slope') + '_slope.tif'))

        #print tile, 
        #src = gdal.Open( string.replace(tile[:-4], 'srtm1', 'landforms') + '_lf300x2000.tif' , gdalconst.GA_ReadOnly)
        
        # landforms
        
        src = gdal.Open( string.replace(tile[:-4], 'srtm1_masked', TYPE + 'NEW') + '_%ssm_classed.tif' % TYPE, gdalconst.GA_ReadOnly)
        LF      = src.GetRasterBand(1).ReadAsArray()
        LNODATA = src.GetRasterBand(1).GetNoDataValue()

        # slopes
        src = gdal.Open( string.replace(tile[:-4], 'srtm1_masked', 'slopeNEW') + '_slope.tif', gdalconst.GA_ReadOnly)
        SLOPE   = src.GetRasterBand(1).ReadAsArray()
        SNODATA = src.GetRasterBand(1).GetNoDataValue()

        # aspects
        src = gdal.Open( string.replace(tile[:-4], 'srtm1_masked', 'aspectNEW') + '_aspect.tif', gdalconst.GA_ReadOnly)
        ASPECT   = src.GetRasterBand(1).ReadAsArray()
        ANODATA = src.GetRasterBand(1).GetNoDataValue()
        
        # classify aspect
        ASPECT2 = np.ma.masked_where(np.asarray(ASPECT, 'i') < 0, np.asarray(ASPECT, 'i'))
        ASPECT2[np.where(( ASPECT >= 315 ) | ((ASPECT>=0) & (ASPECT < 45)))] = 1   # north
        ASPECT2[np.where(( ASPECT >= 45  ) & ( ASPECT < 135))] = 2   # east
        ASPECT2[np.where(( ASPECT >= 135 ) & ( ASPECT < 225))] = 3   # south
        ASPECT2[np.where(( ASPECT >= 225 ) & ( ASPECT < 315))] = 4   # west

        # elevation
        src     = gdal.Open( tile, gdalconst.GA_ReadOnly)
        DEM     = src.GetRasterBand(1).ReadAsArray()
        DNODATA = src.GetRasterBand(1).GetNoDataValue()

        
        # the one and only base mask !!!
        DEMMASK = np.ma.masked_where( (DEM == DNODATA), np.ones_like(DEM))

        # secondary masks (are merged with the water mask of srt1_masked)
        ASPMASK = np.where(ASPECT2 < 0, 1, 0) * np.where(DEMMASK.filled(-9999) == -9999, 1, 0)
        SLOMASK = np.where(SLOPE < 0, 1, 0) # * np.where(DEMMASK.filled(-9999) == -9999, 1, 0)
        LFMASK  = np.where(LF < 0, 1, 0) * np.where(DEMMASK.filled(-9999) == -9999, 1, 0)

        # mask with secondary masks
        DEM0 = np.ma.masked_where( DEM == DNODATA, DEM)
        LF0  = np.ma.masked_where( LFMASK == 1, LF)
        SLO0 = np.ma.masked_where( SLOMASK == 1, SLOPE)
        ASP0 = np.ma.masked_where( ASPMASK == 1, ASPECT2)

        # calculate new landform code
        LF0B    = np.ma.where(((ASP0.filled(-1) > 0) & (np.in1d(LF0, [2,3,5]).reshape(LF0.shape))), (LF0 * 10) + ASP0.filled(-1), LF0 * 10)

        # split into 1x1 into 0.5x0.5 tiles

        for t in range(4):

            # reset output lists
            slopes     = []
            elevations = []


            #print 'subset %d -------' % i
            if t == 0:
                print LF.shape
                LF    = LF0[0:1801,0:1801]
                SLOPE = SLO0[0:1801,0:1801]
                DEM   = DEM0[0:1801,0:1801]
                LFB   = LF0B[0:1801,0:1801]
                AS    = ASP0[0:1801,0:1801]

                print LFB.shape

                lat = blat + 0.75; lon =  blon + 0.25
            if t == 1: 
                LF    = LF0[0:1801,1800:]
                SLOPE = SLO0[0:1801,1800:]
                DEM   = DEM0[0:1801,1800:]
                LFB   = LF0B[0:1801,1800:]
                AS    = ASP0[0:1801,1800:]

                print LFB.shape

                lat = blat + 0.75; lon =  blon + 0.75
            if t == 2: 
                LF    = LF0[1800:,0:1801]
                SLOPE = SLO0[1800:,0:1801]
                DEM   = DEM0[1800:,0:1801]
                LFB   = LF0B[1800:,0:1801]
                AS    = ASP0[1800:,0:1801]

                print LFB.shape

                lat = blat + 0.25; lon =  blon + 0.25
            if t == 3: 
                LF    = LF0[1800:,1800:]
                SLOPE = SLO0[1800:,1800:]
                DEM   = DEM0[1800:,1800:]
                LFB   = LF0B[1800:,1800:]
                AS    = ASP0[1800:,1800:]

                print LFB.shape

                lat = blat + 0.25; lon =  blon + 0.75

            # add dem mask to all arrays
            _mask = np.ma.getmask(DEM)
            if np.ndim(_mask) != 0:
                SLOPE[_mask] = np.ma.masked
                LFB[_mask] = np.ma.masked
                AS[_mask] = np.ma.masked
                LF[_mask] = np.ma.masked

            # create netcdf
            ds = xr.Dataset()
            foutname = 'netcdfs_%s/' % TYPE + rename_lat_lon(lat, lon)
            
            Dunits = {'slope': 'degree', 'elevation': 'm', 'landform': '-', 'aspect': '-'}
            Dlname = {'slope': 'Slope',
                      'elevation': 'Elevation a.s.l.',
                      'landform': 'Landfrom unit',
                      'aspect': 'Aspect'}

            LFC = copy.deepcopy(LFB)

            LATS = np.linspace(0.25,-0.25,1801) + lat
            LONS = np.linspace(-0.25,0.25,1801) + lon

            # ---
            for el_cnt, (lower, upper) in enumerate(zip(ele_breaks[0:-1], ele_breaks[1:])):
                el_mask = np.ma.where((DEM >= lower) & (DEM < upper ), True, False)

                #print el_cnt, np.ma.sum(el_mask)
                LFC = np.ma.where(el_mask == True, LFB + (el_cnt+1)*100, LFC)

            # dump to netcdf
            for name, da in zip(['aspect', 'slope', 'elevation', 'landform'],
                                [AS, SLOPE, DEM, LFC]):

                # extra mask apply (make sure we do not get dodgy slopes
                
                #_mask = np.ma.getmask(DEM)
                #if np.ndim(_mask) != 0:
                #    da[_mask] = np.ma.masked

                da = xr.DataArray(da[:].filled(-9999),
                                  name=name,
                                  coords=[('lat', LATS[:]), ('lon', LONS[:])])
                da.attrs['units'] = Dunits[name]
                da.attrs['long_name'] = Dlname[name]
                da.attrs['missing_value'] = -9999
                da.attrs['_FillValue'] = -9999
                ds[name] = da
            
            ds.to_netcdf(foutname, format='NETCDF4_CLASSIC')


            # calc avg. slope for landform
            if DEM.mask.all():
                pass
            else:

                a = LFC.filled(-1).ravel()
                b = Counter(a)
                bins = [(x, b[x]) for x in lf_full_set]
                bins_sum = sum([n for _, n in bins if n > 0])

                def check_cutoff(x, _sum, _cut):
                    """ check if landform akes the cutoff """
                    if x == 0:
                        return False
                    else:
                        if (100 * float(x) / _sum >= _cut):
                            return True
                        else:
                            return False


                def perc(x, _sum):
                    """ calculate percentage covered"""
                    if x <= 0:
                        return -1
                    else:
                        return round((x/float(_sum))*100.0, 3)


                bins_cut  = [(_, n) if check_cutoff(n, bins_sum, CUTOFF) else (_, 0) for _, n in bins]
                bins_sum2 =  sum([n for _, n in bins_cut if n > 0])
                bins_frac = [perc(n, bins_sum2) if n > -1 else -1 for _, n in bins_cut]

                for el_cnt, (lower, upper) in enumerate(zip(ele_breaks[0:-1], ele_breaks[1:])):
                    el_mask = np.ma.where((DEM >= lower) & (DEM < upper ), True, False)

                    for lf_id in lf_set:
                        lf_mask = np.ma.where(LFB == lf_id, True, False)
                        full_mask = np.invert(el_mask * lf_mask)

                        _SLOPE = copy.deepcopy(SLOPE)
                        _SLOPE[full_mask] = np.ma.masked

                        _DEM = copy.deepcopy(DEM)
                        _DEM[full_mask] = np.ma.masked

                        # iterate over elevation bands (400m each)
                        sl = np.ma.mean(_SLOPE)
                        el = np.ma.mean(_DEM)

                        if type(sl) is np.ma.core.MaskedConstant:
                            sl = -1
                        slopes.append(sl)

                        if type(el) is np.ma.core.MaskedConstant:
                            el = -1
                        elevations.append(el)

                # write data files 

                valid_lfs = sum([1 for x in bins_frac if x > -1])
                out_frac  = [lat, lon, valid_lfs] + bins_frac
                out_slope = [lat, lon, valid_lfs] + [round(x, 2) if (x != -1) and (y != -1) else -1 for x, y in zip(slopes, bins_frac)]
                out_ele   = [lat, lon, valid_lfs] + [int(round(x, 0)) if (x != -1) and (y != -1) else -1 for x, y in zip(elevations, bins_frac)]

                s_out_fr = '\t'.join([str(x) for x in out_frac]) + '\n'
                s_out_sl = '\t'.join([str(x) for x in out_slope]) + '\n'
                s_out_el = '\t'.join([str(x) for x in out_ele]) + '\n'

                f_el.write( s_out_el )
                f_lf.write( s_out_fr )
                f_sl.write( s_out_sl )
                
                f_el.flush()
                f_lf.flush()
                f_sl.flush()

f_el.close()
f_lf.close()
f_sl.close()



# ----------------- PART3 -------------------------------
# (formerly known as create_lpj_lf_netcdf.py)
#
# now base for lgt_createinput.py

# consts and lookups
NODATA = -9999

defaultD = {'missing_value': NODATA, '_FillValue': NODATA}

varD = {'TOTC': ('SOC', 'Soil Organic Carbon', 'percent', 0.1),
        'SDTO': ('SAND', 'Sand', 'percent', 1.0),
        'STPC': ('SILT', 'Silt', 'percent', 1.0),
        'CLPC': ('CLAY', 'Clay', 'percent', 1.0)}

soil_vars = sorted(varD.keys())

# define landform ids
ele_breaks = [-1000] + range(400, 4801, 400) + [10000]
ele_cnt = range(1, len(ele_breaks))
lf_set = "10,21,22,23,24,31,32,33,34,40,51,52,53,54,60".split(',')
lf_full_set = []
for e in ele_cnt:
    lf_full_set += [int(x)+(100*e) for x in lf_set]


def is_3d(ds, v):
    """Check if xr.DataArray has 3 dimensions."""
    dims = ds[v].dims
    if len(dims) == 3:
        return True
    return False


def assign_to_dataarray(data, df, refdata=False):
    """Place value into correct location of data array."""
    for _, r in df.iterrows():
        if refdata:
            data.loc[r.lat, r.lon] = r.lf_cnt
        else:
            v = np.array(r.tolist()[3:])
            data.loc[lf_full_set, r.lat, r.lon] = np.where(v == -1, NODATA, v)

    return data


def build_site_netcdf(soilref, elevref):
    """Build the site netcdf file."""
    ds_soil = xr.open_dataset(soilref)
    ds_ele = xr.open_dataset(elevref)

    lat_slice = slice(-56, -16)
    lon_slice = slice(-76, -66)

    # slice simulation domain
    ds_soil_cl = ds_soil.sel(lat=lat_slice, lon=lon_slice,
                             lev=1.0).squeeze(drop=True)

    del ds_soil_cl['lev']

    ds_ele_cl = ds_ele.sel(latitude=lat_slice,
                           longitude=lon_slice).squeeze(drop=True)

    # identify locations that need filling and use left neighbor
    smask = np.where(ds_soil_cl['TOTC'].to_masked_array().mask, 1, 0)
    emask = np.where(ds_ele_cl['data'].to_masked_array().mask, 1, 0)
    missing = np.where((smask == 1) & (emask == 0), 1, 0)

    ix, jx = np.where(missing == 1)
    for i, j in zip(ix, jx):
        for v in soil_vars:
            ds_soil_cl[v][i, j] = ds_soil_cl[v][i, j-1]

    dsout = xr.Dataset()
    # soil vars
    for v in soil_vars:
        conv = varD[v][-1]
        da = ds_soil_cl[v].copy(deep=True) * conv
        da.name = varD[v][0]

        vattr = {'name': varD[v][0],
                 'long_name': varD[v][1],
                 'units': varD[v][2]}
        vattr.update(defaultD)
        da.attrs.update(vattr)
        da[:] = np.ma.masked_where(emask, da.to_masked_array().filled(NODATA))
        dsout[da.name] = da

    # ele var
    da = xr.full_like(da.copy(deep=True), NODATA)
    da.name = 'ELEVATION'
    vattr = {'name': 'elevation', 'long_name': 'Elevation', 'units': 'meters'}
    vattr.update(defaultD)
    da.attrs.update(vattr)
    
    da[:] = ds_ele_cl['data'].to_masked_array().filled(NODATA)
    dsout[da.name] = da

    return dsout


def build_landform_netcdf(refnc=None):
    """Build landform netcdf based on refnc dims."""
    dsout = xr.Dataset()

    COORDS = [('lf_id', lf_full_set), ('lat', refnc.lat), ('lon', refnc.lon)]
    SHAPE = tuple([len(x) for _, x in COORDS])

    # initiate data arrays
    _blank = np.ones(SHAPE) * NODATA
    da_lfcnt = xr.DataArray(_blank[0,:,:].astype('i'), name='lfcnt', 
                            coords=COORDS[1:])
    da_frac = xr.DataArray(_blank[:], name='frac', coords=COORDS)
    da_slope = xr.DataArray(_blank[:], name='slope', coords=COORDS)
    da_elev = xr.DataArray(_blank[:], name='elevation', coords=COORDS)

    # parse txt input files
    df = pd.read_csv('landforms_full_tpi300.txt', delim_whitespace=True)
    da_lfcnt = assign_to_dataarray(da_lfcnt, df, refdata=True)
    da_frac = assign_to_dataarray(da_frac, df)

    df = pd.read_csv('slopes_full_tpi300.txt', delim_whitespace=True)
    da_slope = assign_to_dataarray(da_slope, df)

    df = pd.read_csv('elevations_full_tpi300.txt', delim_whitespace=True)
    da_elev = assign_to_dataarray(da_elev, df)

    # store arrays in dataset
    dsout[da_lfcnt.name] = da_lfcnt
    dsout[da_frac.name] = da_frac
    dsout[da_slope.name] = da_slope
    dsout[da_elev.name] = da_elev
    for dv in dsout.data_vars:
        dsout[dv].attrs.update(defaultD)

    return dsout


def build_compressed(ds):
    """Build LPJ-Guess 4.0 compatible compressed netcdf file."""
    # identify landforms netcdf
    if 'lfcnt' in ds.data_vars:
        v = 'lfcnt'
    elif 'ELEVATION' in ds.data_vars:
        v = 'ELEVATION'
    else:
        print "Not a valid xr.Dataset (landforms or site only)."
        exit()

    # create id position dataarray
    da_ids = xr.ones_like(ds[v]) * NODATA 

    latL = []
    lonL = []
    d = ds[v].to_masked_array()
    
    # REVIEW: why is 'to_masked_array()'' not working here?
    d = np.ma.masked_where(d == NODATA, d)

    land_id = 0
    D_ids = OrderedDict()
    for j in reversed(range(len(d))):
        for i in range(len(d[0])):
            if d[j, i] is not np.ma.masked:
                lat = float(ds['lat'][j].values)
                lon = float(ds['lon'][i].values)
                latL.append(lat)
                lonL.append(lon)
                da_ids.loc[lat, lon] = land_id
                D_ids[(lat, lon)] = land_id
                land_id += 1 
    LFIDS = range(land_id)    
    
    # create coordinate variables
    _blank = np.zeros(len(LFIDS))
    lats = xr.DataArray(latL, name='lat', coords=[('land_id', LFIDS)])
    lons = xr.DataArray(lonL, name='lon', coords=[('land_id', LFIDS)])

    # create land_id reference array
    # TODO: clip land_id array to Chile country extent?
    da_ids.attrs.update(defaultD)
    ds_ids = da_ids.to_dataset(name='land_id')

    # create xr.Dataset
    dsout = xr.Dataset()
    dsout[lats.name] = lats
    dsout[lons.name] = lons

    # walk through variables, get lat/ lon cells' data 
    for v in ds.data_vars:
        if is_3d(ds, v):
            _shape = (len(LFIDS), len(ds[ds[v].dims[0]]))
            COORDS = [('land_id', LFIDS), ('lf_id', ds.coords['lf_id'])]
        else:
            _shape = (len(LFIDS),)
            COORDS = [('land_id', LFIDS)]
            
        _blank = np.ones( _shape )
        _da = xr.DataArray(_blank[:], name=v, coords=COORDS)
        
        for lat, lon in zip(latL, lonL):
            land_id = D_ids[(lat, lon)]
            vals = ds[v].sel(lat=lat, lon=lon).to_masked_array().filled(NODATA)
            _da.loc[land_id] = vals
        
        _da.attrs.update( ds[v].attrs )
        _da.attrs.update(defaultD)

        dsout[_da.name] = _da

    return (ds_ids, dsout)


def mask_dataset(ds, valid):
    """Mask all values that are not valid/ 1 (2d or 3d)."""
    for v in ds.data_vars:
        dims = ds[v].dims
        if len(dims) > len(valid.shape):
            z = len(ds[v].values)
            valid = np.array(z*[valid])
        ds[v].values = np.ma.masked_where(valid == 0, ds[v].values).filled(NODATA)

    return ds

def create_gridlist(ds):
    """Create LPJ-Guess 4.0 gridlist file."""    
    outL = []
    for j in reversed(range(len(ds['land_id']))):
        for i in range(len(ds['land_id'][0])):
            x = ds['land_id'][j, i].values #to_masked_array()
            if x != NODATA: #p.ma.masked:
                lat = float(ds['lat'][j].values)
                lon = float(ds['lon'][i].values)
                land_id = int(ds['land_id'].sel(lat=lat, lon=lon).values)
                outS = "%3.2f %3.2f %d" % (lat, lon, land_id)
                outL.append(outS)
    
    return '\n'.join(outL) + '\n'
    

def main():
    """Main Script."""    
    
    # TODO: cleanup
    # call the old scripts (that have been reworked)
    call_customtpi()
    call_computeLandformStats5()
    
    # source files
    soilref = os.path.join('soil', 'GLOBAL_WISESOIL_DOM_05deg.nc')
    elevref = os.path.join('elevation_CL.nc')

    # build netcdfs
    print "building 2d netcdf files"
    sitenc = build_site_netcdf(soilref, elevref)
    landformnc = build_landform_netcdf(refnc=sitenc)
    
    # clip to joined mask
    elev_mask = np.where(sitenc['ELEVATION'].values == NODATA, 0, 1)
    landform_mask = np.where(landformnc['lfcnt'].values == NODATA, 0, 1)
    valid_mask = elev_mask * landform_mask
    
    sitenc = mask_dataset(sitenc, valid_mask)
    landformnc = mask_dataset(landformnc, valid_mask)

    # write 2d/ 3d netcdf files
    sitenc.to_netcdf('sites_2d.nc', format='NETCDF4_CLASSIC')
    landformnc.to_netcdf('landforms_2d.nc', format='NETCDF4_CLASSIC')

    # convert to compressed netcdf format
    print "building compressed format netcdf files"
    ids_2d, comp_sitenc = build_compressed(sitenc)
    _, comp_landformnc = build_compressed(landformnc)
    
    # write netcdf files
    ids_2d.to_netcdf("land_ids_2d.nc", format='NETCDF4_CLASSIC')
    comp_landformnc.to_netcdf("landform_data.nc", format='NETCDF4_CLASSIC')
    comp_sitenc.to_netcdf("site_data.nc", format='NETCDF4_CLASSIC')

    # gridlist file
    print "creating gridlist file"
    gridlist = create_gridlist(ids_2d)
    open("gridlist_CL.txt", 'w').write(gridlist)

    print "done"


if __name__ == '__main__':
    """ Global entry point """    
    main()


