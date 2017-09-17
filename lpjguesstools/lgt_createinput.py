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
import glob 
import math
import numpy as np
import os
from osgeo import osr, gdal, gdalconst
import pandas as pd
import rasterio as rio
import string
import xarray as xr

from _tpi import calculate_tpi
from _geoprocessing import process_dem, read_processed_geotiff


# consts and lookups
NODATA = -9999
defaultD = {'missing_value': NODATA, '_FillValue': NODATA}

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




def rename_lat_lon(lat, lon, TYPE):
    latS = 's'
    lonS = 'w'
    
    if lat >= 0: latS = 'n'
    if lon >= 0: lonS = 'e'
    
    outS = '%s%2.2f%s%2.2f_0.5arc_%s.nc' % (latS, math.fabs(lat), lonS, math.fabs(lon), TYPE)
    return outS


def define_landform_classes(step, limit):
    """Define the landform classes"""
    
    # Parameters:
    # - step: elevation interval for landform groups (def: 400m )
    # - limit: elevation limit [inclusive, in m]
    
    ele_breaks = [-1000] + range(step, limit, step) + [10000]
    ele_cnt = range(1, len(ele_breaks))

    # basic set of landforms
    # code: [slopeid<1..6>][aspectid<0,1..4>]
    lf_set = [10,21,22,23,24,31,32,33,34,40,51,52,53,54,60]

    lf_full_set = []
    for e in ele_cnt:
        lf_full_set += [x+(100*e) for x in lf_set]

    return lf_full_set


def compute_landforms(glob_string, shp_mask_dir):
    """Compute landform units based on elevation, slope, aspect and tpi classes."""


    CUTOFF = 1.0        # min percent landmass for LF to be considered
    TYPE = 'tpi300'     # tpi2000

    #for tile in sorted( glob.glob("landforms/*.tif") ):


    lf_full_set = define_landform_classes(400, 5000)
    
    # new landforms
    header = '\t'.join(['lat', 'lon', 'lf_cnt'] + ['LF%d' % x for x in lf_full_set]) + '\n'

    f_el = open('elevations_full_%s.txt' % TYPE, 'w')
    f_lf = open('landforms_full_%s.txt' % TYPE, 'w')
    f_sl = open('slopes_full_%s.txt' % TYPE, 'w')

    f_el.write(header)
    f_lf.write(header)
    f_sl.write(header)
    
    # Rules:
    # - slope and aspect calc should be conducted on water-clipped grid
    # - tpi calculations should be carried out on original data (since they)
    #   involve a potentially larger kernel (otherwise we have NODATA increase)

    #tiles = glob.glob('srtm1_filled/*.tif')
    tiles = sorted(glob.glob(glob_string))
    
    # create list of shp_mask files
    files = []
    
    for tile in tiles:
        fname = os.path.basename(tile)
        fdir  = os.path.dirname(tile)
        str_lat = fname[:3]
        str_lon = fname[4:8]
        
        # create shp file name
        shp_file = glob.glob(shp_mask_dir + '/' + str_lon + str_lat + '*.shp')
        if len(shp_file) == 0:
            tile = (tile, None)
        elif len(shp_file) == 1:
            tile = (tile, shp_file[0])
        else:
            log.error("Too many shape files.")

        print 'processing: ', tile, '(', datetime.datetime.now(), ')'

        DEBUG=False

        if DEBUG:
            rawdata = read_processed_geotiff(tile[0])
        else:
            dumpname = tile[0].replace('srtm1_filled', 'srtm1_dump_gtiff')
            rawdata = process_dem(tile[0], shp_mask=tile[1], dump=dumpname)
            
        DEM       = rawdata.GetRasterBand(1).ReadAsArray()
        SLOPE     = rawdata.GetRasterBand(2).ReadAsArray()
        ASPECT    = rawdata.GetRasterBand(3).ReadAsArray()
        WATERMASK = rawdata.GetRasterBand(4).ReadAsArray()
        
        # make sure we apply the existing nodata values
        DNODATA = rawdata.GetRasterBand(1).GetNoDataValue()        
        SNODATA = rawdata.GetRasterBand(2).GetNoDataValue()
        ANODATA = rawdata.GetRasterBand(3).GetNoDataValue()
        
        NODATA = rawdata.GetRasterBand(1).GetNoDataValue()
        if NODATA is not None:
            DEM    = np.ma.masked_equal(DEM, NODATA)
            SLOPE  = np.ma.masked_equal(SLOPE, NODATA)
            ASPECT = np.ma.masked_equal(ASPECT, NODATA)
            DEMM   = np.ma.masked_where(WATERMASK == 1, DEM) 

        # TODO: 
        # - check if we need to mask before this step
        # - convert to an xarray dataset(s) ?
        # - rename LF (and other stuff)

        dumpname = tile[0].replace('srtm1_filled', 'srtm1_dump_tpi')
        if DEBUG:
            LF = np.fromfile(dumpname)
        else:
            LF = calculate_tpi(DEM, SLOPE, 300, dump=dumpname)

        # we now work with DEM, DEMM, SLOPE, ASPECT, TPI300_CLASSES ()

        # get lower left coord of 1 deg tile, all in the west/south hemisphere
        b = [-int(x[1:]) for x in os.path.basename(tile[0]).split('_')[0:2]]
        blat, blon = b

        # landforms
        # tpi300_classed nodata/ from srtm1_masked
        #LNODATA = src.GetRasterBand(1).GetNoDataValue()
    
        # classify aspect
        ASPECT2 = np.ma.masked_where(np.asarray(ASPECT, 'i') < 0, np.asarray(ASPECT, 'i'))
        ASPECT2[(( ASPECT >= 315 ) | ((ASPECT>=0) & (ASPECT < 45)))] = 1    # north
        ASPECT2[( ASPECT >= 45  ) & ( ASPECT < 135)] = 2                    # east
        ASPECT2[( ASPECT >= 135 ) & ( ASPECT < 225)] = 3                    # south
        ASPECT2[( ASPECT >= 225 ) & ( ASPECT < 315)] = 4                    # west

        # the one and only base mask !!!
        DEMMASK = np.ma.masked_where( (DEM == DNODATA), np.ones_like(DEM))

        # secondary masks (are merged with the water mask of srtm1_masked)
        ASPMASK = np.where(ASPECT2 < 0, 1, 0) * np.where(DEMMASK.filled(NODATA) == NODATA, 1, 0)
        SLOMASK = np.where(SLOPE < 0, 1, 0) # * np.where(DEMMASK.filled(NODATA) == NODATA, 1, 0)
        LFMASK  = np.where(LF < 0, 1, 0) * np.where(DEMMASK.filled(NODATA) == NODATA, 1, 0)

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
            foutname = 'netcdfs_%s/' % TYPE + rename_lat_lon(lat, lon, TYPE)
            
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

                da = xr.DataArray(da[:].filled(NODATA),
                                  name=name,
                                  coords=[('lat', LATS[:]), ('lon', LONS[:])])
                da.attrs['units'] = Dunits[name]
                da.attrs['long_name'] = Dlname[name]
                da.attrs['missing_value'] = NODATA
                da.attrs['_FillValue'] = NODATA
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
    
    print "computing landforms"
    compute_landforms("srtm1_filled/*.tif", "srtm1_shp_mask")
    
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


