"""FILE lgt_createinput.main.py

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

from collections import OrderedDict
import datetime
import glob 
import logging
import math
import numpy as np
import os
import pandas as pd
import string
import time
import xarray as xr

from ._geoprocessing import analyze_filename_dem, \
                            classify_aspect, \
                            classify_landform, \
                            calculate_asp_slope, \
                            compute_spatial_dataset

from ._srtm1 import split_srtm1_dataset
import _xr_geo
import _xr_tile

log = logging.getLogger(__name__)

# import constants
from . import NODATA
from . import ENCODING

# quick helpers
# TODO: move to a dedicated file later

def time_dec(func):
    """A decorator to measure execution time of function"""
    def wrapper(*arg, **kwargs):
        t = time.time()
        res = func(*arg, **kwargs)
        log.debug('DURATION: <%s> : ' % func.func_name + str(time.time()-t))
        return res
    return wrapper


varSoil = {'TOTC': ('soc', 'Soil Organic Carbon', 'soc', 'percent', 0.1),
           'SDTO': ('sand', 'Sand', 'sand', 'percent', 1.0),
           'STPC': ('silt', 'Silt', 'silt', 'percent', 1.0),
           'CLPC': ('clay', 'Clay', 'clay', 'percent', 1.0)}

varLF = {'lfcnt': ('lfcnt', 'Number of landforms', 'lfcnt', '-', 1.0),
         'slope': ('slope', 'Slope', 'slope', 'deg', 1.0),
         'aspect': ('aspect', 'Aspect', 'aspect', 'deg', 1.0),         
         'asp_slope': ('asp_slope', 'Aspect-corrected Slope', 'asp_slope', 'deg', 1.0),         
         'fraction': ('fraction', 'Landform Fraction', 'fraction', '1/1', 1.0),
         'elevation': ('elevation', 'Elevation', 'elevation', 'm', 1.0)}

soil_vars = sorted(varSoil.keys())
lf_vars = sorted(varLF.keys())


def convert_float_coord_to_string(coord, p=2):
    """Convert a (lon,lat) coord to string."""
    lon, lat = round(coord[0], p), round(coord[1], p)
    LA, LO = 'n', 'e'
    if lat < 0: LA = 's'
    if lon < 0: LO = 'w'
    lat_s = "%.2f" % round(abs(lat),2)
    lon_s = "%.2f" % round(abs(lon),2)
    coord_s = '%s%s%s%s' % (LA, lat_s.zfill(p+3), LO, lon_s.zfill(p+4))
    return coord_s


def has_significant_land(ds, min_frac=0.01):
    """Test if land fraction in tile is significant."""
    # min_frac in %, default: 0.001 %
    if (ds['mask'].values.sum() / float(len(ds.lat.values) * len(ds.lon.values))) * 100 > min_frac:
        return True
    return False


def define_landform_classes(step, limit, TYPE='SIMPLE'):
    """Define the landform classes."""
    
    # Parameters:
    # - step: elevation interval for landform groups (def: 400m )
    # - limit: elevation limit [inclusive, in m]
    
    ele_breaks = [-1000] + range(step, limit, step) + [10000]
    ele_cnt = range(1, len(ele_breaks))

    # code system [code position 2 & 3, 1= elevations_tep]
    # code: [slopeid<1..6>][aspectid<0,1..4>]
    #
    #  slope:
    #
    #  Name           SIMPLE     WEISS
    #
    #  hilltop          1           1
    #  upper slope                  2*
    #  mid slope        3*          3*
    #  flats            4           4
    #  lower slope                  5*
    #  valley           6           6
    #
    #  
    #  aspect:
    #
    #  Name           SIMPLE     WEISS
    #
    #  north             1          1
    #  east              2          2
    #  south             3          3
    #  west              4          4
    
        
    if TYPE == 'WEISS':
        lf_set = [10,21,22,23,24,31,32,33,34,40,51,52,53,54,60]
        lf_full_set = []
        for e in ele_cnt:
            lf_full_set += [x+(100*e) for x in lf_set]
    elif TYPE == 'SIMPLE':
        # TYPE: SIMPLE (1:hilltop, 3:midslope, 4:flat, 6:valley)
        lf_set = [10,31,32,33,34,40,60]
        lf_full_set = []
        for e in ele_cnt:
            lf_full_set += [x+(100*e) for x in lf_set]
    else:
        log.error('Currently only classifiation schemes WEISS, SIMPLE supported.')


    return (lf_full_set, ele_breaks)


def tiles_already_processed(TILESTORE_PATH):
    """Check if the tile exists."""
    existing_tiles = glob.glob(os.path.join(TILESTORE_PATH, '*.nc'))
    #existing_tiles = [os.path.basename(x) for x in glob.glob(glob_string)]
    
    processed_tiles = []
    for existing_tile in existing_tiles:
        with xr.open_dataset(existing_tile) as ds:
            source = ds.tile.get('source')
            if source is not None:
                processed_tiles.append(source)
            else:
                log.warn('Source attr not set in file %s.' % existing_tile)
    return processed_tiles


def match_watermask_shpfile(glob_string):
    """Check if the generated shp glob_string exists."""
    found=False
    if len(glob.glob(glob_string)) == 0:
        shp = None
    elif len(glob.glob(glob_string)) == 1:
        shp = glob.glob(glob_string)[0]
        found = True
    else:
        log.error("Too many shape files.")
        exit()
                
    # second try: look for zip file
    if found is False:
        shp = glob_string.replace(".shp", ".zip")
        if len(glob.glob(shp)) == 0:
            shp = None
        elif len(glob.glob(shp)) == 1:
            shp = glob.glob(shp)[0]
        else:
            log.error("Too many shape files.")
            exit()
    return shp


def get_tile_summary(ds, cutoff=0):
    """Compute the fractional cover of the landforms in this tile."""
    
    unique, counts = np.unique(ds['landform_class'].to_masked_array(), return_counts=True)
    counts = np.ma.masked_array(counts, mask=unique.mask)
    unique = np.ma.compressed(unique)
    counts = np.ma.compressed(counts)
    total_valid = float(np.sum(counts))
    
    df = pd.DataFrame({'lf_id': unique.astype('int'), 'cells': counts})
    df['frac'] = (df['cells'] / df['cells'].sum())*100

    df = df[df['frac'] >= cutoff]
    df['frac_scaled'] = (df['cells'] / df['cells'].sum())*100
    
    # also get lf-avg of elevation and slope
    df['elevation'] = -1
    df['slope'] = -1
    df['asp_slope'] = -1
    df['aspect'] = -1

    a_lf = ds['landform_class'].to_masked_array()
    
    # average aspect angles
    def avg_aspect(a):
        x = 0
        y = 0
        for v in a.ravel():
            x += math.sin(math.radians(v))
            y += math.cos(math.radians(v))
        avg = math.degrees(math.atan2(x, y))
        if avg < 0:
            avg += 360
        return avg 
    
    
    # calculate the avg. elevation and slope in landforms
    for i, r in df.iterrows():
        ix = a_lf == int(r['lf_id'])
        lf_slope = ds['slope'].values[ix].mean()
        lf_asp_slope = ds['asp_slope'].values[ix].mean()        
        lf_elevation = ds['elevation'].values[ix].mean()
        lf_aspect = avg_aspect(ds['aspect'].values[ix])
        df.loc[i, 'slope'] = lf_slope
        df.loc[i, 'asp_slope'] = lf_asp_slope
        df.loc[i, 'elevation'] = lf_elevation
        df.loc[i, 'aspect']    = lf_aspect
    return df


def tile_files_compatible(files):
    """Get global attribute from all tile netcdf files and check
    they were created with an identical elevation step.
    """

    fingerprints = []
    for file in files:
        with xr.open_dataset(file) as ds:
            fingerprint = (ds.tile.get('elevation_step'), ds.tile.get('classification'))
        fingerprints.append(fingerprint)
    
    # check if elements are equal    
    if all(x==fingerprints[0] for x in fingerprints):
        # check if there are Nones' in any fingerprint
        if not all(fingerprints):
            return False
        return True
    return False


def create_stats_table(df, var):
    """Create a landform info table for all coords and given var."""

    df_ = df[var].unstack(level=-1, fill_value=NODATA)
    # rename columns and split coord tuple col to lon and lat col
    df_.columns = ['lf' + str(col) for col in df_.columns]
    if 'lf0' in df_.columns:
        del df_['lf0']
    df_ = df_.reset_index()
    df_[['lon', 'lat', 'lf_cnt']] = df_['coord'].apply(pd.Series)
    df_['lf_cnt'] = df_['lf_cnt'].astype(int)
    # cleanup (move lon, lat to front, drop coord col)
    df_.drop('coord', axis=1, inplace=True)
    latloncnt_cols = ['lon', 'lat', 'lf_cnt']
    new_col_order = latloncnt_cols + \
                    [x for x in df_.columns.tolist() if x not in latloncnt_cols]
    return df_[new_col_order]


@time_dec
def convert_dem_files(cfg, lf_ele_levels):
    """Compute landform units based on elevation, slope, aspect and tpi classes."""

    if cfg.SRTMSTORE_PATH is not None:
        
        # if glob_string is a directory, add wildcard for globbing
        glob_string = cfg.SRTMSTORE_PATH
        if os.path.isdir(cfg.SRTMSTORE_PATH):
            glob_string = os.path.join(cfg.SRTMSTORE_PATH, '*')
        dem_files = sorted(glob.glob(glob_string))

        existing_tiles = tiles_already_processed(cfg.TILESTORE_PATH)

        for dem_file in dem_files:
            fname = os.path.basename(dem_file)
            fdir  = os.path.dirname(dem_file)
            
            # SRTM1 default nameing convention
            str_lat = fname[:3]
            str_lon = fname[3:7]
            
            # if tiles don't exist process them
            process_tiles = True
            if cfg.OVERWRITE:
                process_tiles = True
            else:
                _, source_name = analyze_filename_dem(fname)
                if source_name in existing_tiles:
                    process_tiles = False
            
            if process_tiles:                        
                log.info('processing: %s (%s)' % (dem_file, datetime.datetime.now()))

                shp_glob_string = os.path.join(cfg.WATERMASKSTORE_PATH, str_lon + str_lat + '*.shp')
                matched_shp_file = match_watermask_shpfile(shp_glob_string.lower())
                
                ds_srtm1 = compute_spatial_dataset(dem_file, fname_shp=matched_shp_file)
                tiles = split_srtm1_dataset(ds_srtm1)

                for i, tile in enumerate(tiles):

                    # reclass
                    if tile != None and has_significant_land(tile):
                        log.debug("Valid tile %d in file %s." % (i+1, dem_file))

                        classify_aspect(tile)
                        classify_landform(tile, elevation_levels=lf_ele_levels, TYPE=cfg.CLASSIFICATION)            
                        calculate_asp_slope(tile)
                        
                        # store file in tilestore
                        # get tile center coordinate and name
                        lon, lat = tile.geo.center()
                        lonlat_string = convert_float_coord_to_string((lon,lat))
                        tile_name = "srtm1_processed_%s.nc" % lonlat_string
                        tile.to_netcdf(os.path.join(cfg.TILESTORE_PATH, tile_name), \
                                                    format='NETCDF4_CLASSIC')
                    else:
                        log.debug("Empty tile %d in file %s ignored." % (i+1, dem_file))

    
@time_dec
def compute_statistics(cfg):
    """Extract landform statistics from tiles in tilestore."""
    available_tiles = glob.glob(os.path.join(cfg.TILESTORE_PATH, '*.nc'))
    log.debug('Number of tiles found: %d' % len(available_tiles)) 
    if len(available_tiles) == 0:
        log.error('No processed tiles available in directory "%s"' % cfg.TILESTORE_PATH)
        exit()
    tiles = sorted(available_tiles)
    
    if not tile_files_compatible(tiles):
        log.error('Tile files in %s are not compatible.' % cfg.TILESTORE_PATH)
        exit()

    tiles_stats = []
    for tile in tiles:
        log.debug('Computing statistics for tile %s' % tile)
        with xr.open_dataset(tile) as ds:
            lf_stats = get_tile_summary(ds, cutoff=cfg.CUTOFF)
            lf_stats.reset_index(inplace=True)
            number_of_ids = len(lf_stats)
            lon, lat = ds.geo.center()
            coord_tuple = (round(lon,2),round(lat,2), int(number_of_ids)) 
            lf_stats['coord'] = pd.Series([coord_tuple for _ in range(len(lf_stats))])
            lf_stats.set_index(['coord', 'lf_id'], inplace=True)
            tiles_stats.append( lf_stats )

    df = pd.concat(tiles_stats)

    frac_lf = create_stats_table(df, 'frac_scaled')
    elev_lf = create_stats_table(df, 'elevation')
    slope_lf = create_stats_table(df, 'slope')
    asp_slope_lf = create_stats_table(df, 'asp_slope')
    aspect_lf = create_stats_table(df, 'aspect')    
    return (frac_lf, elev_lf, slope_lf, asp_slope_lf, aspect_lf)


def is_3d(ds, v):
    """Check if xr.DataArray has 3 dimensions."""
    dims = ds[v].dims
    if len(dims) == 3:
        return True
    return False


def assign_to_dataarray(data, df, lf_full_set, refdata=False):
    """Place value into correct location of data array."""
    
    if refdata==True:
        data[:] = NODATA
    else:
        data[:] = np.nan
    for _, r in df.iterrows():
        if refdata:
            data.loc[r.lat, r.lon] = r.lf_cnt
        else:
            for lf in r.index[3:]:
                if r[lf] > NODATA:
                    lf_id = int(lf[2:])
                    lf_pos = lf_full_set.index(lf_id)
                    data.loc[dict(lf_id=lf_id, lat=r.lat, lon=r.lon)] = r[lf]
    return data


def spatialclip_df(df, extent):
    """Clip dataframe wit lat lon columns by extent."""
    if any(e is None for e in extent):
        log.warn("SpatialClip: extent passed is None.")
    lon1, lat1, lon2, lat2 = extent

    if ('lon' not in df.columns) or ('lat' not in df.columns):
        log.warn("SpatialClip: lat/ lon cloumn missing in df.")
    return df[((df.lon >= lon1) & (df.lon <= lon2)) & 
              ((df.lat >= lat1) & (df.lat <= lat2))]


def build_site_netcdf(soilref, elevref, extent=None):
    """Build the site netcdf file."""
    
    # extent: (x1, y1, x2, y2)
    ds_soil_orig = xr.open_dataset(soilref)
    ds_ele_orig = xr.open_dataset(elevref)

    if extent is not None:
        lat_min, lat_max = extent[1], extent[3]
        lon_min, lon_max = extent[0], extent[2]

        # slice simulation domain
        ds_soil = ds_soil_orig.sel(lon=((ds_soil_orig.lon >= lon_min) & (ds_soil_orig.lon <= lon_max)),
                                   lat=((ds_soil_orig.lat >= lat_min) & (ds_soil_orig.lat <= lat_max)),
                                   lev=1.0).squeeze(drop=True)
        ds_ele = ds_ele_orig.sel(longitude=((ds_ele_orig.longitude >= lon_min) & (ds_ele_orig.longitude <= lon_max)),
                                 latitude=((ds_ele_orig.latitude >= lat_min) & (ds_ele_orig.latitude <= lat_max))).squeeze(drop=True)
    else:
        ds_soil = ds_soil_orig.sel(lev=1.0).squeeze(drop=True)
        ds_ele = ds_ele_orig.squeeze(drop=True)
    del ds_soil['lev']

    # identify locations that need filling and use left neighbor
    smask = np.where(ds_soil['TOTC'].to_masked_array().mask, 1, 0)
    emask = np.where(ds_ele['data'].to_masked_array().mask, 1, 0)

    # no soil data but elevation: gap-fill wioth neighbors
    missing = np.where((smask == 1) & (emask == 0), 1, 0)
    ix, jx = np.where(missing == 1)
    
    if len(ix) > 0:
        log.debug('Cells with elevation but no soil data [BEFORE GF: %d].' % len(ix))
        
        for i, j in zip(ix, jx):
            for v in soil_vars:
                if np.isfinite(ds_soil[v][i, j-1]):
                    ds_soil[v][i, j] = ds_soil[v][i, j-1].copy(deep=True)
                elif np.isfinite(ds_soil[v][i, j+1]):
                    ds_soil[v][i, j] = ds_soil[v][i, j+1].copy(deep=True)
                else:
                    print 'neighbours have nodata !!!'
                x = ds_soil[v][i, j].to_masked_array()

        smask2 = np.where(ds_soil['TOTC'].to_masked_array().mask, 1, 0)
        missing = np.where((smask2 == 1) & (emask == 0), 1, 0)
        ix, jx = np.where(missing == 1)
        log.debug('Cells with elevation but no soil data [AFTER GF:  %d].' % len(ix))

    dsout = xr.Dataset()
    # soil vars
    for v in soil_vars:
        conv = varSoil[v][-1]
        da = ds_soil[v].copy(deep=True) * conv
        da.name = varSoil[v][0]

        vattr = {'name': varSoil[v][0],
                 'long_name': varSoil[v][1],
                 'standard_name': varSoil[v][2],
                 'units': varSoil[v][3],
                 'coordinates': "lat lon"}
                 
        da.tile.update_attrs(vattr)
        da.tile.update_encoding(ENCODING)
         
        da[:] = np.ma.masked_where(emask, da.to_masked_array())
        dsout[da.name] = da

    # ele var
    da = xr.full_like(da.copy(deep=True), np.nan)
    da.name = 'elevation'
    vattr = {'name': 'elevation', 'long_name': 'Elevation', 
             'units': 'meters',  'standard_name': 'elevation'}
    da.tile.update_attrs(vattr)
    da.tile.update_encoding(ENCODING)

    da[:] = ds_ele['data'].to_masked_array()
    dsout[da.name] = da
    return dsout


@time_dec
def build_landform_netcdf(lf_full_set, frac_lf, elev_lf, slope_lf, asp_slope_lf, aspect_lf, cfg, elevation_levels, refnc=None):
    """Build landform netcdf based on refnc dims and datatables."""
    
    dsout = xr.Dataset()

    COORDS = [('lf_id', lf_full_set), ('lat', refnc.lat), ('lon', refnc.lon)]
    SHAPE = tuple([len(x) for _, x in COORDS])
    
    # initiate data arrays
    _blank = np.empty(SHAPE)
    da_lfcnt = xr.DataArray(_blank.copy()[0,:,:].astype(int), name='lfcnt', 
                            coords=COORDS[1:])
    da_frac = xr.DataArray(_blank.copy(), name='fraction', coords=COORDS)
    da_slope = xr.DataArray(_blank.copy(), name='slope', coords=COORDS)
    da_asp_slope = xr.DataArray(_blank.copy(), name='asp_slope', coords=COORDS)
    da_elev = xr.DataArray(_blank.copy(), name='elevation', coords=COORDS)
    da_aspect = xr.DataArray(_blank.copy(), name='aspect', coords=COORDS)
    
    # check that landform coordinates are in refnc
    df_extent = [frac_lf.lon.min(), frac_lf.lat.min(), frac_lf.lon.max(), frac_lf.lat.max()]
    log.debug('df_extent: %s' % str(df_extent))
    log.debug('contains: %s' % str(refnc.geo.contains(df_extent)))
    
    if refnc.geo.contains(df_extent) == False:
        
        frac_lf = spatialclip_df(frac_lf, refnc.geo.extent)
        slope_lf = spatialclip_df(slope_lf, refnc.geo.extent)
        asp_slope_lf = spatialclip_df(asp_slope_lf, refnc.geo.extent)
        elev_lf = spatialclip_df(elev_lf, refnc.geo.extent)
        aspect_lf = spatialclip_df(aspect_lf, refnc.geo.extent)

    # dump files
    frac_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_frac.csv'), index=False)
    slope_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_slope.csv'), index=False)
    asp_slope_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_asp_slope.csv'), index=False)
    elev_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_elev.csv'), index=False)
    aspect_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_aspect.csv'), index=False)    
 
    # assign dataframe data to arrays
    da_lfcnt = assign_to_dataarray(da_lfcnt, frac_lf, lf_full_set, refdata=True)
    da_frac = assign_to_dataarray(da_frac, frac_lf, lf_full_set)
    da_slope = assign_to_dataarray(da_slope, slope_lf, lf_full_set)
    da_asp_slope = assign_to_dataarray(da_asp_slope, asp_slope_lf, lf_full_set)
    da_elev = assign_to_dataarray(da_elev, elev_lf, lf_full_set)
    da_aspect = assign_to_dataarray(da_aspect, aspect_lf, lf_full_set)

    # store arrays in dataset
    dsout[da_lfcnt.name] = da_lfcnt
    dsout[da_frac.name] = da_frac
    dsout[da_slope.name] = da_slope
    dsout[da_asp_slope.name] = da_asp_slope
    dsout[da_elev.name] = da_elev
    dsout[da_aspect.name] = da_aspect

    for v in dsout.data_vars:
        vattr = {}
        if v in lf_vars:
            vattr = {'name': varLF[v][0],
                     'long_name': varLF[v][1],
                     'standard_name': varLF[v][2],
                     'units': varLF[v][3],
                     'coordinates': "lat lon"}
        dsout[v].tile.update_attrs(vattr)
        dsout[v].tile.update_encoding(ENCODING)
    
    dsout['lat'].tile.update_attrs(dict(standard_name='latitude',
                                        long_name='latitude',
                                        units='degrees_north'))

    dsout['lon'].tile.update_attrs(dict(standard_name='longitude',
                                        long_name='longitude',
                                        units='degrees_east'))

    dsout['lf_id'].tile.update_attrs(dict(standard_name='lf_id',
                                          long_name='lf_id',
                                          units='-'))
    for dv in dsout.data_vars:
        dsout[dv].tile.update_encoding(ENCODING)

    # register the specific landform properties (elevation steps, classfication)
    dsout.tile.set('elevation_step', elevation_levels[1])
    dsout.tile.set('classification', cfg.CLASSIFICATION.lower())

    return dsout


def build_compressed(ds):
    """Build LPJ-Guess 4.0 compatible compressed netcdf file."""
    # identify landforms netcdf
    if 'lfcnt' in ds.data_vars:
        v = 'lfcnt'
    elif 'elevation' in ds.data_vars:
        v = 'elevation'
    else:
        log.error("Not a valid xr.Dataset (landforms or site only).")

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

    lats.tile.update_attrs(dict(standard_name='latitude',
                                long_name='latitude',
                                units='degrees_north'))

    lons.tile.update_attrs(dict(standard_name='longitude',
                                long_name='longitude',
                                units='degrees_east'))

    # create land_id reference array
    # TODO: clip land_id array to Chile country extent?
    da_ids.tile.update_encoding(ENCODING)
    ds_ids = da_ids.to_dataset(name='land_id')

    # create xr.Dataset
    dsout = xr.Dataset()
    dsout[lats.name] = lats
    dsout[lons.name] = lons

    # walk through variables, get lat/ lon cells' data 
    for v in ds.data_vars:
        if is_3d(ds, v):
            _shape = (len(LFIDS), len(ds[ds[v].dims[0]]))
            COORDS = [('land_id', LFIDS), ('lf_id', ds['lf_id'])]
        else:
            _shape = (len(LFIDS),)
            COORDS = [('land_id', LFIDS)]
        _blank = np.ones( _shape )
        _da = xr.DataArray(_blank[:], name=v, coords=COORDS)
        
        for lat, lon in zip(latL, lonL):
            land_id = D_ids[(lat, lon)]
            vals = ds[v].sel(lat=lat, lon=lon).to_masked_array()
            _da.loc[land_id] = vals
        
        _da.tile.update_attrs(ds[v].attrs)
        _da.tile.update_encoding(ENCODING)

        dsout[_da.name] = _da

        if is_3d(ds, v):
            dsout['lf_id'].tile.update_attrs(dict(standard_name='lf_id',
                                                  long_name='lf_id',
                                                  units='-'))

    # copy lgt attributes from ssrc to dst
    dsout.tile.copy_attrs(ds)

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
    

def main(cfg):
    """Main Script."""    
    
    # default soil and elevation data (contained in package)
    import pkg_resources
    SOIL_NC      = pkg_resources.resource_filename(__name__, '../data/GLOBAL_WISESOIL_DOM_05deg.nc')
    ELEVATION_NC = pkg_resources.resource_filename(__name__, '../data/GLOBAL_ELEVATION_05deg.nc')
    
    log.info("Converting DEM files and computing landform stats")

    # define the final landform classes (now with elevation brackets)
    lf_classes, lf_ele_levels = define_landform_classes(200, 6000, TYPE=cfg.CLASSIFICATION)

    # process dem files to tiles (if not already processed)
    convert_dem_files(cfg, lf_ele_levels)

    #sitenc = build_site_netcdf(SOIL_NC, ELEVATION_NC, extent=cfg.REGION)

    # compute stats from tiles

    df_frac, df_elev, df_slope, df_asp_slope, df_aspect = compute_statistics(cfg)
    #print 'reading files'
    #df_frac = pd.read_csv('lfdata.cutoff_1.0p/df_frac.csv')
    #df_asp_slope = pd.read_csv('lfdata.cutoff_1.0p/df_asp_slope.csv')
    #df_slope = pd.read_csv('lfdata.cutoff_1.0p/df_slope.csv')
    #df_aspect = pd.read_csv('lfdata.cutoff_1.0p/df_aspect.csv')
    #df_elev = pd.read_csv('lfdata.cutoff_1.0p/df_elev.csv')

    # build netcdfs
    log.info("Building 2D netCDF files")
    sitenc = build_site_netcdf(SOIL_NC, ELEVATION_NC, extent=cfg.REGION)
    landformnc = build_landform_netcdf(lf_classes, df_frac, df_elev, df_slope, df_asp_slope, df_aspect, cfg, 
                                       lf_ele_levels, refnc=sitenc)
    
    # clip to joined mask
    #elev_mask = np.where(sitenc['elevation'].values == NODATA, 0, 1)
    #landform_mask = np.where(landformnc['lfcnt'].values == NODATA, 0, 1)
    #valid_mask = elev_mask * landform_mask


    elev_mask = ~np.ma.getmaskarray(sitenc['elevation'].to_masked_array())
    sand_mask = ~np.ma.getmaskarray(sitenc['sand'].to_masked_array())
    land_mask = ~np.ma.getmaskarray(landformnc['lfcnt'].to_masked_array())
    valid_mask = elev_mask * sand_mask * land_mask

 
    sitenc = mask_dataset(sitenc, valid_mask)
    landformnc = mask_dataset(landformnc, valid_mask)

    # write 2d/ 3d netcdf files
    sitenc.to_netcdf(os.path.join(cfg.OUTDIR, 'sites_2d.nc'), 
                     format='NETCDF4_CLASSIC')
    landformnc.to_netcdf(os.path.join(cfg.OUTDIR, 'landforms_2d.nc'),
                         format='NETCDF4_CLASSIC')

    # convert to compressed netcdf format
    log.info("Building compressed format netCDF files")
    ids_2d, comp_sitenc = build_compressed(sitenc)
    ids_2db, comp_landformnc = build_compressed(landformnc)
    
    # write netcdf files
    ids_2d.to_netcdf(os.path.join(cfg.OUTDIR, "land_ids_2d.nc"), 
                     format='NETCDF4_CLASSIC')

    ids_2db.to_netcdf(os.path.join(cfg.OUTDIR, "land_ids_2db.nc"),
                     format='NETCDF4_CLASSIC')

    comp_landformnc.to_netcdf(os.path.join(cfg.OUTDIR, "landform_data.nc"), 
                              format='NETCDF4_CLASSIC')
    comp_sitenc.to_netcdf(os.path.join(cfg.OUTDIR, "site_data.nc"), 
                          format='NETCDF4_CLASSIC')

    # gridlist file
    log.info("Creating gridlist file")
    gridlist = create_gridlist(ids_2d)
    open(os.path.join(cfg.OUTDIR, cfg.GRIDLIST_TXT), 'w').write(gridlist)

    log.info("Done")


