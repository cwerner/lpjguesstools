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

import click
from collections import Counter, OrderedDict
import copy
import datetime
import glob 
import logging
import numpy as np
import os
import pandas as pd
import string
import time
import xarray as xr

from _geoprocessing import compute_spatial_dataset, classify_aspect, \
                           classify_landform, get_center_coord, \
                           split_srtm1_dataset, get_global_attr, \
                           analyze_filename_dem, set_global_attr, \
                           update_attrs, update_encoding

log = logging.getLogger(__name__)

# import constants
from . import NODATA
from . import ENCODING
from . import EPILOG

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


class Bunch(object):
    """Simple data storage class."""
    def __init__(self, adict):
        self.__dict__.update(adict)
        

varD = {'TOTC': ('SOC', 'Soil Organic Carbon', 'percent', 0.1),
        'SDTO': ('SAND', 'Sand', 'percent', 1.0),
        'STPC': ('SILT', 'Silt', 'percent', 1.0),
        'CLPC': ('CLAY', 'Clay', 'percent', 1.0)}

soil_vars = sorted(varD.keys())


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
    #  east              2
    #  south             3          3
    #  west              4          
    
        
    if TYPE == 'WEISS':
        lf_set = [10,21,22,23,24,31,32,33,34,40,51,52,53,54,60]
        lf_full_set = []
        for e in ele_cnt:
            lf_full_set += [x+(100*e) for x in lf_set]
    elif TYPE == 'SIMPLE':
        # TYPE: SIMPLE (1:hilltop, 3:midslope, 4:flat, 6:valley)
        lf_set = [10,31,33,40,60]
        lf_full_set = []
        for e in ele_cnt:
            lf_full_set += [x+(100*e) for x in lf_set]
    else:
        log.error('Currently only classifiation schemes WEISS, SIMPLE supported.')


    return (lf_full_set, ele_breaks)


def tile_already_processed(fname, TILESTORE_PATH):
    """Check if the tile exists."""
    existing_tiles = glob.glob(os.path.join(TILESTORE_PATH, '*.nc'))
    #existing_tiles = [os.path.basename(x) for x in glob.glob(glob_string)]
    
    for existing_tile in existing_tiles:
        source_attr = get_global_attr(xr.open_dataset(existing_tile), 'source')
        if source_attr != None:
            # TODO: add second check (version?)
            _, source_name = analyze_filename_dem(fname)
            if source_name == source_attr:
                return True
    return False    


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

    a_lf = ds['landform_class'].to_masked_array()
    
    # calculate the avg. elevation and slope in landforms
    for i, r in df.iterrows():
        ix = a_lf == int(r['lf_id'])
        lf_slope = ds['slope'].values[ix].mean()
        lf_elevation = ds['elevation'].values[ix].mean()
        df.loc[i, 'slope'] = lf_slope
        df.loc[i, 'elevation'] = lf_elevation
        
    #df = df.sort_values(by='cells', ascending=False)
    df.reset_index(inplace=True)
    return df


def tile_files_compatible(files):
    """Get global attribute from all tile netcdf files and check
    they are the same.
    """
    x = [get_global_attr(xr.open_dataset(x), 'lgt.elevation_step') for x in files]
    if all(x):
        if x[0] != None:
            return True
    else:
        return False


def create_stats_table(df, var):
    """Create a landform info table for all coords and given var."""
    
    df_ = df[var].unstack(level=-1, fill_value=NODATA)
    # rename columns and split coord tuple col to lon and lat col
    df_.columns = ['lf' + str(col) for col in df_.columns]
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
        if os.path.isdir(cfg.SRTMSTORE_PATH):
            glob_string = os.path.join(cfg.SRTMSTORE_PATH, '*')
        dem_files = sorted(glob.glob(cfg.SRTMSTORE_PATH))

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
                if tile_already_processed(fname, cfg.TILESTORE_PATH):
                    process_tiles = False
            
            if process_tiles:                        
                log.info('processing: %s (%s)' % (dem_file, datetime.datetime.now()))

                shp_glob_string = os.path.join(cfg.WATERMASKSTORE_PATH, str_lon + str_lat + '*.shp')
                matched_shp_file = match_watermask_shpfile(shp_glob_string.lower())
                
                ds_srtm1 = compute_spatial_dataset(dem_file, fname_shp=matched_shp_file)
                tiles = split_srtm1_dataset(ds_srtm1)

                for tile in tiles:
                    # reclass
                    if tile != None:
                        classify_aspect(tile)
                        classify_landform(tile, elevation_levels=lf_ele_levels, TYPE=cfg.CLASSIFICATION)            
                        
                        # store file in tilestore
                        lon, lat = get_center_coord(tile)
                        lonlat_string = convert_float_coord_to_string((lon,lat))
                        tile.to_netcdf(os.path.join(cfg.TILESTORE_PATH, \
                                       "srtm1_processed_%s.nc" % lonlat_string), 
                                       format='NETCDF4_CLASSIC')

    
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
            number_of_ids = len(lf_stats)
            lon, lat = get_center_coord(ds)
            
            coords = pd.Series([(round(lon,2),round(lat,2), int(number_of_ids)) for x in range(len(lf_stats))])
            lf_stats['coord'] = coords        
            lf_stats.set_index(['coord', 'lf_id'], inplace=True)
            tiles_stats.append( lf_stats )

    df = pd.concat(tiles_stats)

    frac_lf = create_stats_table(df, 'frac_scaled')
    elev_lf = create_stats_table(df, 'elevation')
    slope_lf = create_stats_table(df, 'slope')
    return (frac_lf, elev_lf, slope_lf)


def is_3d(ds, v):
    """Check if xr.DataArray has 3 dimensions."""
    dims = ds[v].dims
    if len(dims) == 3:
        return True
    return False


def assign_to_dataarray(data, df, lf_full_set, refdata=False):
    """Place value into correct location of data array."""

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


def spatialclip_dataframe(df, extent):
    """Clip dataframe wit lat lon columns by extent."""
    lon1, lat1, lon2, lat2 = extent
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
    missing = np.where((smask == 1) & (emask == 0), 1, 0)

    ix, jx = np.where(missing == 1)
    for i, j in zip(ix, jx):
        for v in soil_vars:
            ds_soil[v][i, j] = ds_soil[v][i, j-1]

    dsout = xr.Dataset()
    # soil vars
    for v in soil_vars:
        conv = varD[v][-1]
        da = ds_soil[v].copy(deep=True) * conv
        da.name = varD[v][0]

        vattr = {'name': varD[v][0],
                 'long_name': varD[v][1],
                 'units': varD[v][2]}
                 
        da.pipe(update_attrs, vattr)
        da.pipe(update_encoding, ENCODING)    
        da[:] = np.ma.masked_where(emask, da.to_masked_array())
        dsout[da.name] = da

    # ele var
    da = xr.full_like(da.copy(deep=True), np.nan)
    da.name = 'ELEVATION'
    vattr = {'name': 'elevation', 'long_name': 'Elevation', 'units': 'meters'}
    da.pipe(update_attrs, vattr)
    da.pipe(update_encoding, vattr)

    da[:] = ds_ele['data'].to_masked_array()
    dsout[da.name] = da
    return dsout


@time_dec
def build_landform_netcdf(lf_full_set, frac_lf, elev_lf, slope_lf, cfg, elevation_levels, refnc=None):
    """Build landform netcdf based on refnc dims and datatables."""
    
    dsout = xr.Dataset()

    COORDS = [('lf_id', lf_full_set), ('lat', refnc.lat), ('lon', refnc.lon)]
    SHAPE = tuple([len(x) for _, x in COORDS])
    
    # initiate data arrays
    _blank = np.ones(SHAPE) * NODATA
    da_lfcnt = xr.DataArray(_blank.copy()[0,:,:].astype('i'), name='lfcnt', 
                            coords=COORDS[1:])
    da_frac = xr.DataArray(_blank.copy(), name='frac', coords=COORDS)
    da_slope = xr.DataArray(_blank.copy(), name='slope', coords=COORDS)
    da_elev = xr.DataArray(_blank.copy(), name='elevation', coords=COORDS)
    
    # check that landform coordinates are in refnc
    lat_min, lat_max = frac_lf.lat.min(), frac_lf.lat.max()
    lon_min, lon_max = frac_lf.lon.min(), frac_lf.lon.max()
    
    lats = refnc['lat'].values.tolist()
    lons = refnc['lon'].values.tolist()
    
    if ((lat_min < min(lats)) | (lat_max > max(lats)) |  
        (lon_min < min(lons)) | (lon_max > max(lons))):
        log.warn('DEM tiles not within specified extent. Clipping.')

    # potentially clip dataframes
    frac_lf = spatialclip_dataframe(frac_lf, [min(lons), min(lats), max(lons), max(lats)])
    slope_lf = spatialclip_dataframe(slope_lf, [min(lons), min(lats), max(lons), max(lats)])
    elev_lf = spatialclip_dataframe(elev_lf, [min(lons), min(lats), max(lons), max(lats)])

    # dump files
    frac_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_frac.csv'), index=False)
    slope_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_slope.csv'), index=False)
    elev_lf.to_csv(os.path.join(cfg.OUTDIR, 'df_elev.csv'), index=False)    
    # assign dataframe data to arrays
    da_lfcnt = assign_to_dataarray(da_lfcnt, frac_lf, lf_full_set, refdata=True)
    da_frac = assign_to_dataarray(da_frac, frac_lf, lf_full_set)
    da_slope = assign_to_dataarray(da_slope, slope_lf, lf_full_set)
    da_elev = assign_to_dataarray(da_elev, elev_lf, lf_full_set)

    # store arrays in dataset
    dsout[da_lfcnt.name] = da_lfcnt
    dsout[da_frac.name] = da_frac
    dsout[da_slope.name] = da_slope
    dsout[da_elev.name] = da_elev
    for dv in dsout.data_vars:
        dsout[dv].pipe(update_encoding, ENCODING)

    # register the specific landform properties (elevation steps, classfication)
    set_global_attr(dsout, 'lgt.elevation_step', elevation_levels[1])
    set_global_attr(dsout, 'lgt.classification', cfg.CLASSIFICATION.lower())

    return dsout


def copy_global_lgt_attrs(ds, dsout):
    """Copy global lgt attributes from source to target dataset."""
    source_attrs = dict([(k, v) for k, v in ds.attrs.items() if 'lgt.' in k])
    dsout.attrs.update(source_attrs)

def build_compressed(ds):
    """Build LPJ-Guess 4.0 compatible compressed netcdf file."""
    # identify landforms netcdf
    if 'lfcnt' in ds.data_vars:
        v = 'lfcnt'
    elif 'ELEVATION' in ds.data_vars:
        v = 'ELEVATION'
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

    # create land_id reference array
    # TODO: clip land_id array to Chile country extent?
    da_ids.pipe(update_encoding, ENCODING)
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
            vals = ds[v].sel(lat=lat, lon=lon).to_masked_array()
            _da.loc[land_id] = vals
        
        _da.pipe(update_attrs, ds[v].attrs)
        _da.pipe(update_encoding, ENCODING)

        dsout[_da.name] = _da

    copy_global_lgt_attrs(ds, dsout)

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
    SOIL_NC      = pkg_resources.resource_filename(__name__, 'data/GLOBAL_WISESOIL_DOM_05deg.nc')
    ELEVATION_NC = pkg_resources.resource_filename(__name__, 'data/GLOBAL_ELEVATION_05deg.nc')
    
    log.info("Converting DEM files and computing landform stats")

    # define the final landform classes (now with elevation brackets)
    lf_classes, lf_ele_levels = define_landform_classes(200, 6000, TYPE=cfg.CLASSIFICATION)

    # process dem files to tiles (if not already processed)
    convert_dem_files(cfg, lf_ele_levels)

    sitenc = build_site_netcdf(SOIL_NC, ELEVATION_NC, extent=cfg.REGION)

    # compute stats from tiles
    df_frac, df_elev, df_slope = compute_statistics(cfg)
    
    # build netcdfs
    log.info("Building 2D netCDF files")
    sitenc = build_site_netcdf(SOIL_NC, ELEVATION_NC, extent=cfg.REGION)
    landformnc = build_landform_netcdf(lf_classes, df_frac, df_elev, df_slope, cfg, 
                                       lf_ele_levels, refnc=sitenc)
    
    # clip to joined mask
    elev_mask = np.where(sitenc['ELEVATION'].values == NODATA, 0, 1)
    landform_mask = np.where(landformnc['lfcnt'].values == NODATA, 0, 1)
    valid_mask = elev_mask * landform_mask
    
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
    _, comp_landformnc = build_compressed(landformnc)
    
    # write netcdf files
    ids_2d.to_netcdf(os.path.join(cfg.OUTDIR, "land_ids_2d.nc"), 
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


# command line arguments
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
click.Context.get_usage = click.Context.get_help

@click.command(context_settings=CONTEXT_SETTINGS, epilog=EPILOG)

@click.option('--classfication', type=click.Choice(['SIMPLE', 'WEISS']), 
                    default='SIMPLE', show_default=True,
                    help='classification scheme')

@click.option('--cutoff', default=1.0, show_default=True,
                    help='required area fraction [%]')

@click.option('--dems', metavar='PATH',
                    help='source for DEM files')

@click.option('--extent', nargs=4, type=click.FLOAT, metavar='LON1 LAT1 LON2 LAT2',
                    help='extent of output netcdf files')

@click.option('--force-overwrite', is_flag=True, default=False, 
                    help='overwrite tiles even if they already exists')

@click.option('--gridlist', default='gridlist.txt', 
                    help='name of created gridlist file')

@click.option('--masks', metavar='PATH',
                    help='source for water masks (shp)')

@click.option('--verbose', is_flag=True, 
                    help='increase logging info')

@click.version_option()

@click.argument('storage', type=click.Path(exists=True))
@click.argument('outdir', type=click.Path(exists=True)) 

def cli(cutoff, dems, masks, gridlist, extent, classfication, storage, outdir, verbose, force_overwrite):
    """LPJ-GUESS 4.0 subpixel mode input creation tool
    
    This tools creates site and landform netCDF files and a gridlist file
    from SRTM1 (or potentially other) elevation data.
     
    """
    
    # example:
    #./lgt_createinput.py processed output --dems=srtm1/*.zip --masks=srtm1_shp_mask --extent -76 -56 -66 -16

    if verbose:
        logging.getLogger(__name__).setLevel(logging.DEBUG)
    else:
        logging.getLogger(__name__).setLevel(logging.INFO)
        
    if dems is not None:
        SRTMSTORE_PATH = dems
    
    if masks is not None:
        WATERMASKSTORE_PATH = masks

    REGION = None
    if len(extent) == 4:
        REGION = list(extent)
        
    
    # the setup dictionary to convert into a bunch obj
    config_data=dict(SRTMSTORE_PATH=dems,
                     WATERMASKSTORE_PATH=masks,
                     TILESTORE_PATH=storage,
                     REGION=REGION,
                     CLASSIFICATION=classfication,
                     CUTOFF=cutoff,
                     OUTPUT_PATH=outdir,
                     GRIDLIST_TXT=gridlist,
                     OUTDIR=outdir,
                     OVERWRITE=force_overwrite)
    
    # TODO: change logging level based on verbose flag
    cfg = Bunch(config_data)

    main(cfg)
    