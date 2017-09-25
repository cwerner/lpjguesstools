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
import logging
import math
import numpy as np
import os
import pandas as pd
import string
import xarray as xr

from _geoprocessing import compute_spatial_dataset, classify_aspect, \
                           classify_landform, get_center_coord, \
                           split_srtm1_dataset, get_global_attr, \
                           analyze_filename

log = logging.getLogger(__name__)

# import constants
from . import NODATA
from . import defaultAttrsDA

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


def define_landform_classes(step, limit, TYPE='SIMPLIFIED'):
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
    #  Name         SIMPLIFIED     WEISS
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
    #  Name         SIMPLIFIED     WEISS
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
    elif TYPE == 'SIMPLIFIED':
        # TYPE: SIMPLIFIED (1:hilltop, 3:midslope, 4:flat, 6:valley)
        lf_set = [10,31,33,40,60]
        lf_full_set = []
        for e in ele_cnt:
            lf_full_set += [x+(100*e) for x in lf_set]
    else:
        log.error('Currently only classifiation schemes WEISS, SIMPLIFIED supported.')


    return (lf_full_set, ele_breaks)


def tile_already_processed(fname, TILESTORE_PATH):
    """Check if the tile exists."""
    existing_tiles = glob.glob(os.path.join(TILESTORE_PATH, '*.nc'))
    #existing_tiles = [os.path.basename(x) for x in glob.glob(glob_string)]
    
    for existing_tile in existing_tiles:
        source_attr = get_global_attr(xr.open_dataset(existing_tile), 'source')
        if source_attr != None:
            # TODO: add second check (version?)
            _, source_name = analyze_filename(fname)
            if source_name == source_attr:
                return True
    return False    


def match_watermask_shpfile(glob_string):
    """Check if the generated shp glob_string exists."""
    if len(glob.glob(glob_string)) == 0:
        shp = None
    elif len(glob.glob(glob_string)) == 1:
        shp = glob.glob(glob_string)[0]
    else:
        log.error("Too many shape files.")
    return shp


def get_tile_summary(ds, cutoff=0):
    """Compute the fractional cover of the landforms in this tile."""
    
    unique, counts = np.unique(ds['landform_class'].to_masked_array().astype('i'), return_counts=True)
    remove_ix = np.where(unique == NODATA)
    unique = np.delete(unique, remove_ix)
    counts = np.delete(counts, remove_ix)
    total_valid = float(np.sum(counts))
    
    df = pd.DataFrame({'lf_id': unique, 'cells': counts})
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
    
    # cleanup (move lon, lat to front, drop coord col)
    df_.drop('coord', axis=1, inplace=True)
    latloncnt_cols = ['lon', 'lat', 'lf_cnt']
    new_col_order = latloncnt_cols + \
                    [x for x in df_.columns.tolist() if x not in latloncnt_cols]
    return df_[new_col_order]


def compute_landforms(glob_string, shp_mask_dir, tilestore_path='tiles', cutoff=1.0):
    """Compute landform units based on elevation, slope, aspect and tpi classes."""

    dem_files = sorted(glob.glob(glob_string))

    # define the final landform classes (now with elevation brackets)
    # SIMPLIFIED: reduced number of lf classes, more elevation differences
    lf_classes, lf_ele_levels = define_landform_classes(200, 6000, TYPE='SIMPLIFIED')

    for dem_file in dem_files:
        fname = os.path.basename(dem_file)
        fdir  = os.path.dirname(dem_file)
        
        # SRTM1 default nameing convention
        str_lat = fname[:3]
        str_lon = fname[3:7]
        
        # if tiles don't exist process them
        if not tile_already_processed(fname, tilestore_path):        
            log.info('processing: %s (%s)' % (dem_file, datetime.datetime.now()))

            shp_glob_string = os.path.join(shp_mask_dir, str_lon + str_lat + '*.shp')
            matched_shp_file = match_watermask_shpfile(shp_glob_string.lower())
            
            ds_srtm1 = compute_spatial_dataset(dem_file, fname_shp=matched_shp_file)
            tiles = split_srtm1_dataset(ds_srtm1)

            for tile in tiles:
                # reclass
                if tile != None:
                    classify_aspect(tile)
                    classify_landform(tile, elevation_levels=lf_ele_levels, TYPE='SIMPLIFIED')            
                    
                    # store file in tilestore
                    lon, lat = get_center_coord(tile)
                    lonlat_string = convert_float_coord_to_string((lon,lat))
                    tile.to_netcdf(os.path.join(tilestore_path, \
                                   "srtm1_processed_%s.nc" % lonlat_string)) 
                               
  
    # section 2
    log.info("START OF SECTION2")
    tiles = sorted(glob.glob(os.path.join(tilestore_path, '*.nc')))
    
    if not tile_files_compatible(tiles):
        log.error('Tile files in %s are not compatible.' % tilestore_path)

    tiles_stats = []
    for tile in tiles:
        ds = xr.open_dataset(tile, decode_cf=False)
        lf_stats = get_tile_summary(ds, cutoff=cutoff)
        number_of_ids = len(lf_stats)
        lon, lat = get_center_coord(ds)
        
        coords = pd.Series([(round(lon,2),round(lat,2), number_of_ids) for x in range(len(lf_stats))])
        lf_stats['coord'] = coords        
        lf_stats.set_index(['coord', 'lf_id'], inplace=True)
        tiles_stats.append( lf_stats )

    df = pd.concat(tiles_stats)

    frac_lf = create_stats_table(df, 'frac_scaled')
    elev_lf = create_stats_table(df, 'elevation')
    slope_lf = create_stats_table(df, 'slope')
    
    # return the dataframes and the list of all possible landform units
    return (frac_lf, elev_lf, slope_lf, lf_classes)


varD = {'TOTC': ('SOC', 'Soil Organic Carbon', 'percent', 0.1),
        'SDTO': ('SAND', 'Sand', 'percent', 1.0),
        'STPC': ('SILT', 'Silt', 'percent', 1.0),
        'CLPC': ('CLAY', 'Clay', 'percent', 1.0)}

soil_vars = sorted(varD.keys())


def is_3d(ds, v):
    """Check if xr.DataArray has 3 dimensions."""
    dims = ds[v].dims
    if len(dims) == 3:
        return True
    return False


def assign_to_dataarray(data, df, lf_full_set, refdata=False):
    """Place value into correct location of data array."""
    for _, r in df.iterrows():
        if refdata:
            data.loc[r.lat, r.lon] = r.lf_cnt
        else:
            data[:] = NODATA
            for lf in r.index[3:]:
                if r[lf] > NODATA:
                    lf_id = int(lf[2:])
                    lf_pos = lf_full_set.index(lf_id)
                    data.loc[lf_id, r.lat, r.lon] = r[lf]

    return data


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
        vattr.update(defaultAttrsDA)
        da.attrs.update(vattr)
        da[:] = np.ma.masked_where(emask, da.to_masked_array().filled(NODATA))
        dsout[da.name] = da

    # ele var
    da = xr.full_like(da.copy(deep=True), NODATA)
    da.name = 'ELEVATION'
    vattr = {'name': 'elevation', 'long_name': 'Elevation', 'units': 'meters'}
    vattr.update(defaultAttrsDA)
    da.attrs.update(vattr)
    
    da[:] = ds_ele['data'].to_masked_array().filled(NODATA)
    dsout[da.name] = da

    return dsout


def build_landform_netcdf(lf_full_set, frac_lf, elev_lf, slope_lf, refnc=None):
    """Build landform netcdf based on refnc dims and datatables."""
    
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
    
    # check that landform coordinates are in refnc
    lat_min, lat_max = frac_lf.lat.min(), frac_lf.lat.max()
    lon_min, lon_max = frac_lf.lon.min(), frac_lf.lon.max()
    
    lats = refnc['lat'].values.tolist()
    lons = refnc['lon'].values.tolist()
    
    if (((lat_min < min(lats)) or (lat_max > max(lats))) or 
       (((lon_min < min(lons)) or (lon_max > max(lons))))):
       log.error('DEM tiles not within specified extent.')
       exit()   
    
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
        dsout[dv].attrs.update(defaultAttrsDA)

    return dsout


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
    da_ids.attrs.update(defaultAttrsDA)
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
        _da.attrs.update(defaultAttrsDA)

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
    
    # defaults (will move to cli)
    SRTMSTORE_STRING = "srtm1/*.zip"
    WATERMASKSTORE_PATH = "srtm1_shp_mask"
    GRIDLIST_TXT = 'gridlist_CL.txt'
    # [6, 40, 12, 55] Deutschland
    REGION = #[-76, -56, -66, -16]    # lon1, lat1, lon2, lat2
    TILESTORE_PATH = 'processed'
    CUTOFF = 1.0    # % area required to keep landform

    # default soil and elevation data (contained in package)
    import pkg_resources
    SOIL_NC      = pkg_resources.resource_filename(__name__, 'data/GLOBAL_WISESOIL_DOM_05deg.nc')
    ELEVATION_NC = pkg_resources.resource_filename(__name__, 'data/GLOBAL_ELEVATION_05deg.nc')
    
    # current assumptions (make them cli options later):
    # - soil (use global ISRIC-WISE dataset in data dir)
    # - elevation (lat,lon vars: "latitude", "longitude"; name: "data")
    
    # section 1:
    # compute the 0.5x0.5 deg tiles
    log.info("computing landforms")
    
    # TODO: find a better way to access lf_full_set (instead of passing it around)
    df_frac, df_elev, df_slope, lf_full_set = compute_landforms(SRTMSTORE_STRING,
                                                                WATERMASKSTORE_PATH,
                                                                tilestore_path=TILESTORE_PATH,
                                                                cutoff=CUTOFF)
    
    # section 2:
    # build netcdfs
    log.info("building 2d netcdf files")
    sitenc = build_site_netcdf(SOIL_NC, ELEVATION_NC, extent=REGION)
    landformnc = build_landform_netcdf(lf_full_set, df_frac, df_elev, df_slope, refnc=sitenc)
    
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
    log.info("building compressed format netcdf files")
    ids_2d, comp_sitenc = build_compressed(sitenc)
    _, comp_landformnc = build_compressed(landformnc)
    
    # write netcdf files
    ids_2d.to_netcdf("land_ids_2d.nc", format='NETCDF4_CLASSIC')
    comp_landformnc.to_netcdf("landform_data.nc", format='NETCDF4_CLASSIC')
    comp_sitenc.to_netcdf("site_data.nc", format='NETCDF4_CLASSIC')

    # gridlist file
    log.info("creating gridlist file")
    gridlist = create_gridlist(ids_2d)
    open(GRIDLIST_TXT, 'w').write(gridlist)

    log.info("done")


if __name__ == '__main__':
    """ Global entry point """    
    main()


