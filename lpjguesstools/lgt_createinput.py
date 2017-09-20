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
                           split_srtm1_dataset, get_global_attr

log = logging.getLogger(__name__)

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


def define_landform_classes(step, limit):
    """Define the landform classes."""
    
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

    return (lf_full_set, ele_breaks)


def tile_already_processed(fname, TILESTORE_PATH):
    """Check if the tile exists."""
    existing_tiles = glob.glob(os.path.join(TILESTORE_PATH, '*.nc'))
    #existing_tiles = [os.path.basename(x) for x in glob.glob(glob_string)]
    
    for existing_tile in existing_tiles:
        print existing_tile
        source_attr = get_global_attr(xr.open_dataset(existing_tile), 'source')
        if source_attr != None:
            # TODO: add second check (version?)
            if fname == source_attr:
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


def compute_landforms(glob_string, shp_mask_dir):
    """Compute landform units based on elevation, slope, aspect and tpi classes."""

    # config settings (that eventually should go to cli or conf file)
    TILESTORE_PATH = "processed"        # location for final 0.5x0.5 deg tiles
    CUTOFF = 1.0                        # min percent covered by LF to be considered
    
    dem_files = sorted(glob.glob(glob_string))

    # define the final landform classes (now with elevation brackets)
    lf_classes, lf_ele_levels = define_landform_classes(400, 6000)

    for dem_file in dem_files:
        fname = os.path.basename(dem_file)
        fdir  = os.path.dirname(dem_file)
        str_lat = fname[:3]
        str_lon = fname[4:8]
        
        # if tiles don't exist process them
        if not tile_already_processed(fname, TILESTORE_PATH):        
            log.info('processing: %s (%s)' % (dem_file, datetime.datetime.now()))

            shp_glob_string = shp_mask_dir + '/' + str_lon + str_lat + '*.shp'
            matched_shp_file = match_watermask_shpfile(shp_glob_string)
            
            ds_srtm1 = compute_spatial_dataset(dem_file, fname_shp=matched_shp_file)
            tiles = split_srtm1_dataset(ds_srtm1)

            for tile in tiles:
                # reclass
                classify_aspect(tile)
                classify_landform(tile, elevation_levels=lf_ele_levels)            
                
                # store file in tilestore
                lon, lat = get_center_coord(tile)
                lonlat_string = convert_float_coord_to_string((lon,lat))
                tile.to_netcdf(os.path.join(TILESTORE_PATH, \
                               "srtm1_processed_%s.nc" % lonlat_string)) 
                               
                               
    # section 2
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
        df['dem'] = -1
        df['slope'] = -1

        a_lf = ds['landform_class'].to_masked_array()
        
        # calculate the avg. elevation and slope in landforms
        for i, r in df.iterrows():
            print r['lf_id'], type(r['lf_id'])
            
            ix = a_lf == int(r['lf_id'])
            lf_slope = ds['slope'].values[ix].mean()
            lf_elevation = ds['dem'].values[ix].mean()
            df.loc[i, 'slope'] = lf_slope
            df.loc[i, 'dem'] = lf_elevation
            
        #df = df.sort_values(by='cells', ascending=False)
        df.reset_index(inplace=True)
        return df
    

    def tile_files_compatible(files):
        """Get global attribute from all tile netcdf files and check
        they are the same.
        """
        x = [get_global_attr(xr.open_dataset(x), 'landform_elevation_step') for x in files]
        if all(x):
            if x[0] != None:
                return True
        else:
            return False

    tiles = sorted(glob.glob(os.path.join(TILESTORE_PATH, '*.nc')))
    
    if not tile_files_compatible(tiles):
        log.error('Tile files in %s are not compatible.' % TILESTORE_PATH)

    tiles_stats = []
    for tile in tiles:
        ds = xr.open_dataset(tile, decode_cf=False)
        lf_stats = get_tile_summary(ds, cutoff=CUTOFF)
        lon, lat = get_center_coord(ds)
        
        coords = pd.Series([(round(lon,2),round(lat,2)) for x in range(len(lf_stats))])
        lf_stats['coord'] = coords
        lf_stats.set_index(['coord', 'lf_id'], inplace=True)
        tiles_stats.append( lf_stats )

    df = pd.concat(tiles_stats)

    for col in ['frac_scaled', 'slope', 'dem']:
        dfx = df[col]
        dfx = dfx.unstack(level=-1, fill_value=-9999)
        
        # rename columns and split coord tuple col to lon and lat col
        dfx.columns = ['lf' + str(col) for col in dfx.columns]
        dfx = dfx.reset_index()
        dfx[['lon', 'lat']] = dfx['coord'].apply(pd.Series)
        
        # cleanup (move lon, lat to front, drop coord col)
        dfx.drop('coord', axis=1, inplace=True)
        latlon_cols = ['lon', 'lat']
        new_col_order = ['lon', 'lat'] + \
                        [x for x in dfx.columns.tolist() if x not in latlon_cols]
        dfx = dfx[new_col_order]
        print dfx
        print '----------------'


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
    
    SRTMSTORE_STRING = "srtm1_filled/*.tif"
    WATERMASKSTORE_PATH = "srtm1_shp_mask"
    
    # section 1:
    # compute the 0.5x0.5 deg tiles
    print "computing landforms"
    compute_landforms(SRTMSTORE_STRING, WATERMASKSTORE_PATH)
    
    # section 2:
    # build the actual LPJ-Guess 4.0 subpixel input files
    # ...
    exit()
    
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


