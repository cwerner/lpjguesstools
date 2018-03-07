"""FILE lgt_convert.main.py

This script converts standard and SubPixel 
LPJ-GUESS .out files to netCDF files

Christian Werner, SENCKENBERG Biodiversity and Climate Research Centre (BiK-F)
email: christian.werner@senkenberg.de
2017/11/07
"""

from collections import OrderedDict
import logging
import numpy as np
import pandas as pd
import xarray as xr
import os
import sys

from .extra import set_config, get_config, parse_config #, RefDataBuilder

__version__ = "0.0.2"

log = logging.getLogger(__name__)

# default attributes for netCDF variable of dataarrays
NODATA = -9999
defaultAttrsDA = {'_FillValue': NODATA, 'missing_value': NODATA}

# functions

def is_monthly(df):
    """Check if df is a monthly dataframe."""
    col_names = df.columns.values
    return (('Jan' in col_names) and ('Dec' in col_names))

def is_pft(df):
    """Check if df is a per-pft dataframe."""
    col_names = df.columns.values
    return 'Total' in col_names

def has_stand(df):
    """Check if df is a subpixel dataframe."""
    return 'Stand' in df.columns.values

def get_data_column_index(df):
    """Determine the start and end position of data columns."""

    col_names = df.columns.values.tolist()

    # determine start position of data
    if 'Patch' in col_names:
        cid_start = col_names.index('Patch') + 1
    elif 'Stand' in col_names:
        cid_start = col_names.index('Stand') + 1
    elif 'Year' in col_names:
        cid_start = col_names.index('Year') + 1
    else:
        cid_start = col_names.index('Lat') + 1

    # determine end position of data
    cid_end   = len(col_names)
    if is_monthly(df):
        cid_end   = col_names.index('Dec')
    if is_pft(df):
        cid_end = col_names.index('Total') - 1
    return (cid_start, cid_end)

def _read_global_info(cfg):
    """Read general project info from cfg file."""
    info = parse_config(cfg, section='info')
    project = parse_config(cfg, section='project')
    all_info = OrderedDict()
    if info:
        for k in info.keys():
            all_info[k] = info[k]
    else:
        log.warn("No <info> data found in config")
    if project:
        for k in project.keys():
            all_info[k] = project[k]
    else:
        log.warn("No <project> data found in config")
    return all_info

def _read_data_info(cfg):
    """Read data section from cfg file."""
    return parse_config(cfg, section='data')

def _read_refdata_info(cfg):
    rd = parse_config(cfg, section='refdata')
    xmin, ymin, xmax, ymax = tuple(rd['bbox'])
    step = rd['res']
    hstep = step * 0.5
    lats = np.arange(ymin+hstep, ymax, hstep)
    lons = np.arange(xmin+hstep, xmax, hstep)
    return (lats, lons)

def update_attrs(obj, *args, **kwargs):
    """Update netcdf attributes."""
    obj.attrs.update(*args, **kwargs)
    return obj


class IndexMapper():
    def __init__(self, lats, lons):
        if ((type(lats) == list) and (type(lons) == list)):
            self._lons = lons
            self._lats = lats
        elif ((type(lats) == np.ndarray) and (type(lons) == np.ndarray)):
            if ((len(lons) == 1) and (len(lats) == 1)):
                self._lons = lons.tolist()
                self._lats = lats.tolist()
            else:
                log.critical("Only 1D numpy arrays supported")
                log.critical(lats)
                log.critical(lons)
                exit(1)
        elif ((type(lats) == xr.DataArray) and (type(lons) == xr.DataArray)):
            self._lons = lons.values.tolist()
            self._lats = lats.values.tolist()
        else:
            log.critical("Only 1D numpy arrays or lists supported")
            log.critical(lats)
            log.critical(lons)
            exit(1)

    def __call__(self, lat, lon):
        if ((lon >= self._lons[0]) and (lon <= self._lons[-1])):
            ix = self._lons.index(lon)
        else:
            log.critical("Lon %s out of range [%f...%f]" % (lon,
                                                     self._lons[0],
                                                     self._lons[-1]))
            exit(1)

        if ((lat >= self._lats[0]) and (lat <= self._lats[-1])):
            jx = self._lats.index(lat)
        else:
            log.critical("Lon %s out of range [%f...%f]" % (lon,
                                                     self._lats[0],
                                                     self._lats[-1]))
            exit(1)
        return (jx, ix)


def get_annual_data(var, landforms, df_frac, cfg,
                    sel_var=None, lf_ids=None):
    """ Parse variable and return DataArrays (data, total_data) """

    if lf_ids is None:
        KEEP_LF_DIM = False
    else:
        KEEP_LF_DIM = True

    # derive dimensions and landform info from landforms file
    lats = landforms.coords['lat']
    lons = landforms.coords['lon']

    # map lat, lon to indices
    mapper = IndexMapper(lats, lons)

    var_str = var
    if 'sp_' in var:
        var_str = var.replace('sp_', '')

    # read either regular or compressed data files
    try:
        df = pd.read_csv(os.path.join(cfg.INDIR, "%s.out" % var),
            delim_whitespace=True)
    except:
        df = pd.read_csv(os.path.join(cfg.INDIR, "%s.out.gz" % var),
            delim_whitespace=True)

    if 'Stand' in df.columns.values:
        if cfg.DEFAULT:
            log.info('Selecting only lfid 0 data rows (for now).')
            df = df[df.Stand == 0]
            KEEP_LF_DIM = False
        else:
            df = df[df.Stand != 0]
    
    # special transformation for FireRT variable
    if sel_var == 'FireRT':
        log.debug("  Inverting fire return time")
        df['FireRT'] = 1.0 / df['FireRT']
    
    # add a dummy Total column for identification purposes
    if var_str == 'height':
        df['Total'] = 1.
    
    log.debug("  Total number of data rows in file (raw data): %d" % len(df))

    # limit df to years (either specified years or use last nyears)
    if cfg.LAST_NYEARS is not None:
        max_yr = df.Year.max()
        cfg.overwrite(dict(YEARS=range(max_yr - (cfg.LAST_NYEARS-1), max_yr+1)))

    if len(cfg.YEARS) > 0:
        log.debug("  Limiting years")
        df = df[ (df.Year >= cfg.YEARS[0]) & (df.Year <= cfg.YEARS[-1])]

    if len(df) == 0:
        log.critical("Requested years not in data.")
        exit(1)
    
    log.debug("  Total number of data rows in file (year sel): %d" % len(df))
    # determine z dimension
    #nyears = df.Year.max() - df.Year.min()
    outyears = range(df.Year.min(), df.Year.max()+1)

    # safety check: if years > 2000 and not site mode
    #               possible user error (and memory kill likely) 
    if len(outyears) > 2000 and cfg.SMODE == False:
        log.warn("Large number of years (%d) in data and region mode" % (len(outyears)))
        log.warn("selected. Are you sure you do not want -S (site mode)? This")
        log.warn("will likely crash due to insufficient memory.")
        
    if var == 'sp_fpc':
        log.info("Capping to 100% - total")
        
        cid_start, cid_end = get_data_column_index(df)
        col_names = df.columns.values.tolist()
        a = col_names[cid_start]
        b = col_names[cid_end]    
        mask = (df['Total'] > 1)    # all rows that need scaling
        df.loc[mask,a:'Total'] = df.loc[mask,a:'Total'].mul(1/df['Total'], axis=0)

    # calc mean over patches
    if 'Patch' in df.columns.values:
        groupcols = ['Lon','Lat','Year','Stand']
        df = df.groupby(groupcols).mean().reset_index()
        del df['Patch']

    # return to proper column sortin Lon,Lat,Year,(Stand)
    if has_stand(df):
        # fix index before merge
        df.set_index(['Lon','Lat','Stand'], inplace=True)# set index for join
        
        # if default (lfid=0), we do not join but give a default fraction
        if cfg.DEFAULT:
            log.info(df_frac.head())
            df['fraction'] = 1
        else:
            df = df.merge(df_frac, left_index=True, right_index=True)

    df_yrs = []

    df.reset_index(inplace=True)
    if has_stand(df):
        df.set_index(['Lon','Lat','Year','Stand'], inplace=True)
        df.reset_index(inplace=True)
        
    if cfg.AVG:
        # average years
        log.debug("  Averaging years over selected timespan.")
        if 'Stand' in df.columns.values:
            groupcols = ['Lon','Lat','Stand']
            df = df.groupby(groupcols).mean().reset_index()
        else:
            groupcols = ['Lon','Lat']
            df = df.groupby(groupcols).mean().reset_index()
        del df['Year']
        df['Year'] = 0
        # place this into yearly list with tuple index value: None 
        df_yrs = [(0, df)]
        outyears = [0]


    log.debug("  Total number of data rows in file (annual avg): %d" % len(df))

    data_cols = [c for c in df.columns.values if c not in ['Lat', 'Lon', 'Year', 'Stand']]


    # do the stand/ lf_id avg

    def weighed_average(grp):
        return grp._get_numeric_data().multiply(grp['fraction'], axis=0).sum()/grp['fraction'].sum()

    if 'Stand' in df.columns.values and KEEP_LF_DIM == False:
        # if average is requested, do (weighted) average over the year column, too    

        df_yrs = []
        # split by year for performance considerations
        for cnt, yr in enumerate(outyears):
            yr_mask = df.Year == yr
            df_yr = df[yr_mask]        
            df_yrs.append((yr, df_yr))

        # loop over years
        new_dfs = []
        groupcols = ['Lon','Lat','Year'] #,'Stand']
        for yr, df_yr in df_yrs:
            log.debug("  Processing yr %s" % str(yr))
            #new_df = df.groupby(groupcols).apply(wavg).reset_index() #.compute()
            
            new_df = df_yr.groupby(groupcols).apply(weighed_average)
            del new_df['Lat']
            del new_df['Lon']
            del new_df['Year']
            del new_df['Stand']
            
            new_dfs.append(new_df)
        df = pd.concat(new_dfs)
        df.columns = data_cols
        df.reset_index(inplace=True)

    
    if sel_var == 'FireRT':
        log.debug("  Fire Return Interval not inverted: do so with 1/FireRT yourself")
    #    df['FireRT'] = 1.0 / df['FireRT']

    log.debug("  Total number of data rows in file (final):    %d" % len(df))

    cid_start, cid_end = get_data_column_index(df)
    if is_pft(df):
        PFTS = df.columns.values.tolist()[cid_start:cid_end+1]
    else:
        PFTS = []

    class DataStore(object):
        def __init__(self, years, **kwargs):
            """ Initialice with list of requested dims """
            
            if ('lat' not in kwargs.keys()) or ('lon' not in kwargs.keys()):
                log.error("We need 'lat' and 'lon' to constuct array.")
                exit()
                 
            # dimensions order
            full_set_dims = ['time', 'time_m', 'month', 'lf_id', 'pft', 'lat', 'lon']
            
            if cfg.AVG:
                self.years = [0]
            
            self.dim_names = []
            self.dim_sizes = []
            self.data = None
            self.data_total = None
            
            if cfg.AVG:
                self.years = [0]
            else:
                self.years = years
            # special year list for monthly time axis
            self.yearsm = [item for item in self.years for i in range(12)]
            self.lfids = []
            self.pfts = []
            
            if 'lf_id' in kwargs.keys():
                self.lfids = kwargs['lf_id'].values.tolist()
            
            if 'pft' in kwargs.keys():
                self.pfts = kwargs['pft'].values.tolist()
                
            # build coords
            da_coords = []
            da_coords2 = []
            for d in full_set_dims:
                if d in kwargs.keys():
                    value = kwargs[d]
                    da_dim = (d, value)
                    da_coords.append(da_dim)
                    if d != 'pft':
                        da_coords2.append(da_dim)
                    self.dim_names.append(d)
                    self.dim_sizes.append(len(value))    

            self.data = xr.DataArray(np.ones( self.dim_sizes ) * np.nan,
                                     coords=da_coords)

            self.data.attrs['units'] = '-'
            self.data.attrs['_FillValue'] = NODATA

            if 'pft' in kwargs.keys():
                p = self.dim_names.index('pft')
                dim_names = self.dim_names[:]
                dim_sizes = self.dim_sizes[:]
                dim_names.pop(p)
                dim_sizes.pop(p)
                self.data_total = xr.DataArray(np.ones( dim_sizes ) * np.nan,
                                               coords=da_coords2)

                self.data_total.attrs['units'] = '-'
                self.data_total.attrs['_FillValue'] = NODATA

        def add(self, year, jx, ix, values, lf_id=None, sel_var=None):
            if cfg.AVG:
                year_p = 0
            else:
                if 'time_m' in self.dim_names:
                    year_p = self.yearsm.index(int(year)) 
                else:
                    year_p = self.years.index(int(year))
                

            if 'pft' in self.dim_names:
                if 'lf_id' in self.dim_names:
                    lfid_p = self.lfids.index(lf_id)
                    self.data[year_p, lfid_p, :, jx, ix] = values
                    self.data_total[year_p, lfid_p, jx, ix] = sum(values)                
                else:
                    self.data[year_p, :, jx, ix] = values
                    self.data_total[year_p, jx, ix] = sum(values)
            elif 'month' in self.dim_names or 'time_m' in self.dim_names:
                if 'month' in self.dim_names:
                    add_p = self.dim_names.index('month')
                    if 'lf_id' in self.dim_names:
                        lfid_p = self.lfids.index(lf_id) 
                        self.data[year_p, :, lfid_p, jx, ix] = values
                    else:
                        self.data[year_p, :, jx, ix] = values
                        
                elif 'time_m' in self.dim_names:
                    add_p = self.dim_names.index('time_m')
                    if 'lf_id' in self.dim_names:
                        lfid_p = self.lfids.index(lf_id) 
                        self.data[year_p:year_p+12, lfid_p, jx, ix] = values
                    else:
                        self.data[year_p:year_p+12, jx, ix] = values
            else:
                # individual variable
                if 'lf_id' in self.dim_names:
                    lfid_p = self.lfids.index(lf_id)
                    self.data[year_p, lfid_p, jx, ix] = values[sel_var]
                else:
                    self.data[year_p, jx, ix] = values[sel_var]

        def fetch(self):
            return self.data, self.data_total
        
        def mask(self):
            dims = list(self.data.dims)
            dims.pop(dims.index('lat'))
            dims.pop(dims.index('lon'))
            x = self.data.sum(dim=dims, skipna=True).values
            return np.where(x==0,0,1)


    # setup output array
    mcoords = None  # extra month axis
    lfcoords = None # extra lf_ids axis

    # build the coordinate dictionary
    if cfg.SMODE:
        # single site mode, simply take lat, lon from first line of data
        site_lon = df.loc[0, 'Lon']
        site_lat = df.loc[0, 'Lat']
        coords = dict(lat=[site_lat], lon=[site_lon])
    else:
        # default mode, take lat, lon from landforms_2d file
        coords = dict(lat=lats, lon=lons)

        # modify coords for single site mode


    if cfg.AVG:
        len_outyears = 1
    else:
        len_outyears = len(outyears)

    # second time axis (monthly substeps)
    zcoords2 = None

    if KEEP_LF_DIM:
        # lf_id axis
        lfcoords = xr.DataArray( lf_ids, name='lf_id') 
        lfcoords.attrs['units'] = '-'
        coords['lf_id'] = lfcoords

    # add time axis
    zcoords = xr.DataArray( outyears, name='time')
    zcoords.attrs['units'] = 'yearly'

    if is_pft(df):
        log.debug("  We have a pft variable.")
        coords['pft'] = xr.DataArray( PFTS, name='pft')

    # create 1 or 2 time dims (yr, yr+month)
    if cfg.USE_MONTH_DIM:
        mcoords = xr.DataArray( range(12), name='month')
        mcoords.attrs['units'] = 'month'

    if is_monthly(df):
        log.debug("  We have a monthly variable.")
        if not cfg.USE_MONTH_DIM:
            zcoords2 = xr.DataArray( np.linspace(outyears[0], outyears[-1],
                                                 num=len(outyears)*12,
                                                 endpoint=False), name='time')
            
            # another time dim
            coords['time_m'] = zcoords2
        else:
            coords['month'] = mcoords
            coords['time'] = zcoords
    else:
        coords['time'] = zcoords

    # build the datastore
    dstore = DataStore( outyears, **coords )

    for _, row in df.iterrows():
        # map lat lon position to index
        if cfg.SMODE:
            # site mode: only one lat, lon position
            jx, ix = 0, 0
        else:
            jx, ix = mapper(row.Lat, row.Lon)
        rdata = row[cid_start: cid_end+1]
        if cfg.AVG:
            _year = 0
        else:
            _year = row.Year
        if 'Stand' in row.index:
            dstore.add(_year, jx, ix, rdata, lf_id=row.Stand, sel_var=sel_var)
        else:            
            dstore.add(_year, jx, ix, rdata, sel_var=sel_var)    

    data, data2 = dstore.fetch()
    mask = dstore.mask()
    return (data, data2, mask)


def main(cfg):
    """Main routine."""

    cfgdata = get_config(cfgFile=cfg.CONFIG)

    if cfg.STORECONFIG:
        set_config(cfgdata)

    def use_cli_refdata():
        return cfg.REFINFO is not None

    def get_fractions(ds_lf, refvar):
        """Parse landform fractions from reference file"""
        
        if refvar not in ds_lf.data_vars:
            log.critical("Var <%s> not in %s" % ('fraction', reffile))
            exit(1)
        else:
            log.info("Get fractional landform cover. Reading var %s from %s..." % (refvar, reffile))
            da = ds_lf[refvar]
            # iterate over lf_ids
            dfs = []
            for lf_id in ds_lf.lf_id:
                # fix nan
                x = da.sel(lf_id=lf_id).to_pandas().fillna(-9999)
                df = da.sel(lf_id=lf_id).to_pandas().stack().reset_index()
                
                # replace -9999 with nan again
                df.replace(to_replace=-9999, value=np.nan, inplace=True)
                df.columns = ['Lat','Lon','fraction']
                if len(df) > 0:
                    df['Stand'] = lf_id.values
                    dfs.append(df)            

            df = pd.concat(dfs, ignore_index=True)
            df.set_index(['Lon', 'Lat', 'Stand'], inplace=True)

        return df

    if use_cli_refdata():
        rdata = cfg.REFINFO
        reffile, refvar = rdata

        if os.path.isfile(reffile):
            with xr.open_dataset(reffile).load() as ds_lf:
                if refvar is not None:
                    df_frac = get_fractions(ds_lf, refvar)
                else:
                    refvar = 'fraction'
                    df_frac = get_fractions(ds_lf, refvar)
                lf_ids = ds_lf.lf_id.values
                da_frac = ds_lf['fraction']
                da_lfcnt = ds_lf['lfcnt']
    else:
        log.error("We currently need a refdata file!")
        exit(-1)
                        
    # derive data from cli or config file
    global_info = _read_global_info(cfgdata)

    if not os.path.isdir(cfg.INDIR):
        log.critical("Specified input directory does not exist.\n")
        exit(1)
    
    outname = cfg.OUTNAME

    data_files = [x[0] for x in _read_data_info(cfgdata)]

    # output netcdf file
    ds = xr.Dataset()

    masks = []

    for file in data_files:

        # use plain filename or also derive selected vars and renamed vars
        named_vars={}
        if type(file) == tuple:
            file, named_vars = file
        log.debug("Processing file %s ..." % file)
        
        # create two extra sets for north / south differences
        # if we have caused by exposition

        def is_subpixel_file(x):
            return 'sp_' in x

        # TODO: split into seperate netCDF files
        ##if is_subpixel_file(file) and args.north_south:
        ##    log.debug("  Calculating south, north and full landform average.")
        ##    sset_list = ['north', 'south', '']
        ##else:
        log.debug("  Calculating full landform average.")
        sset_list = ['']

        for sset in sset_list:
            suffix = ''
            if sset != '':
                suffix = '_%s' % sset

            # use filename as variable, also selected column/ var if requested    
            # use last_nyears instead of actual year range
            if cfg.LAST_NYEARS is not None:
                log.debug("  Using last nyears %d" % cfg.LAST_NYEARS)
            
            if len(named_vars.keys()) > 0:
                sel_vars = named_vars.keys()
                for var in sel_vars:
                    da1, _, ma = get_annual_data(file, ds_lf, df_frac, cfg,
                                            sel_var=var,
                                            lf_ids=list(lf_ids))
                    # rename variable
                    ds[named_vars[var] + suffix] = da1
            else:
                da1, da2, ma = get_annual_data(file, ds_lf, df_frac, cfg,
                                        lf_ids=list(lf_ids))
            
                ds[file + suffix] = da1
                # if we have a second (total) dataarray
                if type(da2) == xr.DataArray:
                    ds['tot_%s%s' % (file, suffix)] = da2

            masks.append(ma)

    mask = np.sum(masks, axis=0)
    mask = np.where(mask==0,np.nan,1)
    da_mask = xr.DataArray(mask, coords=[da1.coords['lat'], da1.coords['lon']], dims=['lat','lon'])

    # all potential coordinate variables
    all_c_vars = ['time','time_m','month','pft','lat','lon']
    c_vars = [x for x in all_c_vars if x in ds.coords.keys()]

    # order variables (coords first, then by name)
    d_vars = [x for x in sorted([x for x in ds.data_vars]) if x not in c_vars]

    # fraction variable
    ds['fraction'] = da_frac * da_mask
    ds['lfcnt'] = da_lfcnt * da_mask
    d_vars = ['lfcnt', 'fraction'] + d_vars

    # add site_lat, site_lon to file
    if cfg.SMODE:
        site_lon = ds.coords['lon'].values[0]        
        site_lat = ds.coords['lat'].values[0]
        ds.attrs['site_lon'] = site_lon
        ds.attrs['site_lat'] = site_lat
        
    ds = ds[c_vars + d_vars].squeeze(drop=True)
    
    # add compression
    comp = dict(zlib=True, complevel=5)
    
    if cfg.DEFAULT == False:
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(cfg.OUTNAME[:-3] + '_lfid.nc', format='NETCDF4_CLASSIC', unlimited_dims='time', encoding=encoding)

    # create lf average version of file
    dvars = [x for x in ds.data_vars if x not in ['lfcnt', 'fraction']] 
    dsout = xr.Dataset()
    
    scaled_fraction = ds['fraction']/ds['fraction'].sum(dim='lf_id')
    
    for dv in dvars:
        if cfg.DEFAULT == False:
            dsout[dv] = (ds[dv] * scaled_fraction).sum(dim='lf_id').squeeze(drop=True).where(ds['lfcnt']>0)
        else:
            dsout[dv] = (ds[dv]).where(ds['lfcnt']>0)
    
    # add site_lat, site_lon to file
    if cfg.SMODE:
        dsout.attrs['site_lon'] = site_lon
        dsout.attrs['site_lat'] = site_lat
    
    encoding = {var: comp for var in dsout.data_vars}    
    dsout.to_netcdf(cfg.OUTNAME, format='NETCDF4_CLASSIC', unlimited_dims='time', encoding=encoding)

    log.debug("Done.")
