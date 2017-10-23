# -*- coding: utf-8 -*-
#
# lgt_convert.py
# ==============
#
# Christian Werner
# christian.werner@senckenberg.de
"""lpjguesstools.lgt_convert: provides entry point main()."""

from collections import OrderedDict
import logging
import numpy as np
import pandas as pd
import xarray as xr
import os
import sys
import dask.dataframe as dd

from .cli import cli
from .extra import set_config, get_config, parse_config #, RefDataBuilder

__version__ = "0.0.1"

log = logging.getLogger(__name__)

# default attributes for netCDF variable of dataarrays
NODATA = -9999
defaultAttrsDA = {'_FillValue': NODATA, 'missing_value': NODATA}

# standard columns
basecols = ['id', 'year', 'julianday']

YEARS = range(1960, 1990)

# functions

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

#ds.pipe(update_attrs, foo='bar')


def enum(*sequential, **named):
    """ Python2 enum type """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

File = enum('UNKNOWN', 'PFT', 'MONTHLY', 'OTHER')


class LPJFile():
    """ Holds information about the LPJ file parsed """
    def __init__(self, var, inpath):
        self._filepath = os.path.join(inpath, "%s.out" % var)
        self._type = File.UNKNOWN       # type of LPJ output file
        self._years = []                # year range of data in file
        self._load_data(var, inpath)
        self.name = os.path.basename(self._filepath)

    def _load_data(self, var, inpath):
        pass

    def _determine_type():
        pass


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


def get_annual_data(var, landforms, df_frac, args, inpath='', 
                    years=[], use_month_dim=False, subset='',
                    sel_var=None):
    """ Parse variable and return DataArrays (data, total_data) """

    # derive dimensions and landform info from landforms file
    lats = landforms.coords['lat']
    lons = landforms.coords['lon']

    mapper = IndexMapper(lats, lons)

    var_str = var
    if 'sp_' in var:
        var_str = var.replace('sp_', '')

    # read either regular or compressed data files
    try:
        df = pd.read_csv(os.path.join(inpath, "%s.out" % var),
            delim_whitespace=True)
    except:
        df = pd.read_csv(os.path.join(inpath, "%s.out.gz" % var),
            delim_whitespace=True)
    
    
    # convert to dask
        
    log.debug("  Total number of data rows in file (raw data): %d" % len(df))
    
    # limit df to years
    if len(years) > 0:
        log.debug("  Limiting years.")
        df = df[ (df.Year >= years[0]) & (df.Year < years[-1])]
    if len(df) == 0:
        log.critical("Requested years not in data.")
        exit(1)
    
    log.debug("  Total number of data rows in file (year sel): %d" % len(df))

    # determine z dimension
    nyears = df.Year.max() - df.Year.min()
    print nyears
    outyears = range(df.Year.min(), df.Year.max()+1)

    def is_monthly(df):
        """Check if df is a montlhy dataframe."""
        col_names = df.columns.values
        return (('Jan' in col_names) and ('Dec' in col_names))

    def is_pft(df):
        """Check if df is a per-pft dataframe."""
        col_names = df.columns.values
        return (('Total' in col_names) and ('C3G' in col_names))

    def has_stand(df):
        """Check if df is a subpixel dataframe."""
        return 'Stand' in df.columns.values

    # TODO: cut testing runtime, remove this later
    # df = df.head(150000)


    # calc mean over patches
    if 'Patch' in df.columns.values:
        groupcols = ['Lon','Lat','Year','Stand']
        df = df.groupby(groupcols).mean().reset_index()
        del df['Patch']
        print df.head()

    if has_stand(df):
        # fix index before merge
        df.set_index(['Lon','Lat','Stand'], inplace=True)# set index for join
        df = df.merge(df_frac, left_index=True, right_index=True)
        df.reset_index(inplace=True)
    
    # calc
    if args.avg:
        log.info("Averaging years over selected timespan.")
        if 'Stand' in df.columns.values:
            groupcols = ['Lon','Lat','Stand']
        else:
            groupcols = ['Lon','Lat']
        df = df.groupby(groupcols).mean().reset_index()
        del df['Year']

    log.debug("  Total number of data rows in file (annual avg):    %d" % len(df))

    data_cols = [c for c in df.columns.values if c not in ['Lat', 'Lon', 'Year', 'Stand']]

    wtavg = lambda x: np.average(x.ix[:, data_cols], weights = x.fraction, axis=0)

    # do the stand/ lf_id avg
    if 'Stand' in df.columns.values:
        # if average is requested, do (weighted) average over the year column, too    
        #df = df.groupby(groupcols).apply(wtavg).reset_index()

        wavg = lambda x: pd.Series([sum(x[v]*x.fraction)/sum(x.fraction) for v in data_cols])
        
        df = dd.from_pandas(df, npartitions=100)
        
        df = df.groupby(groupcols).apply(wavg).compute()
        df.reset_index(inplace=True)

    log.debug("  Total number of data rows in file (final):    %d" % len(df))

    
    # determine start column position
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

    cid_start, cid_end = get_data_column_index(df)
    if is_pft(df):
        PFTS = df.columns.values.tolist()[cid_start:cid_end+1]
    else:
        PFTS = []
    
    # setup output array
    #ylen, xlen = landforms.values.shape
    xcoords = lons
    ycoords = lats
    lcoords = None  # pft axis
    mcoords = None  # extra month axis

    # time axis
    zcoords = xr.DataArray( range(len(outyears)), name='time')
    zcoords.attrs['units'] = 'yearly'

    # second time axis (monthly substeps)
    zcoords2 = None

    if use_month_dim:
        mcoords = xr.DataArray( range(12), name='month')
        mcoords.attrs['units'] = 'month'

    if is_pft(df):
        lcoords = xr.DataArray( PFTS, name='PFT')

    if is_monthly(df):
        if not args.use_month_dim:
            zcoords2 = xr.DataArray( np.linspace(0, len(outyears),
                                                 num=len(outyears)*12,
                                                 endpoint=False), name='time')

    if is_pft(df):
        # pft var
        DATA = np.ones( (len(zcoords), len(lcoords), len(ycoords), len(xcoords) ) ) * NODATA
        da_coords = [('time', zcoords), ('pft', lcoords), ('lat', ycoords), ('lon', xcoords)]
        data = xr.DataArray(DATA, name=var_str, coords=da_coords)

        DATA2 = np.ones( (len(zcoords), len(ycoords), len(xcoords) ) ) * NODATA
        da_coords2 = [('time', zcoords), ('lat', ycoords), ('lon', xcoords)]
        data2 = xr.DataArray(DATA2, name=var_str, coords=da_coords2)
    else:
        # monthly var ?
        if is_monthly(df):
            if args.use_month_dim:
                DATA = np.ones( (len(zcoords), len(mcoords),
                                 len(ycoords), len(xcoords) ) ) * NODATA
                da_coords = [('time', zcoords), ('month', mcoords),
                             ('lat', ycoords), ('lon', xcoords)]
                data = xr.DataArray(DATA, name=var_str, coords=da_coords)
            else:
                DATA = np.ones( (len(zcoords2), len(ycoords), len(xcoords) ) ) * NODATA
                da_coords = [('time_m', zcoords2), ('lat', ycoords), ('lon', xcoords)]
                data = xr.DataArray(DATA, name=var_str, coords=da_coords)

        else:
            # do we have a file with individual var columns and a selection
            # was requested?
            if sel_var is not None:
                DATA = np.ones( (len(zcoords), len(ycoords), len(xcoords) ) ) * NODATA
                da_coords = [('time', zcoords), ('lat', ycoords), ('lon', xcoords)]                
                data = xr.DataArray(DATA, name=var_str, coords=da_coords)

        data2 = None

    # add variable attributes
    data.attrs['units'] = '-'
    data.attrs['missing_value'] = NODATA
    data.attrs['_FillValue'] = NODATA
    # optional (second) DataArray for Total (PFT) data

    # landform lookup
    # TODO: optimize that we only need to do this once
    #       maybe wrap into a class (?)

    # processing for stand/ patch output
    if has_stand(df) and 'lf_id' in landforms.data_vars:

        lookup = {}
        
        reffile, lf_fraction_var = args.refinfo
        lfids = landforms['lf_id'].values
        # add stand weight column (using lanform fraction)
        df['lffrac'] = -1.0

        done_coords = []

        if subset=='north':
            # aspect bit 1
            subset_stands = [x for x in lfids if str(int(x))[-1] == '1']
        elif subset == 'south':
            # aspect bit 3
            subset_stands = [x for x in lfids if str(int(x))[-1] == '3']
        else:
            subset_stands = lfids

        for rx, row in df.iterrows():
            # get coord and add to lookup table, do not run if coord is present
            jx, ix = mapper(row.Lat, row.Lon)
            stand = int(row.Stand)

            if (jx, ix, stand) not in lookup.keys():
                if (jx, ix) not in done_coords:
                    done_coords.append((jx, ix))
                frac = landforms[lf_fraction_var][:, jx, ix].sel(lf_id=stand)
                lookup[(jx, ix, stand)] = frac.values * 0.01

            # assign value (-1 values will be removed in the next step)
            if stand in subset_stands:
                df.loc[rx, 'lffrac'] = lookup[(jx, ix, stand)]


        # drop rows if they are not in subset_stands
        df = df[df['lffrac'] != -1]

        # aggregate (subset) of landforms

        def weighted(x, cols, w="lffrac"):
            """Weighted-average over a group of rows."""
            return pd.Series(np.average(x[cols], weights=x[w], axis=0), cols)

        groupcols = ['Lon','Lat','Year']
        exclcols = groupcols + ['lffrac']
        datacols = [x for x in list(df.columns.values) if x not in exclcols]
        df = df.groupby(groupcols).apply(weighted, datacols).reset_index()

    for _, row in df.iterrows():
        jx, ix = mapper(row.Lat, row.Lon)
        rdata = row[cid_start: cid_end+1]
        if is_monthly(df):
            if args.use_month_dim:
                zpos = outyears.index(row['Year'])
                data[zpos, :, jx, ix] = rdata
            else:
                zpos = outyears.index(row['Year'])*12
                data[zpos:zpos+12, jx, ix] = rdata
        elif is_pft(df):
            zpos = outyears.index(row['Year'])
            data[zpos,:,jx,ix] = rdata
            data2[zpos,jx,ix] = row['Total']
        elif sel_var in rdata.index:
            zpos = outyears.index(row['Year'])
            data[zpos,jx,ix] = rdata[sel_var]
        else:
            zpos = outyears.index(row['Year'])
            if len(rdata) == 1:
                data[zpos,jx,ix] = rdata
            else:
                log.critical("Mixed output file not handled yet:")
                log.critical("  This should create a single DataArray")
                log.critical("  for each column (with filename prefix")
                log.critical(var)
                exit(1)

    # add attributes for total variable
    if type(data2) == xr.DataArray:
        data2.attrs['units'] = '-'
        data2.attrs['missing_value'] = NODATA
        data2.attrs['_FillValue'] = NODATA

    return data, data2


def main():
    """Main routine."""

    # INPUT
    parser, args = cli()
    cfg = get_config(cfgFile=args.config)

    if args.storeconfig:
        set_config(cfg)

    def use_cli_refdata():
        return args.refinfo is not None

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
                df = da.sel(lf_id=lf_id).to_pandas().stack().reset_index()
                df.columns = ['Lat','Lon','fraction']
                if len(df) > 0:
                    df['Stand'] = lf_id.values
                    dfs.append(df)

            df = pd.concat(dfs)
            df.set_index(['Lon', 'Lat', 'Stand'], inplace=True)
            print df.head()
        return df

    if use_cli_refdata():
        rdata = args.refinfo
        reffile, refvar = rdata

        if os.path.isfile(reffile):
            with xr.open_dataset(reffile).load() as ds_lf:
                if refvar is not None:
                    df_frac = get_fractions(ds_lf, refvar)
                else:
                    refvar = 'fraction'
                    df_frac = get_fractions(ds_lf, refvar)
    else:
        log.error("We currently need a refdata file!")
        exit(-1)
                        
    #else:
    #    # get domain info from conf file
    #    # build a dummy array with the specified coords but without
    #    # landform info
    #    if 'refdata' in cfg.keys():
    #        lats, lons = _read_refdata_info(cfg)
    #        dummy = xr.DataArray(np.ones((len(lats), len(lons))), 
    #                               coords=[('lat', lats ),('lon', lons)])
    #        landforms = xr.Dataset()
    #        landforms['dummy'] = dummy
    #    else:
    #        log.critical("No refdata section in conf file found nor -r flag specified.\n")
    #        log.critical(parser.print_help())
    #        exit(1)
            
    # derive data from cli or config file
    global_info = _read_global_info(cfg)

    inpath = args.indir
    
    if not os.path.isdir(inpath):
        log.critical("Specified input directory does not exist.\n")
        log.critical(parser.print_help())
        exit(1)
    
    outname = args.outname #'lpjguess_spatial_200_mdim.nc'
    use_month_dim = args.use_month_dim

    data_files = [x[0] for x in _read_data_info(cfg)]

    log.debug("TODO: We need to produce proper time attributes")

    # output netcdf file
    ds = xr.Dataset()

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
        if is_subpixel_file(file) and args.north_south:
            log.debug("  Calculating south, north and full landform average.")
            sset_list = ['north', 'south', '']
        else:
            log.debug("  Calculating full landform average.")
            sset_list = ['']

        for sset in sset_list:
            suffix = ''
            if sset != '':
                suffix = '_%s' % sset

            # use filename as variable, also selected column/ var if requested
            if len(named_vars.keys()) > 0:
                sel_vars = named_vars.keys()
                for var in sel_vars:
                    da1, _ = get_annual_data(file, ds_lf, df_frac, args,
                                            inpath=inpath,
                                            use_month_dim=use_month_dim,
                                            subset=sset,
                                            sel_var=var)
                    # rename variable
                    ds[named_vars[var] + suffix] = da1
            else:
                da1, da2 = get_annual_data(file, ds_lf, df_frac, args,
                                        inpath=inpath,
                                        use_month_dim=use_month_dim,
                                        subset=sset)
            
                ds[file + suffix] = da1
                # if we have a second (total) dataarray
                if type(da2) == xr.DataArray:
                    ds['tot_%s%s' % (file, suffix)] = da2


    # all potential coordinate variables
    all_c_vars = ['time','time_m','month','pft','lat','lon']
    c_vars = [x for x in all_c_vars if x in ds.coords.keys()]

    # order variables (coords first, then by name)
    d_vars = [x for x in sorted([x for x in ds.data_vars]) if x not in c_vars]

    ds[c_vars + d_vars].to_netcdf(outname, format='NETCDF4_CLASSIC', unlimited_dims='time')

    log.debug("Done.")
