# -*- coding: utf-8 -*-
#
# ldndc2nc.py
# ==================
#
# Christian Werner
# christian.werner@senckenberg.de
"""ldndc2nc.ldndc2nc: provides entry point main()."""

import calendar
from collections import OrderedDict
import datetime as dt
import glob
import logging
import os
import re

import numpy as np
import pandas as pd
import xarray as xr

from .cli import cli
from .extra import set_config, get_config, parse_config, RefDataBuilder

__version__ = "0.0.3"

log = logging.getLogger(__name__)

# default attributes for netCDF variable of dataarrays
NODATA = -9999
defaultAttrsDA = {'_FillValue': NODATA, 'missing_value': NODATA}

# standard columns
basecols = ['id', 'year', 'julianday']


# functions
def _split_colname(colname):
    """ split colname into varname and units

        :param str colname: original ldndc outputfile colname
        :return: varname and unit
        :rtype: tuple
    """
    out = (colname, "unknown")
    if '[' in colname:
        name, var_units = colname.split('[')
        units = var_units[:-1]
        out = (name, units)
    return out


def _daterange(start_date, end_date):
    """ create timeseries

        :param str start_date: start date
        :param str end_date: start date
        :return: list of dates
        :rtype: iterator
    """
    for n in range(int((end_date - start_date).days)):
        yield start_date + dt.timedelta(n)


def _gapfill_df(df, dates, ids):
    """ gap-fill data.frame """
    basecoldata = [(0, d.year, d.timetuple().tm_yday) for d in dates]
    df_ref_all = []
    for id in ids:
        df_ref = pd.DataFrame(basecoldata, columns=basecols)
        df_ref.id = id
        df_ref_all.append(df_ref)
    df_template = pd.concat(df_ref_all, axis=0)
    df = pd.merge(df_template, df, how='left', on=basecols)
    df = df.fillna(0.0)
    df.reset_index()
    return df


def _ndays(yr):
    """ return the number of days in year """
    ndays = 365
    if calendar.isleap(yr):
        ndays = 366
    return ndays


def _is_composite_var(v):
    return type(v) == tuple


def _all_items_identical(x):
    return x.count(x[0]) == len(x)


def _build_id_lut(array):
    """ create lookup table to query cellid by i,j value """
    Dlut = {}
    for j in range(len(array)):
        for i in range(len(array[0])):
            if not np.isnan(array[j, i]):
                Dlut[int(array[j, i])] = (len(array) - j - 1, i)
    return Dlut


def _extract_fileno(fname):
    """ extract file iterator

        :param str fname: ldndc txt filename
        :return: file number
        :rtype: int

        example: GLOBAL_002_soilchemistry-daily.txt -> 002 -> 2
    """
    fname = os.path.basename(fname)
    fileno = 0
    # find fileno in string (must be 2-6 digits long)
    x = re.findall(r"[0-9]{2,6}(?![0-9])", fname)
    if len(x) == 0:
        pass
    elif len(x) == 1:
        fileno = int(x[0])
    else:
        log.critical("Multiple matches! fname: %s" % fname)
        exit(1)
    return fileno


def _select_files(inpath, ldndc_file_type, limiter=""):
    """ find all ldndc outfiles of given type from inpath (limit using limiter)

        :param str inpath: path where files are located
        :param str ldndc_file_type: LandscapeDNDC txt filename pattern
                   (i.e. soilchemistry-daily.txt)
        :param str limiter: (optional) limit selection using this expression
        :return: list of matching LandscapeDNDC txt files in indir
        :rtype: list
    """
    infile_pattern = os.path.join(inpath, "*" + ldndc_file_type)
    infiles = glob.glob(infile_pattern)

    if limiter != "":
        infiles = [x for x in infiles if limiter in os.path.basename(x)]

    infiles.sort()

    if len(infiles) == 0:
        msg = "No LandscapeDNDC input files of type <%s>\n" % ldndc_file_type
        msg += "Input dir:    %s\n" % inpath
        msg += "Pattern used: %s" % infile_pattern
        if limiter != "":
            msg += "\nFilter used:  %s" % limiter
        log.critical(msg)
        exit(1)

    return infiles


def _limit_df_years(years, df, yearcol='year'):
    """ limit data.frame to specified years """
    if (years[-1] - years[0] == len(years) - 1) and (len(years) > 1):
        df = df[(df[yearcol] >= years[0]) & (df[yearcol] <= years[-1])]
    else:
        df = df[df[yearcol].isin(years)]
    if len(df) == 0:
        if len(years) == 1:
            log.critical('Year %d not in data' % years[0])
        else:
            log.critical('Year range %d-%d not in data' %
                         (years[0], years[-1]))
        exit(1)
    df = df.sort_values(by=basecols)
    return df


def _read_global_info(cfg):
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


def read_ldndc_txt(inpath, varData, years, limiter=''):
    """ parse ldndc txt output files and return dataframe """

    ldndc_file_types = varData.keys()

    varnames = []  # (updated) column names
    Dids = {}  # file ids

    df_all = []

    for ldndc_file_type in ldndc_file_types:

        dfs = []
        datacols = []

        infiles = _select_files(inpath, ldndc_file_type, limiter=limiter)

        # special treatment for tuple entries in varData
        for v in varData[ldndc_file_type]:
            if _is_composite_var(v):
                varnames.append(v[0])
                datacols += v[1]
            else:
                varnames.append(v)
                datacols.append(v)

        # iterate over all files of one ldndc file type
        for fcnt, fname in enumerate(infiles):
            fno = _extract_fileno(fname)
            df = pd.read_csv(fname,
                             delim_whitespace=True,
                             error_bad_lines=False,
                             usecols=basecols + datacols)

            Dids.setdefault(fno, sorted(list(set(df['id']))))

            df = _limit_df_years(years, df)

            dates = list(_daterange(dt.date(years[0], 1, 1),
                         dt.date(years[-1]+1, 1, 1)))

            len_full_df = len(Dids[fno]) * len(dates)
            len_this_df = len(df)

            if len_this_df < len_full_df:
                df = _gapfill_df(df, dates, Dids[fno])

            df = df.sort_values(by=basecols)
            dfs.append(df)

        # we don't have any dataframes, return
        # TODO: the control flow here should be more obvious
        if len(dfs) == 0:
            log.warn("No data.frame filetype %s!" % ldndc_file_type)
            continue

        df = pd.concat(dfs, axis=0)
        df.set_index(basecols, inplace=True)

        # sum columns if this was requested in the conf file
        for v in varData[ldndc_file_type]:
            if _is_composite_var(v):
                new_colname, src_colnames = v
                drop_colnames = []

                df[new_colname] = df[src_colnames].sum(axis=1)

                # drop original columns if they are not explicitly requested
                for v2 in varData[ldndc_file_type]:
                    if not _is_composite_var(v2):
                        if v2 in src_colnames:
                            drop_colnames.append(v2)

                df.drop(drop_colnames, axis=1)

        # TODO check why we have this line
        df = df[~df.index.duplicated(keep='first')]
        df_all.append(df)

    # check if all tables have the same number of rows
    if _all_items_identical([len(x) for x in df_all]):
        log.debug("All data.frames have the same length (n=%d)" %
                  len(df_all[0]))
    else:
        log.debug("Rows differ in data.frames: %s" %
                  ''.join([str(len(x)) for x in df_all]))

    df = pd.concat(df_all, axis=1)
    df.reset_index(inplace=True)

    return (varnames, df)


def main():
    # parse args
    args = cli()

    # read config
    cfg = get_config(args.config)

    # write config
    if args.storeconfig:
        set_config(cfg)

    # read or build refdata array
    def use_cli_refdata():
        return args.refinfo is not None

    if use_cli_refdata():
        reffile, refvar = args.refinfo
        if os.path.isfile(reffile):
            with (xr.open_dataset(reffile)) as refnc:
                if refvar not in refnc.data_vars:
                    log.critical("Var <%s> not in %s" % (refvar, reffile))
                    exit(1)
                cell_ids = np.flipud(refnc[refvar].values)
                lats = refnc['lat'].values
                lons = refnc['lon'].values
        else:
            log.critical("Specified reffile %s not found" % reffile)
            exit(1)
    else:
        rdb = RefDataBuilder(cfg)
        cell_ids, lats, lons = rdb.build()

    # get general info
    global_info = _read_global_info(cfg)

    # create lut for fast id-i,j matching
    Dlut = _build_id_lut(cell_ids)

    # read source output from ldndc
    varinfos, df = read_ldndc_txt(args.indir,
                                  cfg['variables'],
                                  args.years,
                                  limiter=args.limiter)

    # process data in yearly chunks
    ds_all = []
    for yr, yr_group in df.groupby('year'):

        # loop over variables and add those the netcdf file
        times = pd.date_range('%s-01-01' % yr,
                              freq='D',
                              periods=_ndays(yr),
                              tz=None)

        new_shape = (_ndays(yr),) + cell_ids.shape

        # init datasets
        ds = xr.Dataset(attrs=global_info)
        for varinfo in varinfos:
            name, units = _split_colname(varinfo)
            varAttrs = defaultAttrsDA.copy()
            varAttrs.update({'units': units})
            blank_array = np.ma.array(np.ones(new_shape)*NODATA, mask=True)
            ds[name] = xr.DataArray(blank_array,
                                    name=name,
                                    coords=[('time', times),
                                            ('lat', lats),
                                            ('lon', lons)],
                                    attrs=varAttrs,
                                    encoding={'complevel': 5,
                                              'zlib': True,
                                              'chunksizes': (10, 40, 20),
                                              'shuffle': True})

        # iterate over cellids and variables
        for id, id_group in yr_group.groupby('id'):
            jdx, idx = Dlut[id]
            for varinfo in varinfos:
                name, _ = _split_colname(varinfo)
                if len(id_group[varinfo]) < len(times):
                    missingvals = _ndays(yr) - len(id_group[varinfo])
                    dslize = np.concatenate(id_group[varinfo],
                                            np.asarray([NODATA] * missingvals))
                    log.warn("Data length encountered shorter than expected!")
                else:
                    dslize = id_group[varinfo]
                ds[name][:, jdx, idx] = dslize

        if args.split:
            outfilename = args.outfile[:-3] + '_%d' % yr + '.nc'
            ds.to_netcdf(
                os.path.join(args.outdir, outfilename),
                format='NETCDF4_CLASSIC')
            ds.close()
        else:
            ds_all.append(ds)

    if not args.split:
        ds = xr.concat(ds_all, dim='time')
        ds.to_netcdf(
            os.path.join(args.outdir, args.outfile),
            format='NETCDF4_CLASSIC')
        ds.close()
