# -*- coding: utf-8 -*-
#
# _lpjguess2nc.py
# ====================
#
# Christian Werner
# christian.werner@senckenberg.de
"""lpjguess2nc.lpjguess2nc: provides entry point main()."""

import logging
import numpy as np
import pandas as pd
import xarray as xr
import os
import sys

#from .cli import cli
#from .extra import set_config, get_config, parse_config, RefDataBuilder

__version__ = "0.0.1"

log = logging.getLogger(__name__)

# default attributes for netCDF variable of dataarrays
NODATA = -9999
defaultAttrsDA = {'_FillValue': NODATA, 'missing_value': NODATA}

# standard columns
basecols = ['id', 'year', 'julianday']



YEARS = range(1950, 1990)

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
                print 'Only 1D numpy arrays supported'
                print lats
                print lons
                exit()
        elif ((type(lats) == xr.DataArray) and (type(lons) == xr.DataArray)):
            self._lons = lons.values.tolist()
            self._lats = lats.values.tolist()
        else:
            print 'Only 1D numpy arrays or lists supported'
            print lats
            print lons
            exit()

    def __call__(self, lat, lon):
        if ((lon >= self._lons[0]) and (lon <= self._lons[-1])):
            ix = self._lons.index(lon)
        else:
            print 'Lon %s out of range [%f...%f]' % (lon, 
                                                     self._lons[0],
                                                     self._lons[-1])
            exit()
            
        if ((lat >= self._lats[0]) and (lat <= self._lats[-1])):
            jx = self._lats.index(lat)
        else:
            print 'Lon %s out of range [%f...%f]' % (lon, 
                                                     self._lats[0],
                                                     self._lats[-1])
            exit()
        return (jx, ix)


def get_annual_data(var, refarray, inpath='', years=[], use_month_dim=False):
    """ Parse variable and return DataArrays (data, total_data) """

    lats = refarray.coords['lat']
    lons = refarray.coords['lon']
    
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
        
    # limit df to years
    if len(years) > 0:
        print years
        df = df[ (df.Year >= years[0]) & (df.Year < years[-1])]
    if len(df) == 0:
        print 'Requested years not in data.'
        exit()
    
    # determine z dimension
    nyears = max(df.Year) - min(df.Year)
    outyears = range(min(df.Year), max(df.Year)+1)
    
    ismonthly = False
    ispft     = False
    PFTS = []

    # determine start column position
    if 'Patch' in df.columns.values:
        cid_start = df.columns.values.tolist().index('Patch') + 1
    elif 'Stand' in df.columns.values:
        cid_start = df.columns.values.tolist().index('Stand') + 1
    else:
        cid_start = df.columns.values.tolist().index('Year') + 1

    # determine end column position
    # default: last column
    cid_end   = len(df.columns.values.tolist()) 
    
    # determine if monthly file
    if 'Jan' in df.columns.values:
        ismonthly = True
        # time axis years x12
        cid_end   = df.columns.values.tolist().index('Dec')
        
    # determine if PFT file
    if (('Total' in df.columns.values) and ('C3G' in df.columns.values)):
        ispft = True
        # time axis years 
        # level axis PFT S4

        cid_end = df.columns.values.tolist().index('Total') - 1
        PFTS = df.columns.tolist()[cid_start:cid_end+1]
        
    # setup output array
    ylen, xlen = refarray.values.shape
    xcoords = refarray.coords['lon']
    ycoords = refarray.coords['lat']
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
    
    if ispft:
        lcoords = xr.DataArray( PFTS, name='PFT')

    if ismonthly:
        if not use_month_dim:
            zcoords2 = xr.DataArray( np.linspace(0, len(outyears), 
                                                 num=len(outyears)*12, 
                                                 endpoint=False), name='time')
            

    if ispft is False:
        # monthly var ?
        if ismonthly:
            if use_month_dim:
                DATA = np.ones( (len(zcoords), len(mcoords), 
                                 len(ycoords), len(xcoords) ) ) * NODATA       
                da_coords = [('time', zcoords), ('month', mcoords), 
                             ('lat', ycoords), ('lon', xcoords)]
                data = xr.DataArray(DATA, name=var_str, coords=da_coords)            
            else:
                DATA = np.ones( (len(zcoords2), len(ycoords), len(xcoords) ) ) * NODATA       
                da_coords = [('time_m', zcoords2), ('lat', ycoords), ('lon', xcoords)]
                data = xr.DataArray(DATA, name=var_str, coords=da_coords)

        data2 = None

    else:
        # pft var
        DATA = np.ones( (len(zcoords), len(lcoords), len(ycoords), len(xcoords) ) ) * NODATA
        da_coords = [('time', zcoords), ('pft', lcoords), ('lat', ycoords), ('lon', xcoords)]
        data = xr.DataArray(DATA, name=var_str, coords=da_coords)
            
        DATA2 = np.ones( (len(zcoords), len(ycoords), len(xcoords) ) ) * NODATA
        da_coords2 = [('time', zcoords), ('lat', ycoords), ('lon', xcoords)]
        data2 = xr.DataArray(DATA2, name=var_str, coords=da_coords2)


    data.attrs['units'] = '-'
    data.attrs['missing_value'] = NODATA
    data.attrs['_FillValue'] = NODATA    
    # optional (second) DataArray for Total (PFT) data
    
    for _, row in df.iterrows():
        jx, ix = mapper(row.Lat, row.Lon)
        rdata = row[cid_start: cid_end+1]

        if ismonthly:
            if use_month_dim:
                zpos = outyears.index(row['Year'])
                data[zpos, :, jx, ix] = rdata
            else:
                zpos = outyears.index(row['Year'])*12
                data[zpos:zpos+12, jx, ix] = rdata 
        elif ispft:
            zpos = outyears.index(row['Year'])            
            data[zpos,:,jx,ix] = rdata
            data2[zpos,jx,ix] = row['Total']
        else:
            zpos = outyears.index(row['Year'])
            if len(rdata) == 1:        
                data[zpos,jx,ix] = rdata
            else:
                print 'Mixed output file not handled yet.'
                print '  this should create a single DataArray'
                print '  for each column (with filename prefix)'
                print var
                exit()
                
    if type(data2) == xr.DataArray:
        data2.attrs['units'] = '-'
        data2.attrs['missing_value'] = NODATA
        data2.attrs['_FillValue'] = NODATA

    return data, data2


def main():
    """Main routine."""
    
    # INPUT
    inpath = sys.argv[1]
    outname = sys.argv[2] #'lpjguess_spatial_200_mdim.nc'
    use_month_dim = True
        
    print 'TODO: We need to produce proper time attributes'
    reffile = os.path.join(inpath, '..', 'lfdata', 'sites_2d.nc')
    print reffile
    refarray = xr.open_dataset(reffile)['elevation']
    
    vars = ['lai','mlai','fpc','cmass','mrunoff','mwcont_upper','mnpp','mgpp'] #'dens', 'mnpp','mgpp','firert', 'cmass','mlai','fpc','lai']
    ds = xr.Dataset()
    
    for var in vars:
        print 'parsing %s ...' % var
        da1, da2 = get_annual_data(var, refarray, 
                                   inpath=inpath, 
                                   use_month_dim=use_month_dim) #, years=YEARS)
        
        ds[var] = da1
        if type(da2) == xr.DataArray:
            ds['tot_%s' % var] = da2
    
    if use_month_dim:
        c_vars = ['time','month','pft','lat','lon']
    else:
        c_vars = ['time','time_m','pft','lat','lon']
    
    d_vars = [x for x in sorted([x for x in ds.data_vars]) if x not in c_vars]
    print c_vars
    print d_vars
    
    print c_vars + d_vars
    
    ds[c_vars + d_vars].to_netcdf(outname, format='NETCDF4_CLASSIC', unlimited_dims='time')

