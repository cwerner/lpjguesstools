"""FILE lgt_biomize.main.py

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

#from .cli import cli
#from .extra import set_config, get_config, parse_config #, RefDataBuilder

__version__ = "0.0.1"

log = logging.getLogger(__name__)

# default attributes for netCDF variable of dataarrays
NODATA = -9999

def enum(*sequential, **named):
    enums = OrderedDict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['content'] = [reverse[x] for x in enums.values()]
    enums['reverse_mapping'] = reverse
    return type('Enum', (), enums)



class Biomization( object ):
    def __init__(self, biome_type='CHILE_ES_NEW'):
        self._pfts = None
        self._btype = biome_type

        if biome_type not in ['SMITH', 'CHILE', 'CHILE_ES', 'CHILE_ES_NEW']:
            log.error("Biome classification not valid.")
            exit(-1)

        if self._btype == 'CHILE_ES_NEW':
            from biomes_earthshape import biomes
            from biomes_earthshape import classification
            self._biomes = biomes
            self._lut = classification

    def __call__(self, da):
        """Execute classification on supplied DataArray."""
        if self._pfts is None:
            self._pfts = da['pft'].values
        return self._classify(da)

    def _classify(self, da, class_var='sp_lai'):
        """Call classification."""
        return self._lut(self._biomes, self._pfts, da)

    def biomes(self):
        """Return the biome names of this classification."""
        return self._biomes.content

    def pfts(self):
        """Return the pft names used in this biome classification."""
        return self._pfts



def main(cfg):
    """Main Script."""    

    # load data
    with xr.open_dataset(cfg.INFILE) as ds:
        da_lai  = ds['sp_lai'].load()
        da_frac = ds['fraction'].load()
        lf_ids  = ds['lf_id'].values
        pfts    = ds['pft'].values


    # init the classifier
    biomizer = Biomization(biome_type=cfg.CLASSIFICATION)

    # destination xr.DataArray
    if cfg.SMODE:
        # target layout dims: 2D (lf, time)
        log.info("Single site mode.")
        log.debug("  attention: file time dim required (no time_m !)")
        time_len = len(da_lai['time'])
        da_biomes = xr.DataArray(np.ones( (len(lf_ids), time_len) ) * np.nan, 
                            coords=[da_lai.coords['lf_id'], da_lai.coords['time']],
                            dims=['lf_id', 'time'], name='biome_lf')
        
        # get coordinate from global attribute
        if ('site_lon' in ds.attrs.keys()) and ('site_lat' in ds.attrs.keys()):
            site_lon = float(ds.attrs['site_lon'])
            site_lat = float(ds.attrs['site_lat'])
            #ix = list(ds.coords['lon'].values).index(site_lon)
            #jx = list(ds.coords['lat'].values).index(site_lat)
        else:
            log.error('INFILE does not contain global attributes: site_lon/ site_lat.')
            exit(-1)
    
        for tx, _year in enumerate(da_lai.coords['time']):
            if tx % 1000 == 0: print tx
            for zx, _lf_id in enumerate(lf_ids):
                if da_frac[zx] > 0:
                    b = biomizer( da_lai.sel(time=_year, lf_id=_lf_id) )
                    da_biomes[zx, tx] = b        

    else:
        # target layout dims: 3D (lf, lat, lon)
        log.info("Spatial mode.")
        da_biomes = xr.ones_like(da_frac).rename('biome_lf') * np.nan

        # loop over lf_ids and compute biomization in this layer
        for jx, _lat in enumerate(da_frac.coords['lat']):
            print jx
            for ix, _lon in enumerate(da_frac.coords['lon']):
                for zx, _lf_id in enumerate(lf_ids):
                    #if da_frac.sel(lat=_lat, lon=_lon, lf_id=_lf_id) > 0:
                    if da_frac[zx, jx, ix] > 0:
                        b = biomizer( da_lai.sel(lat=_lat, lon=_lon, lf_id=_lf_id) )
                        da_biomes[zx, jx, ix] = b



    # switch to remove aspect landforms from dataset
    drop_aspect_lf = False
    if drop_aspect_lf:                
        # kick out all aspect landforms for now
        log.warn("Deleting all sloped landforms (aspect in 1,2,3,4) for now!")
        da_biomes = da_biomes.sel(lf_id=[x for x in da_biomes.lf_id.values if x%10==0])

    da_biomes.attrs['_FillValue'] = NODATA
    da_biomes.attrs['units'] = 'biome_id'
    ds = da_biomes.to_dataset()
    ds['fraction'] = da_frac
    
    # sum fractions containing a given biome 
    # returns: stacked array (z: axis biome id (0...8)
    agg = []
    for b, bname in enumerate(biomizer.biomes()):
        fr=ds['fraction'].where(ds['biome_lf']==b)
        u=fr.sum(dim='lf_id')
        agg.append(u)

    da2 = xr.concat(agg, dim='biome').argmax(dim='biome').where(ds['fraction'].sum(dim='lf_id')>0)
    
    ds['biome'] = da2
    ds['biome'].attrs['_FillValue'] = NODATA
    ds['biome'].attrs['units'] = 'biome_id'

    ds.to_netcdf(cfg.OUTFILE[:-3] + '_biome.nc', format='NETCDF4_CLASSIC')

    # --- section 2 ---
    # aggregate to elevation levels
    # this is only for spatial mode
    
    if cfg.SMODE == False:
    
        ds['elevclass'] = ds['lf_id'].astype('int')/100

        eles = np.unique(ds['elevclass'])

        blank = np.ones((len(eles), len(ds.coords['lat']), len(ds.coords['lon']))) * np.nan
        ds['elevfraction'] = xr.DataArray(blank, coords=[('ele', eles), ('lat', ds.coords['lat']), ('lon', ds.coords['lon'])], dims=['ele', 'lat', 'lon'])

        # aggregate to lat/ele levels
        ele_L = []
        # values of lf_ids of from this ele class
        A = list(ds['fraction'].groupby( ds['elevclass'] ))
        B = list(ds['biome_lf'].groupby( ds['elevclass'] ))

        # iterate over elevation levels
        for e_id in range(len(A)):
            _, e_frac  = A[e_id]
            _, e_biome = B[e_id]

            agg = []
            for b, bname in enumerate(biomizer.biomes()):
                fr=e_frac.where(e_biome==b)
                u=fr.sum(dim=['lf_id', 'lon'])
                u=u.where(u>0)
                
                # add a dummy layer to differentiate all-nan from biome 0
                if b==0:
                    agg.append(xr.ones_like(u) * -1)
                agg.append(u)

            e_da1 = xr.concat(agg, dim='biome').argmax(dim='biome')
            
            # if argmax returned 0 it was all-nan: skip cell
            # also, due to the dummy layer we need to subtract 1 to get true 
            # biome_id number
            e_da1 = e_da1.where(e_da1 > 0) - 1 
            e_da2 = e_da1.where(e_frac.sum(dim=['lf_id', 'lon'])>0)
            ele_L.append(e_da2)

        da_ele = xr.concat(ele_L, pd.Index( np.arange(100, 30*200+100, 200, dtype=np.int32), name='ele')).astype(np.int32)
        da_ele.name='biome'
        da_ele.attrs['units'] = 'biome_id'

        da_ele['ele'] = da_ele['ele'].astype(np.int32)
        da_ele['ele'].attrs['axis'] = 'X'
        da_ele['lat'].attrs['axis'] = 'Y'
        da_ele.T.to_netcdf(cfg.OUTFILE[:-3] + '_biome_latele.nc', format='NETCDF4_CLASSIC')

