# -*- coding: utf-8 -*-
#
# lgt_convert.py
# ==============
#
# Christian Werner
# christian.werner@senckenberg.de
"""lpjguesstools.lgt_biomize: provides entry point main()."""

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

def fr(v, tot):
    """Calculate fraction."""
    if tot > 0:
        return v/float(tot)
    else:
        return 0


def LUT_ES_CHILE_NEW(biome, all_pfts, data):
    """EarthShape Custom Biomization"""

    #biome  = enum('TeBEF',  #= 1,   # BIOME 1 - Temperate Broadleaved Evergreen Forest
    #              'TeMF',   #= 2,   # BIOME 2 - Temperate Mixed Forest
    #              'SclWF',  #= 3,   # BIOME 3 - Sclerophyllous Woodland/ Forest
    #              'TeW',    #= 4    Temperate Woodland
    #              'ArT',    #= 7,   # BIOME 6 - Arctic/alpine Tundra
    #              'Gr',     #= 9,   # BIOME 8 - Grass Tundra
    #              'AShr',   #= 10,  # BIOME 9 - Arid Shrubland
    #              'D')      #= 10}) # BIOME 10 - Desert   


    def check_if_pfts_dominate(pfts):
        """ Check if any of the given PFTs dominate the tree PFTs """

        _TreePFT = [p for p in pfts if p in TreePFT]
        if len(_TreePFT) == 0:
            return False

        _lai = np.array([data.sel(pft=p).values for p in _TreePFT])
        ix = np.argmax(_lai)
        if _TreePFT[ix] in pfts:
            return True
        else:
            return False

    def fr(pfts, woody=False):
        """ Fraction of the sum of given PFTs with reference to woody or total LAI """

        # reference: woody or total
        if woody:
            lai_ref = lai_w
        else:
            lai_ref = lai_tot

        if type(pfts) == list:
            _pfts = [p for p in pfts if p in all_pfts]
            if lai_ref > 0:
                return sum([data.sel(pft=p).values for p in _pfts]) / lai_ref
            else:
                return 0.0
        else:
            if pfts in all_pfts:
                return data.sel(pft=pfts).values / lai_ref
            else:
                return 0.0


    # Categorize PFTs

    # classify shrubs
    SclShrubPFT = ['TeE_s', 'TeR_s']        # sclerophyllous shrubs
    CldShrubPFT = ['BE_s', 'BS_s']          # cold zone shrubs
    ShrubPFT = SclShrubPFT + CldShrubPFT    # total shrubs

    # grass PFTs
    GrassPFT = ['C3G']

    # tree PFTs
    TreePFT  = [x for x in all_pfts if x not in GrassPFT+ShrubPFT]
    
    # woody PFTs
    WoodyPFT = [x for x in all_pfts if x not in GrassPFT]


    # Special tree PFT groups
    tebe_trees = ['TeBE_tm', 'TeBE_itm']    # temperate broadleaf evergreen
    tebs_trees = ['TeBS_tm', 'TeBS_itm']    # temperate decidious (Notofagus)

    # Abbreviations
    max_t   = check_if_pfts_dominate    # shorthand for the check_if_pfts_dominate function

    # Aggregated LAI values
    lai_t   = sum([data.sel(pft=p).values for p in TreePFT])   # total lai of trees
    lai_s   = sum([data.sel(pft=p).values for p in ShrubPFT])  # total lai of shrubs
    lai_w   = sum([data.sel(pft=p).values for p in WoodyPFT])  # total lai of woody (tree+shrub)
    lai_g   = sum([data.sel(pft=p).values for p in GrassPFT])  # total lai of grasses
    lai_tot = lai_t + lai_w                                 # lai of all PFTs

    def sum_lai(pfts):
        return sum([data.sel(pft=p).values for p in pfts])


    # Thresholds
    FT = 2.5   # forest/woodland threshold
    WT = 1.0  # woodland/grassland threshold

    # Decision tree (new)

    # temperate rainforest (TeBEF) - a forest, TeBE >30%, TeBE dominate trees
    if (lai_t >= FT and fr(tebe_trees) >= 0.3 and max_t(tebe_trees)):
        b = biome.TeBEF   # temperate rainforest

    # temperate mixed forest (TeMF) - a forest, TeBE <30%
    elif (lai_t >= WT and fr(tebe_trees) < 0.3 and fr(tebs_trees) <= 0.15):
        b = biome.TeMF   # temperate mixed forest

    # notofagus decidious woodland (TeBS) - a forest/ woodland, TeBS > 30%, TeBS dominante trees
    elif (lai_w >= WT and fr(tebe_trees) < 0.3 and fr(tebs_trees) > 0.15): # and max_t(tebs_trees)):
        b = biome.TeBS    # decidious forest (Notofagus)

    # high-altitude forests TeBS dominante trees
    elif (lai_w >= WT * 0.5 and fr(['TeBE_itscl']+SclShrubPFT, woody=True) < 0.5) and sum_lai(CldShrubPFT) < 0.3:   # and max_t(tebs_trees)):
        b = biome.TeMF    # high-alt mixed

    # sclerophyllous (SclWF) - woodland, sclerophyllous PFTs >= 50% of woody, cold shrubs < 0.05 LAI
    elif (lai_w >= WT * 0.5 and fr(['TeBE_itscl']+SclShrubPFT, woody=True) >= 0.5): # and sum_lai(CldShrubPFT) < 0.05:
        b = biome.SclWF   # sclerophyllous woodland

    # matorral (Mat) - less than a woodland, sclerophyllous PFTs >= 50%
    elif (lai_w < WT * 0.5 and fr(['TeBE_itscl']+SclShrubPFT, woody=True) >= 0.5 and lai_w > 0.1) :
        b = biome.Mat     # matorral

    # high-altitude forests TeBS dominante trees
    elif (lai_w < WT * 0.5 and fr(['TeBE_itscl']+SclShrubPFT, woody=True) < 0.5 and lai_w > 0.1) :
        b = biome.TeMF     # high-alt mixed forest

    # temperate mixed forest/ bogs - a woodland, TeBE >30%, TeBE dominate trees
    elif (lai_w >= WT * 0.5 and fr(tebe_trees + CldShrubPFT, woody=True) > 0.75): # and lai_s < 0.2:
        b = biome.MaF    # magellanic forest, woodland and bog

    # cheat 1: all other woodlands with TeNE go to mixed forest
    elif (lai_w >= WT) and (max_t('TeNE') or sum_lai(CldShrubPFT) > 0.3):
        b = biome.TeMF     # high-alt mixed forest

    # cheat 2: all other woodlands are MaF
    elif (lai_w >= WT): 
        b = biome.MaF    # magellanic forest, woodland and bog

    # patagonian steppe - lai_woody < woodland, lai hrubs < 0.1, cold shrubs > 10% woody, total lai > 0.2
    elif (lai_w < WT and fr(CldShrubPFT, woody=True) > 0.5) and lai_tot > 0.1:
        b = biome.PSt     # B=14

    # arid shrubland
    elif (lai_tot > 0.1): # and lai_w > 0.1):
        b = biome.AShr    # B=16 (b)

    # desert - total lai < 0.2
    elif (lai_tot <= 0.1):
        b = biome.D       # B=17
    else:
        print 'UNDEFINED'
        print lai_tot, lai_s, lai_w, lai_t, lai_g
        b = 99
        exit()
    return b



chile_es_biome  = enum(# patagonian steppe
                       'PSt',    # cold, grasses low boreal shrubs 000
                       # magellanic forest (cold, shrubs some trees)
                       'MaF',     #= 9,   # BIOME 8 - Grass Tundra 001
                       # temperate rainforest
                       'TeBEF',  # temperate evergreen 002
                       # temperate mixed forest
                       'TeMF',   # mixed, valdivian forest: evergreen + conifers 003
                       # decidious transition zone
                       'TeBS',   # Notofagus parklands 004
                       # sclerophyllous zone
                       'SclWF',  # Sclerophyllous Woodland/ forest 005
                       'Mat',    # Matorral 006
                       # desert zone
                       'AShr',   #= 10,  # BIOME 9 - Arid Shrubland 007
                       'D')      #= 10}) # BIOME 10 - Desert    008


class Biomization( object ):
    def __init__(self, biome_type='CHILE_ES_NEW'):
        self._pfts = None
        self._btype = biome_type

        if biome_type not in ['SMITH', 'CHILE', 'CHILE_ES', 'CHILE_ES_NEW']:
            exit()

        if self._btype == 'CHILE_ES_NEW':
            self._biomes = chile_es_biome

    def __call__(self, da):
        """Execute classification on supplied DataArray."""
        if self._pfts is None:
            self._pfts = da['pft'].values
        return self._classify(da)

    def _classify(self, da, class_var='sp_lai'):
        """Call classification."""
        return LUT_ES_CHILE_NEW(self._biomes, self._pfts, da)

    def biomes(self):
        """Return the biome names of this classification."""
        return self._biomes.content

    def pfts(self):
        """Return the pft names used in this biome classification."""
        return self._pfts



def main():
    """Main routine."""

    # INPUT
    #parser, args = cli()
    #cfg = get_config(cfgFile=args.config)

    #if args.storeconfig:
    #    set_config(cfg)

    fname = sys.argv[1]

    # load data
    with xr.open_dataset(fname) as ds:
        da_lai  = ds['sp_lai'].load()
        da_frac = ds['fraction'].load()
        lf_ids  = ds['lf_id'].values
        pfts    = ds['pft'].values


    # init the classifier
    biomizer = Biomization()

    print biomizer.biomes()
    print biomizer.pfts()

    # destination xr.DataArray
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


    da_biomes.attrs['_FillValue'] = NODATA
    da_biomes.attrs['units'] = 'biome_id'
    ds = da_biomes.to_dataset()
    ds['fraction'] = da_frac
    

    agg = []
    for b, bname in enumerate(biomizer.biomes()):
        fr=ds['fraction'].where(ds['biome_lf']==b)
        u=fr.sum(dim='lf_id')
        agg.append(u)

    da2 = xr.concat(agg, dim='biome').argmax(dim='biome').where(ds['fraction'].sum(dim='lf_id')>0)
    ds['biome'] = da2
    ds['biome'].attrs['_FillValue'] = NODATA
    ds['biome'].attrs['units'] = 'biome_id'

    ds.to_netcdf(fname[:-3] + '_biome.nc', format='NETCDF4_CLASSIC')

    # aggregate to elevation levels
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
            agg.append(u)

        e_da2 = xr.concat(agg, dim='biome').argmax(dim='biome').where(e_frac.sum(dim=['lf_id', 'lon'])>0)
        ele_L.append(e_da2)

    da_ele = xr.concat(ele_L, pd.Index( np.arange(100, 30*200+100, 200, dtype=np.int32), name='ele')).astype(np.int32)
    da_ele.name='biome'
    da_ele.attrs['units'] = 'biome_id'

    da_ele['ele'] = da_ele['ele'].astype(np.int32)
    da_ele['ele'].attrs['axis'] = 'X'
    da_ele['lat'].attrs['axis'] = 'Y'
    da_ele.T.to_netcdf(fname[:-3] + '_biome_latele.nc', format='NETCDF4_CLASSIC')

