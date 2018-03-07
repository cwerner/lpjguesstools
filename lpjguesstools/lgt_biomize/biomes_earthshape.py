# -*- coding: utf-8 -*-
#
# biomes_earthshape.py
# ====================
#
# Christian Werner
# christian.werner@senckenberg.de

# NOTE:
# to be valid, this file requires:
# - the biomes enum
# - the classification function


import numpy as np

from ..tools import enum

# D:     Desert
# ASh:   Arid shrubland
# Mat:   matorral/ arid shrubland (disabled)
# SclW:  Sclerophyllous woodland/ forest (incl. matorral)
# DMF:   Decidious Maule forest
# MixF:  Mixed forest (leftovers?)
# NPL:   Nothofagus Parkland (paleo landscape not present today) ???
# VRF:   Valdivian Rainforest
# MeW:   Mesic woodland
# MFW:   Magellanic Forest/ Woodland
# CDF:   Cool decidious forest
# T:     Tundra
# St:    Steppe

biomes = enum('D', 'ASh', 'Mat', 'St', 'SclW', 'DMF', 'MixF', 'NPL', 'VRF', 'MeW', 'CDF', 'MFW') # 'T'

biome_color = {biomes.D    : '#8c8c8c', # 0 Desert
               biomes.ASh  : '#e2db84', # 1 Arid Shrubland
               biomes.Mat  : '#c4a727', # 2 Mattoral
               biomes.St   : '#bd8321', # 3 Steppe
               biomes.SclW : '#817821', # 4 Sclerophyllous Woodland
               biomes.DMF  : '#87f8e3', # 5 Decidious 'Maule' forest
               biomes.MixF : '#e57d7d', # 6 Mixed Forest
               biomes.NPL  : '#7a4d7e', # 7 Nothofagus Parkland
               biomes.VRF  : '#009b3d', # 8 Valdivian Rainforest
               biomes.MeW  : '#9bc372', # 9 Mesic Woodland
               biomes.CDF  : '#3f6cb0', # 10 Cold Decidious Forest               
               biomes.MFW  : '#0f3c17', # 11 Magellanic Forest/ Woodland
               #biomes.T    : '#0fa17e', # 11 zu wenig (!)
               biomes.unclassified : '#db087d'}

class DataLaiPFT( object ):
    def __init__(self, da_lai, grass=None, shrubs=None, trees=None):
        self._da = da_lai
        self._tree_pfts = []
        self._shrub_pfts = []
        self._grass_pfts = []
        self._all_pfts = self._da.pft.values
        self.lai_t = 0      # tree
        self.lai_w = 0      # woody
        self.lai_g = 0      # grass
        self.lai_tot = 0    # total

        if grass is not None:
            self.set_grass( grass )
        if shrubs is not None:
            self.set_shrubs( shrubs )
        if trees is not None:
            self.set_trees( trees )

        self.lai_tot = self.lai_w + self.lai_g


    def set_trees(self, x):
        # identify tree pfts (optional)
        if type(x) != list:
            x = [x]
        self._tree_pfts += x
        self.lai_t = sum([self._da.sel(pft=p).values for p in self._tree_pfts])
        self.lai_w = self.lai_t + self.lai_s

    def set_shrubs(self, x):
        # identify tree pfts
        if type(x) != list:
            x = [x]
        self._shrub_pfts += x
        self.lai_s = sum([self._da.sel(pft=p).values for p in self._shrub_pfts])

        # derive tree pfts from total - (shrub + grass)
        if len(self._grass_pfts) > 0 and len(self._tree_pfts) == 0:
            self._tree_pfts = [p for p in self._all_pfts if p not in self._grass_pfts + self._shrub_pfts]
            self.lai_t = sum([self._da.sel(pft=p).values for p in self._tree_pfts])
            self.lai_w = self.lai_t + self.lai_s

    def set_grass(self, x):
        # identify tree pfts
        if type(x) != list:
            x = [x]
        self._grass_pfts += x
        self.lai_g = sum([self._da.sel(pft=p).values for p in self._grass_pfts])

    def fr(self, pfts, woody=False):
        """Fraction of the sum of given PFTs with reference to woody or total LAI """

        # reference: woody or total
        if woody:
            lai_ref = self.lai_w
        else:
            lai_ref = self.lai_tot

        if type(pfts) == list:
            # list of values
            _pfts = [p for p in pfts if p in self._all_pfts]
            if lai_ref > 0:
                return sum([self._da.sel(pft=p).values for p in _pfts]) / float(lai_ref)
        else:
            # single value
            if pfts in all_pfts:
                return self._da.sel(pft=pfts).values / lai_ref
        return 0.0


    def max_t(self, pfts):
        """Check if any of the given PFTs dominate the tree PFTs """

        _TreePFT = [p for p in pfts if p in self._tree_pfts]
        if len(_TreePFT) > 0:
            _lai = np.array([self._da.sel(pft=p).values for p in _TreePFT])
            ix = np.argmax(_lai)
            if _TreePFT[ix] in pfts:
                return True
        return False


    def sum_lai(self, pfts):
        """Sum the lai for provided pfts."""
        return sum([self._da.sel(pft=p).values for p in pfts])



def classification(biome, all_pfts, data):
    """EarthShape Custom Biomization"""


    # Categorize PFTs ------------------

    # Classify shrubs
    SclShrubPFT = ['TeE_s', 'TeR_s']        # sclerophyllous shrubs
    CldShrubPFT = ['BE_s'] #, 'BS_s']          # cold zone shrubs
    ShrubPFT = SclShrubPFT + CldShrubPFT    # total shrubs

    # Special tree PFT groups
    mesic_tebe_trees  = ['TeBE_itm', 'TeBE_tm'] # scleophyllous trees
    xeric_tebe_trees  = ['TeBE_itscl'] # scleophyllous trees
    tebe_trees = xeric_tebe_trees + mesic_tebe_trees
    
    tebe_woody = tebe_trees + SclShrubPFT
    
    boreal_trees = ['BBE_itm','BBS_itm']
    tebs_trees = ['TeBS_tm', 'TeBS_itm']    # temperate decidious (Notofagus)
        
    mesic_woody = mesic_tebe_trees + CldShrubPFT + ['TeNE'] + boreal_trees + tebs_trees
    xeric_woody = xeric_tebe_trees + SclShrubPFT
    
    mesic_decidious = tebs_trees + ['TeNE'] + ['BBS_itm'] #+ ['BS_s']

    be_woody = ['TeBE_itm', 'TeBE_tm', 'TeBE_itscl', 'TeE_s', 'BBE_tm', 'BE_ts']
    bs_woody = ['TeBS_itm', 'TeBS_tm', 'BBS_itm', 'TeR_s'] # BS_S
    
    boreal_pfts = boreal_trees + CldShrubPFT
    
    mesic_trees = mesic_tebe_trees + tebs_trees + boreal_trees + ['TeNE']
        
    # Thresholds
    FT = 2.0  # forest/woodland threshold (lai_t)
    WT = 1.0  # woodland/grassland threshold (lai_w)

    # instance of data object
    d = DataLaiPFT( data, grass=['C3G'], shrubs=SclShrubPFT+CldShrubPFT )

    # Decision tree (new) -------------
    # biome_color = {biomes.D    : '#eae4e6', # 0 
    #                biomes.SclW : '#9a8e15', # 1 ok
    #                biomes.TDF  : '#87f8e3', # 2 nix
    #                biomes.MixF : '#gfe1ac', # 3 wenige
    #                biomes.NPL  : '#78201e', # 4 wenige (high-altitude, falsche klasse?, auch uebergang zu wueste!!!)
    #            biomes.VRF  : '#7dfc9c', # 5 ok
    #            biomes.CDF  : '#ccffff', # 6 nix
    #            biomes.MFW  : '#749079', # 7 ok (rel. viel!!!)
    #            biomes.T    : '#95a0b5', # 8 zu wenig (!)
    #            biomes.St   : '#7b566d', # 9 zu viel
    #            biomes.unclassified : '#ff62b9'}
    
    # default: missing
    b = 98
    
    # TODO FIX deserts also in LPJ-GUESS (!) 75mm threshold (Elbert Nature Geoscience 2012)

    if (d.lai_tot < 0.2):
        # desert
        b = biome.D
    else:
        # forest threshold
        if (d.lai_t >= FT ):
            # forests
            if d.fr(tebe_trees) >= 0.5:
                # Valdivian Rainforest
                b = biome.VRF
            elif d.fr(tebs_trees) > 0.5:
                # Temperate Broadleaf Decidious Forest
                b = biome.DMF   # decidious Maule forest
            elif d.fr(boreal_pfts, woody=True) >= 0.5 and d.fr(['BBE_itm']) > d.fr(['BBS_itm']): #d.fr(be_woody, woody=True) > d.fr(bs_woody, woody=True):
                # Magellanic Forest/ Woodland 
                b = biome.MFW
            elif d.fr(boreal_pfts, woody=True) >= 0.5 and d.fr(['BBE_itm']) <= d.fr(['BBS_itm']): #d.fr(be_woody, woody=True) <= d.fr(bs_woody, woody=True):
                # Cold decidious Forest/ Woodland 
                b = biome.CDF
            else:
                # mixed (check)
                b = biome.MixF

        # N. parkland ???
        elif d.lai_t >= 1.0 and (d.lai_g / d.lai_tot) >= 0.5:
            b = biome.NPL

        # intermediate classes
        elif (d.lai_t >= 1.0 and d.fr(boreal_pfts, woody=True) >= 0.5 and d.fr(['BBE_itm']) > d.fr(['BBS_itm'])):  #d.fr(be_woody, woody=True) > d.fr(bs_woody, woody=True)):
            # Magellanic Forest/ Woodland (2nd chance)
            b = biome.MFW
        elif (d.lai_t >= 1.0 and d.fr(boreal_pfts, woody=True) >= 0.5 and d.fr(['BBE_itm']) <= d.fr(['BBS_itm'])): #d.fr(be_woody, woody=True) <= d.fr(bs_woody, woody=True)):
            # Cold Decidious Forest/ Woodland (2nd chance)
            b = biome.CDF
            
        # woodlands
        elif d.lai_w > WT:
            if d.fr(xeric_woody, woody=True) >= 0.5:
                # Sclerophyllous Woodland/ Matorral
                b = biome.SclW
            elif d.fr(mesic_woody, woody=True) > 0.5 and d.fr(boreal_pfts, woody=True) >= 0.5 and d.fr(['BBE_itm']) <= d.fr(['BBS_itm']):
                b = biome.CDF # mesic woodland
            elif d.fr(mesic_woody, woody=True) > 0.5 and d.fr(boreal_pfts, woody=True) >= 0.5 and d.fr(['BBE_itm']) > d.fr(['BBS_itm']):
                b = biome.MFW # magellanic forest woodland
            else:
                b = biome.MeW
                
        # matorral (low shrub)
        elif d.lai_tot >= 0.2 and (d.lai_g / d.lai_tot) < 0.7 and d.fr(xeric_woody) >= d.fr(mesic_woody):
            # arid shrubland or mattoral
            if d.lai_tot >= 0.5:
                b = biome.Mat
            else:
                b = biome.ASh            
                
        elif d.lai_tot >= 0.2  and (d.lai_g / d.lai_tot) < 0.7 and d.fr(xeric_woody) < d.fr(mesic_woody):
            # mesic woodland
            #if d.sum_lai(['BBE_itm']) >= 0.1:
            #    b = biome.MFW
            #else:
            b = biome.MeW
    
        elif (d.lai_g / d.lai_tot) >= 0.7:
            b = biome.St
        else:
            b = 100
            
    if b >= 90:
        print('UNCLASSIFIED', b)
        b = 99
    return b
