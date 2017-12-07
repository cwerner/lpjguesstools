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


from collections import OrderedDict
import numpy as np

# TODO move this into a tools file
def enum(*sequential, **named):
    # default enum unclassified (99)
    named.update(dict(unclassified=99))
    enums = OrderedDict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['content'] = [reverse[x] for x in enums.values()]
    enums['reverse_mapping'] = reverse
    return type('Enum', (), enums)

# D:     Desert
# Mat:   matorral/ arid shrubland (disabled)
# SclW:  Sclerophyllous woodland/ forest (incl. matorral)
# TDF:   Temperate broadleaved decidious forest
# MixF:  Mixed forest (leftovers?)
# NPL:   Nothofagus Parkland (paleo landscape not present today) ???
# VRF:   Valdivian Rainforest
# MFW:   Magellanic Forest/ Woodland
# CDF:   Cool decidious forest
# T:     Tundra
# St:    Steppe
biomes = enum('D', 'ASh', 'Mat', 'SclW', 'TDF', 'MixF', 'NPL', 'VRF', 'MeW', 'MFW', 'T', 'St')

# Mat='#f0d8b6', 

biome_color = {biomes.D    : '#eae4e6', # 0 
               biomes.ASh  : '#b889d4', # 1 
               #biomes.Mat  : '#e29688', # 2
               biomes.SclW : '#9a8e15', # 3 ok
               biomes.TDF  : '#87f8e3', # 4 nix
               biomes.MixF : '#gfe1ac', # 5 wenige
               biomes.NPL  : '#78201e', # 6 wenige (high-altitude, falsche klasse?, auch uebergang zu wueste!!!)
               biomes.VRF  : '#7dfc9c', # 7 ok
               biomes.MeW  : '#ccffff', # 8 nix
               biomes.MFW  : '#749079', # 9 ok (rel. viel!!!)
               biomes.T    : '#95a0b5', # 10 zu wenig (!)
               biomes.St   : '#7b566d', # 11 zu viel
               biomes.unclassified : '#ff62b9'}

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
    CldShrubPFT = ['BE_s', 'BS_s']          # cold zone shrubs
    ShrubPFT = SclShrubPFT + CldShrubPFT    # total shrubs

    # Special tree PFT groups
    mesic_tebe_trees  = ['TeBE_itm', 'TeBE_tm'] # scleophyllous trees
    xeric_tebe_trees  = ['TeBE_itscl', 'TeBE_tscl'] # scleophyllous trees
    tebe_trees = xeric_tebe_trees + mesic_tebe_trees
    
    tebe_woody = tebe_trees + SclShrubPFT
    
    boreal_trees = ['BBE_itm','BBS_itm']
    tebs_trees = ['TeBS_tm', 'TeBS_itm']    # temperate decidious (Notofagus)
        
    mesic_woody = mesic_tebe_trees + CldShrubPFT + ['TeNE'] + boreal_trees + tebs_trees
    xeric_woody = xeric_tebe_trees + SclShrubPFT
    
    
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
    b = 99
    
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
            elif d.fr(tebs_trees + ['TeNE']) >= 0.5:
                # Temperate Broadleaf Decidious Forest
                b = biome.TDF   # broadleaf decidious forest
            elif d.fr(boreal_trees) >= 0.5:
                # Magellanic Forest/ Woodland 
                b = biome.MFW
            else:
                # mixed (check)
                b = biome.MixF
        
        # N. parkland ???
        elif (d.lai_t >= 0.5 and (d.lai_g / d.lai_tot) > 0.66):
            b = biome.NPL
            
        # woodlands
        elif d.lai_w > WT:
            if d.fr(xeric_woody, woody=True) > 0.5:
                # Sclerophyllous Woodland/ Matorral
                b = biome.SclW
            elif d.fr(mesic_woody, woody=True) > 0.5:
                # Cold Broadleaf Decidious Forest
                b = biome.MeW # cold decidious forest/ woodland    
            else:
                # Tundra
                b = biome.T
                
        # matorral (low shrub)
        elif d.lai_tot > 0.2 and (d.lai_g / d.lai_tot) < 0.66:
            b = biome.ASh
        else:
            b = biome.St
            
            
    if b == 99:
        print 'UNCLASSIFIED'

    return b
