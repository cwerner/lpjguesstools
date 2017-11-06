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
    enums = OrderedDict(zip(sequential, range(len(sequential))), **named)
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['content'] = [reverse[x] for x in enums.values()]
    enums['reverse_mapping'] = reverse
    return type('Enum', (), enums)


biomes  = enum(# patagonian steppe
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

    #biome  = enum('TeBEF',  #= 1,   # BIOME 1 - Temperate Broadleaved Evergreen Forest
    #              'TeMF',   #= 2,   # BIOME 2 - Temperate Mixed Forest
    #              'SclWF',  #= 3,   # BIOME 3 - Sclerophyllous Woodland/ Forest
    #              'TeW',    #= 4    Temperate Woodland
    #              'ArT',    #= 7,   # BIOME 6 - Arctic/alpine Tundra
    #              'Gr',     #= 9,   # BIOME 8 - Grass Tundra
    #              'AShr',   #= 10,  # BIOME 9 - Arid Shrubland
    #              'D')      #= 10}) # BIOME 10 - Desert   

    # Categorize PFTs ------------------

    # Classify shrubs
    SclShrubPFT = ['TeE_s', 'TeR_s']        # sclerophyllous shrubs
    CldShrubPFT = ['BE_s', 'BS_s']          # cold zone shrubs
    ShrubPFT = SclShrubPFT + CldShrubPFT    # total shrubs

    # Special tree PFT groups
    tebe_trees = ['TeBE_tm', 'TeBE_itm']    # temperate broadleaf evergreen
    tebs_trees = ['TeBS_tm', 'TeBS_itm']    # temperate decidious (Notofagus)

    # Thresholds
    FT = 2.5   # forest/woodland threshold
    WT = 1.0  # woodland/grassland threshold

    # instance of data object
    d = DataLaiPFT( data, grass=['C3G'], shrubs=SclShrubPFT+CldShrubPFT )


    # Decision tree (new) -------------
    # temperate rainforest (TeBEF) - a forest, TeBE >30%, TeBE dominate trees
    if (d.lai_t >= FT and d.fr(tebe_trees) >= 0.3 and d.max_t(tebe_trees)):
        b = biome.TeBEF   # temperate rainforest

    # temperate mixed forest (TeMF) - a forest, TeBE <30%
    elif (d.lai_t >= WT and d.fr(tebe_trees) < 0.3 and d.fr(tebs_trees) <= 0.15):
        b = biome.TeMF   # temperate mixed forest

    # notofagus decidious woodland (TeBS) - a forest/ woodland, TeBS > 30%, TeBS dominante trees
    elif (d.lai_w >= WT and d.fr(tebe_trees) < 0.3 and d.fr(tebs_trees) > 0.15): # and max_t(tebs_trees)):
        b = biome.TeBS    # decidious forest (Notofagus)

    # high-altitude forests TeBS dominante trees
    elif (d.lai_w >= WT * 0.5 and d.fr(['TeBE_itscl']+SclShrubPFT, woody=True) < 0.5) and d.sum_lai(CldShrubPFT) < 0.3:   # and max_t(tebs_trees)):
        b = biome.TeMF    # high-alt mixed

    # sclerophyllous (SclWF) - woodland, sclerophyllous PFTs >= 50% of woody, cold shrubs < 0.05 LAI
    elif (d.lai_w >= WT * 0.5 and d.fr(['TeBE_itscl']+SclShrubPFT, woody=True) >= 0.5): # and sum_lai(CldShrubPFT) < 0.05:
        b = biome.SclWF   # sclerophyllous woodland

    # matorral (Mat) - less than a woodland, sclerophyllous PFTs >= 50%
    elif (d.lai_w < WT * 0.5 and d.fr(['TeBE_itscl']+SclShrubPFT, woody=True) >= 0.5 and d.lai_w > 0.1) :
        b = biome.Mat     # matorral

    # high-altitude forests TeBS dominante trees
    elif (d.lai_w < WT * 0.5 and d.fr(['TeBE_itscl']+SclShrubPFT, woody=True) < 0.5 and d.lai_w > 0.1) :
        b = biome.TeMF     # high-alt mixed forest

    # temperate mixed forest/ bogs - a woodland, TeBE >30%, TeBE dominate trees
    elif (d.lai_w >= WT * 0.5 and d.fr(tebe_trees + CldShrubPFT, woody=True) > 0.75): # and lai_s < 0.2:
        b = biome.MaF    # magellanic forest, woodland and bog

    # cheat 1: all other woodlands with TeNE go to mixed forest
    elif (d.lai_w >= WT) and (d.max_t('TeNE') or d.sum_lai(CldShrubPFT) > 0.3):
        b = biome.TeMF     # high-alt mixed forest

    # cheat 2: all other woodlands are MaF
    elif (d.lai_w >= WT): 
        b = biome.MaF    # magellanic forest, woodland and bog

    # patagonian steppe - lai_woody < woodland, lai hrubs < 0.1, cold shrubs > 10% woody, total lai > 0.2
    elif (d.lai_w < WT and d.fr(CldShrubPFT, woody=True) > 0.5) and d.lai_tot > 0.1:
        print "PATAGONIAN STEPPE!!!"
        b = biome.PSt     # B=14

    # arid shrubland
    elif (d.lai_tot > 0.1): # and lai_w > 0.1):
        b = biome.AShr    # B=16 (b)

    # desert - total lai < 0.2
    elif (d.lai_tot <= 0.1):
        b = biome.D       # B=17
    else:
        print 'UNDEFINED'
        b = 99
        exit()
    return b
