from collections import Sequence

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.legend_handler import HandlerPatch
from mpl_toolkits.axes_grid1 import AxesGrid

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.geoaxes import GeoAxes



__all__ = ['Map', 'MapContainer']


fg_color = 'black'


countries = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='50m',
    facecolor='none')

ocean = cfeature.NaturalEarthFeature(
        category='physical', 
        name='ocean', 
        scale='50m',
        edgecolor='face',
        facecolor=cfeature.COLORS['water'])


def get_country_outline(country_name):
    """Return the country border polygone
    """

    shpfilename = shpreader.natural_earth(resolution='50m',
                                      category='cultural',
                                      name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries = reader.records()

    country_outline = [c.geometry for c in countries if c.attributes['adm0_a3'] == country_name]
    
    if len(country_outline) > 1:
        return country_outline[0]
    return country_outline


def check_data_limits(mc, var, kwargs):
    """Determine joined vmin, vmax values for MapContainer
    """
    if var:
        VMIN = min([x.data[var].min() for x in mc])
        VMAX = max([x.data[var].max() for x in mc])
    else:
        VMIN = min([x.data.min() for x in mc])
        VMAX = max([x.data.max() for x in mc])
        
    if 'center' in kwargs.keys():
        if abs(VMIN) < abs(VMAX):
            if 'vmin' not in kwargs.keys():
                kwargs['vmin'] = -VMAX             
        else:
            if 'vmin' not in kwargs.keys():
                kwargs['vmin'] = VMIN
        if 'vmax' in kwargs.keys():
            kwargs['vmin'] = -kwargs['vmax']
            del kwargs['vmax'] 
    else:
        if 'vmin' not in kwargs.keys():
            kwargs['vmin'] = VMIN
        if 'vmax' not in kwargs.keys():
            kwargs['vmax'] = VMAX

    if 'center' not in kwargs.keys():
        # check valid vmin/ vmax
        if kwargs['vmin'] > kwargs['vmax']:
            print 'Specified vmin > vmax (using vmin, vmax from data)'
            kwargs['vmin'] = VMIN
            kwargs['vmax'] = VMAX

    return kwargs



def drawmap(ax, orientation=None, country=None, name=None, 
            left_label=True, bottom_label=True, **kwargs):
    """Draw map decorations
    """
    
    country_fill = []
    country_name = []
    
    # only single entries for the moment
    if country:
        if type(country) == dict:
            for k, v in country.iteritems():
                country_name.append(k)
                if 'fill' in v.keys():
                    country_fill.append( v['fill'] )
                else:
                    country_fill.append( None )
        elif type(country) == list:
            country_name = country
            country_fill = [None] * len(country)
        else:
            country_name = [country]
            country_fill = [None]
        country_outline = [get_country_outline(c) for c in country_name]

    # backdrop (ocean, countries)
    ax.add_feature(ocean)                                    # ocean
    ax.add_feature(countries, edgecolor='gray', zorder=9998) # country borders
    
    # individual country (if given [background with fill])
    if country:
        for c, c_fill in zip(country_outline, country_fill):
            if c_fill:
                ax.add_geometries(c,  ccrs.PlateCarree(), 
                    edgecolor='black', facecolor=c_fill)
                  
    ax.coastlines(resolution='50m')                 # coastlines

    # overlay country outline (ontop, black outline)
    if country:
        for c in country_outline:
            ax.add_geometries(c,  ccrs.PlateCarree(), 
                  edgecolor='black', facecolor='none', zorder=9999)

    # gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                  linestyle='--', linewidth=1, zorder=9999, color='gray')
    gl.xlocator = mticker.FixedLocator([-80, -75, -70, -65, -60])
    gl.ylocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabels_bottom = False
    gl.ylabels_left = False
    if left_label:
        gl.ylabels_left = True
    if bottom_label:
        gl.xlabels_bottom = True

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER    
    gl.xlabel_style = {'size': 10, 'color': fg_color}
    gl.ylabel_style = {'size': 10, 'color': fg_color}

    # panel label
    if name:
        ax.annotate(name, xy=(0, 1), xycoords='axes fraction', fontsize=12,
            xytext=(5, -5), textcoords='offset points', ha='left', va='top', zorder=10000,
            bbox={'facecolor':'white', 'alpha':1, 'pad':5}, clip_on=True)


def draw_locations(ax, locations, label_locations=False, **kwargs):
    for location, args in locations.iteritems():

        ax.scatter(args['lon'], args['lat'], facecolor='black', edgecolor='black', 
                   s=25, lw=.5, zorder=10000)

        if label_locations:
            ax.annotate(location, (args['lon'] + .75, args['lat']),
                        horizontalalignment=args['just'],
                        verticalalignment=args['vert'], 
                        bbox={'boxstyle': 'square', 'fc':'w', 'ec':'k'}, 
                        zorder=10001,
                        fontsize=8,
                        clip_on=True)


class MapContainer( Sequence ):
    """Map Container class
    """

    default_figsize = (7,7)
    default_orientation = 'horizontal'
    default_legend_position = 'right'
    default_names = ['untitled']

    def __init__(self, names=None, orientation=None, figsize=None, nrows_ncols=None, **kwargs):
        self.orientation = orientation if orientation else self.default_orientation
        self.names = names if names else self.default_names
        self.figsize = figsize if figsize else self.default_figsize
        projection = ccrs.PlateCarree()
        axes_class = (GeoAxes, dict(map_projection=projection))

        if self.orientation == 'horizontal':
            simple_nrows_ncols = (1, len(names))
        else:
            simple_nrows_ncols = (len(names), 1)
        self.nrows_ncols = nrows_ncols if nrows_ncols else simple_nrows_ncols

        self.fig = plt.figure(figsize=self.figsize)
        self.axes = AxesGrid(self.fig, 111, axes_class=axes_class,
                             nrows_ncols=self.nrows_ncols,
                             axes_pad = 0.2,
                             ngrids=len(self.names),
                             #share_all=True,
                             label_mode = "",
                             cbar_location = "right",
                             cbar_mode="single",
                             cbar_size='10%')


        self._plots = []
        self.map = []

        for i in range(len(self.names)):
            m_ = Map(ax=self.axes[i], name=self.names[i], **kwargs)
            self.map.append( m_ )

        # handle labels in multi-plots
        bottom_ids = []
        left_ids = []
        
        if self.nrows_ncols[0] == 1:
            # one row only, all xlabels
            bottom_ids = range(len(self.axes))
        elif self.nrows_ncols[1] == 1:
            # one column only, all ylabels
            left_ids = range(len(self.axes))
        else:
            # identify left and bottom axis in multi row, col
            n = self.nrows_ncols[1]
            ids = range(len(self.names))
            ids_2d = [ids[i:i+n] for i in xrange(0, len(ids), n)]
            left_ids = [x[0] for x in ids_2d]
            bottom_ids = ids_2d[-1]
        
        for i, ax in enumerate(self.axes):
            left_label=True
            bottom_label=True
            if i not in left_ids:
                left_label=False
            if i not in bottom_ids:
                bottom_label=False
            self.map[i].drawmap(left_label=left_label, bottom_label=bottom_label, **kwargs)
        
    def __getitem__(self, i):
        return self.map[i]
    
    def __len__(self):
        return len(self.map)

    def add(self, ix, data):
        self.map[ix].data = data

    def to_file(self, filename, dpi=None):
        if dpi:
            self.fig.savefig(filename, bbox_inches='tight', dpi=dpi)
        else:
            self.fig.savefig(filename, bbox_inches='tight')


    def plot_data(self, variable='', **kwargs):
        # joined limits

        kwargs = check_data_limits(self.map, variable, kwargs)
        # get joined data limits
        
        for i in range(len(self.map)):
            if variable:
                data = self.map[i].data[variable]
            else:                 
                data = self.map[i].data            
            p = data.plot(ax=self.map[i].ax, zorder=1000, 
                        add_colorbar=True, cbar_ax=self.axes.cbar_axes[0], **kwargs)
            self._plots.append(p)


class Map( object ):
    """Map class
    """
    
    default_projection = ccrs.PlateCarree()
    default_extent = [-76, -64, -57, -16]

    def __init__(self, ax=None, projection=None, extent=None, data=None, name=None, drawmap=False, **kwargs):
        self.projection = projection if projection else self.default_projection
        self.extent = extent if extent else self.default_extent
        self.data = data
        self.name = name

        self.ax = ax if ax else plt.axes(projection=self.projection, **kwargs)
        self.ax.set_extent(self.extent)

        if drawmap:
            self.drawmap(**kwargs)

    def __getattr__(self, item):
        return getattr(self.ax, item)

    def __dir__(self):
        self.ax.__dir__() + ['drawmap']

    def drawmap(self, country=None, locations=None, left_label=True, bottom_label=True, **kwargs):
        drawmap(self.ax, name=self.name, country=country, 
                left_label=left_label, bottom_label=bottom_label, **kwargs)

        if locations:
            draw_locations(self.ax, locations, **kwargs)
    
    def drawlocations(self, locations, **kwargs):
            draw_locations(self.ax, locations, **kwargs)
        