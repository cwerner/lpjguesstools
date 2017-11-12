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


def drawmap(ax, orientation=None, country=None, name=None, **kwargs):
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

    if orientation:
        if orientation == 'vertical':
            gl.xlabels_bottom = False
        elif orientation == 'horizontal':
            gl.ylabels_left = False


    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER    
    gl.xlabel_style = {'size': 12, 'color': fg_color}
    gl.ylabel_style = {'size': 12, 'color': fg_color}

    # panel label
    if name:
        ax.annotate(name, xy=(0, 1), xycoords='axes fraction', fontsize=12,
            xytext=(5, -5), textcoords='offset points', ha='left', va='top', zorder=10000,
            bbox={'facecolor':'white', 'alpha':1, 'pad':5}, clip_on=True)


def draw_locations(ax, locations, label_locations=False):
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

    def __init__(self, names=None, orientation=None, figsize=None, **kwargs):
        self.orientation = orientation if orientation else self.default_orientation
        self.names = names if names else self.default_names
        self.figsize = figsize if figsize else self.default_figsize
        projection = ccrs.PlateCarree()
        axes_class = (GeoAxes, dict(map_projection=projection))

        self.fig = plt.figure(figsize=self.figsize)
        self.axes = AxesGrid(self.fig, 111, axes_class=axes_class,
                             nrows_ncols=(1, len(self.names)),
                             axes_pad = 0.2,
                             share_all=True,
                             label_mode = "",
                             cbar_location = "right",
                             cbar_mode="single",
                             cbar_size='10%')


        self._plots = []
        self.map = []

        for i, ax in enumerate(self.axes):
            m_ = Map(ax=ax, name=self.names[i], **kwargs)
            self.map.append( m_ )

        # handle labels in multi-plots
        if self.orientation == 'horizontal':
            self.map[0].drawmap(**kwargs)
            if len(self.map) > 1:
                for i in range(1, len(self.map)):
                    self.map[i].drawmap(orientation=self.orientation, **kwargs)
        elif self.orientation == 'vertical':
            if len(self.map) > 1:
                for i in range(len(self.map)-1):
                    self.map[i].drawmap(orientation=self.orientation, **kwargs)
            self.map[-1].drawmap(**kwargs)
        else:
             for i in range(len(self.map)):
                self.map[i].drawmap(**kwargs)
        
        # init abc
        #super().__init__()

    def __getitem__(self, i):
        return self.map[i]
    
    def __len__(self):
        return len(self.map)

    def add(self, ix, data):
        self.map[ix].data = data


    def _get_joined_limits(self, var):
        min_val = min([x.data[var].min() for x in self.map])
        max_val = max([x.data[var].max() for x in self.map])
        return (min_val, max_val)

    def to_file(self, filename, dpi=None):
        if dpi:
            self.fig.savefig(filename, bbox_inches='tight', dpi=dpi)
        else:
            self.fig.savefig(filename, bbox_inches='tight')


    def plot_data(self, variable='', **kwargs):
        # joined limits

        VMIN, VMAX = self._get_joined_limits(variable)
        for i in range(len(self.map)):
            p = self.map[i].data[variable].plot(ax=self.map[i].ax, zorder=1000, vmin=VMIN, vmax=VMAX, 
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

    def drawmap(self, orientation=None, country=None, locations=None, **kwargs):
        drawmap(self.ax, name=self.name, orientation=orientation, country=country, **kwargs)

        if locations:
            draw_locations(self.ax, locations, **kwargs)



"""
if __name__ == '__main__':

    LVL = [-20, -10,-5,-1.5,1.5, 5, 10, 15, 20]
    data1 = xr.open_dataset('../../2017-11-07/data.nstd4_000_avg.nc').mean(dim='month', skipna=False)
    data2 = xr.open_dataset('../../2017-11-07/data.nstd4_155_avg.nc').mean(dim='month', skipna=False)
    data3 = xr.open_dataset('../../2017-11-07/data.nstd4_220_avg.nc').mean(dim='month', skipna=False)

    data_list = [data1, data2, data3]

    # create container
    
    mc = MapContainer(['LGM', 'MH', 'PD'], orientation='horizontal')

    # add data
    for i, data in enumerate(data_list):
        mc.add(i, data)

    # plot data
    mc.plot_data(variable='sp_mtemp', levels=LVL, cmap='RdBu_r')

    # save
    plt.savefig('dummy2.png', bbox_inches='tight', dpi=150) 
"""

