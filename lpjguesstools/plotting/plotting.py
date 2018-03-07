from collections import Sequence

import numpy as np
import xarray as xr

import math
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

from ..tools import enum


__all__ = ['ElevationTransect', 'Map', 'MapContainer']


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

# reduce clutter in contour labels
def clean_contours(contours, style='%.2f', zorder=None):
    """Clean the contour labels (max 2 per item)"""
    def path_length(path):
        v = path.vertices
        dv = np.diff(v, axis=0)
        return np.sum(np.sqrt(np.sum(dv**2, axis=-1)))

    # remain the longest path and remove all others
    deleted_path = []
    for c in contours.collections:
        paths = c.get_paths()
        if len(paths) > 2:
            paths.sort(key=path_length, reverse=True)
            for p in paths[2:]:
                deleted_path.append((c, p))
            del paths[2:]

    # create labels
    if zorder:
        r = plt.clabel(contours, contours.levels, fmt=style, 
                        inline=True, fontsize=8, zorder=zorder)
    else:
        r = plt.clabel(contours, contours.levels, fmt=style,
                        inline=True, fontsize=8)
        
    # restore all removed paths
    for c, p in deleted_path:
        c.get_paths().append(p)


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
            print('Specified vmin > vmax (using vmin, vmax from data)')
            kwargs['vmin'] = VMIN
            kwargs['vmax'] = VMAX

    return kwargs



def draw_map_layout(ax, orientation=None, country=None, name=None, 
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


def draw_elevationtransect_layout(ax, orientation=None, name=None, 
                                  left_label=True, bottom_label=True, **kwargs):
    """Draw map decorations
    """
    
    # gridlines
    gl = ax.grid(which='major', linestyle='--', linewidth=1, zorder=9999, color='gray')

    ax.set_yticks([-60, -50, -40, -30, -20, -10])
    ax.set_xticks([0, 2000, 4000, 6000])

    def y_fmt(x, y):
        hemisphere = 'N'
        if x <= 0:
            hemisphere = 'S'
        return '%d$\degree$%s' % (abs(x), hemisphere)

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(y_fmt))
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    
    #gl.ylocator = mticker.FixedLocator([-60, -50, -40, -30, -20, -10])
    #gl.xlabels_top = False
    #gl.ylabels_right = False
    #gl.xlabels_bottom = False
    #gl.ylabels_left = False
    
    #print ax, left_label, bottom_label
    if left_label == False:
        ax.tick_params(labelleft='off', left='off')
        #ax.set_yticks([])
    if bottom_label == False:
        ax.tick_params(labelbottom='off', bottom='off')  
    #    ax.set_xticks([])

    #gl.yformatter = LATITUDE_FORMATTER    
    #gl.xlabel_style = {'size': 10, 'color': fg_color}
    #gl.ylabel_style = {'size': 10, 'color': fg_color}

    # panel label
    do_subscript = False
    
    if name:
        # hack to create subscript
        if do_subscript:
            x = name
            if '_itm' in x:
                x = x.replace('_itm', '$_{itm}$')
            elif '_tm' in x:
                x = x.replace('_tm', '$_{tm}$')
            elif 'TeBE_itscl' in x:
                x = x.replace('_itscl', '$_{itscl}$')
            elif 'TeB_s' in x:
                x = 'MB$_s$'
            elif 'TeE_s' in x:
                x = 'ME$_s$'
            elif 'TeE_s' in x:
                x = 'B$_s$'
            name = x
            
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
    default_cbar_size = 10

    def __init__(self, names=None, orientation=None, figsize=None, 
                 nrows_ncols=None, spatial=True, cbar_size=None, **kwargs):
        self.orientation = orientation if orientation else self.default_orientation
        self.names = names if names else self.default_names
        self.figsize = figsize if figsize else self.default_figsize        
        self.cbar_size = cbar_size if cbar_size else self.default_cbar_size        

        if self.orientation == 'horizontal':
            simple_nrows_ncols = (1, len(names))
        else:
            simple_nrows_ncols = (len(names), 1)
        self.nrows_ncols = nrows_ncols if nrows_ncols else simple_nrows_ncols

        self.fig = plt.figure(figsize=self.figsize)

        axisgrid_kwargs = dict(nrows_ncols=self.nrows_ncols,
                               axes_pad = 0.2,
                               ngrids=len(self.names),
                               label_mode = "",
                               cbar_location = "right",
                               cbar_mode="single",
                               cbar_size=str(self.cbar_size)+'%',
                               share_all=False)
        if spatial:
            projection = ccrs.PlateCarree()
            axes_class = (GeoAxes, dict(map_projection=projection))
            axisgrid_kwargs['axes_class'] = axes_class
        else:
            axisgrid_kwargs['aspect'] = False
            axisgrid_kwargs['axes_pad'] = 0.3


        self.axes = AxesGrid(self.fig, 111, **axisgrid_kwargs)

        self._plots = []
        self.map = []

        for i in range(len(self.names)):
            if spatial:
                m_ = Map(ax=self.axes[i], name=self.names[i], **kwargs)
            else:
                m_ = ElevationTransect(ax=self.axes[i], name=self.names[i], **kwargs)
           
            self.map.append( m_ )                

        # handle labels in multi-plots
        bottom_ids = []
        left_ids = []

        n = self.nrows_ncols[1]
        ids = range(len(self.names))
        ids_2d = [ids[i:i+n] for i in xrange(0, len(ids), n)]

        def label_in_second_last_row(ids_2d):
            """Identify second-last row ids that need bottom labels"""
            ids = []
            if len(ids_2d) > 1:
                second_last_row = ids_2d[-2]
                last_row = ids_2d[-1]
                ids = second_last_row[len(last_row):]
            return ids
        
        if self.nrows_ncols[0] == 1:
            # one row only, all xlabels
            left_ids = [0]
            bottom_ids = range(len(self.axes))
        elif self.nrows_ncols[1] == 1:
            # one column only, all ylabels
            left_ids = range(len(self.axes))
            bottom_ids = range(len(self.axes))[-1:]
        else:
            # identify left and bottom axis in multi row, col
            left_ids = [x[0] for x in ids_2d]
            bottom_ids = ids_2d[-1] + label_in_second_last_row(ids_2d)
        

        for i, ax in enumerate(self.axes):
            left_label=True
            bottom_label=True
            if i not in left_ids:
                left_label=False
            if i not in bottom_ids:
                bottom_label=False
            
            self.map[i].drawlayout(left_label=left_label, bottom_label=bottom_label, **kwargs)
        
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

            if np.isnan(data.values).all():
                self._plots.append(None)
            else:                        
                p = data.plot(ax=self.map[i].ax, zorder=1000, 
                            add_colorbar=True, cbar_ax=self.axes.cbar_axes[0], **kwargs)
                self._plots.append(p)

    def plot_data_discrete(self, variable='', **kwargs):
        # joined limits
        
        levels = []
        colors = []
    
        # required arguments: levels, colors 
        if 'biome' in kwargs.keys():
            print(type(kwargs['biome']))
            levels = kwargs['biome'].items
            del kwargs['biome']
             # get keys from enum
            
        if 'color_mapping' in kwargs.keys():
            for l in levels:
                print(l)
                print(kwargs['color_mapping'])
                colors.append(kwargs['color_mapping'][l])
            del kwargs['color_mapping']

        kwargs = check_data_limits(self.map, variable, kwargs)
        # get joined data limits
        
        for i in range(len(self.map)):
            if variable:
                data = self.map[i].data[variable]
            else:                 
                data = self.map[i].data

            if np.isnan(data.values).all():
                self._plots.append(None)
            else:                        
                p = data.plot(ax=self.map[i].ax, zorder=1000, 
                            levels=levels, colors=colors,
                            add_colorbar=False, **kwargs)
                self.map[i].ax.set_ylabel("")
                disable_xlabel = True
                if len(self.map)<=2:
                    disable_xlabel = False
                elif i==int(math.ceil(len(self.map)/2)) and len(self.map)%2 != 0:
                    disable_xlabel = False
                else:
                    pass

                if disable_xlabel:
                    self.map[i].ax.set_xlabel("")
                else:
                    self.map[i].ax.set_xlabel("Elevation [m]")


                self.map[i].ax.set_ylabel("")
                self._plots.append(p)

        # create discrete legend
        # self.axes.cbar_axes[0]
        self.fig.delaxes(self.axes.cbar_axes[0]) 
        
    
    def plot_data_contour(self, variable='', clevels=None, **kwargs):
        # joined limits

        kwargs = check_data_limits(self.map, variable, kwargs)
        # get joined data limits
        
        # use second set of levels for contour lines, otherwise use same
        if 'levels' in kwargs.keys():
            clevels = clevels if clevels else kwargs['levels']
        
        for i in range(len(self.map)):
            if variable:
                data = self.map[i].data[variable]
            else:                 
                data = self.map[i].data
            
            # upsample data for contour by 5
            import scipy.ndimage
            n = 3   # multiplier (should be odd)
            lats = data.coords['lat'].values.tolist()
            eles = data.coords['ele'].values.tolist()
            step_lat = lats[1] - lats[0]
            step_ele = eles[1] - eles[0]
            
            # this is 0.01 for now (to avail covering everything in zeros)
            data_min = max(data.to_masked_array().min(), 0.01)
            # this needs adjusting for larger multipliers
            
            def calc_mod(step, n):
                """Calulate the step modifications"""
                x = (n-1)/2
                return (np.arange(-x,x+1)*(step/float(n))).tolist()
            
            mod_lat = calc_mod(step_lat, n)
            mod_ele = calc_mod(step_ele, n)
            mask = np.ma.getmaskarray(data.to_masked_array()).repeat(n, 1).repeat(n, 0)
            
            array_z = scipy.ndimage.zoom(data.to_masked_array().filled(0.0), n)
            lat_z = np.repeat(lats, n) + np.array(mod_lat * len(lats))
            ele_z = np.repeat(eles, n) + np.array(mod_ele * len(eles))
            array_z = np.where(array_z < data_min, data_min, array_z)
            array_zm = np.ma.array(array_z, mask=mask)
            data_z = xr.DataArray(array_zm, coords=[('lat', lat_z), ('ele', ele_z)], dims=['lat','ele'])
            
            try:
                p = data_z.plot.contourf(ax=self.map[i].ax, zorder=1000, 
                        add_colorbar=True, cbar_ax=self.axes.cbar_axes[0], **kwargs)

                # contour lines                
                cs = data_z.plot.contour(ax=self.map[i].ax, zorder=1000,
                            levels=clevels, colors=('k',),
                            linewidths=(0.5,), add_colorbar=False,
                            linestyles=('dashed',))
                self._plots.append(p)
            except:
                self._plots.append(None)


            # clean labels
            #
            # clabel: %d labels as decimals, %.1f labels as 1-digit float
            #         default is %.2f
            #if 'clabel' in kwargs.keys():
            #    clean_contours(cs, style=clabel, zorder=1001)
            #else:
            #    clean_contours(cs, zorder=1001)
            #
            self.map[i].set_xlabel('')
            self.map[i].set_ylabel('')
            self.map[i].set_xlim(left=0)
            



class Map( object ):
    """Map class
    """
    
    default_projection = ccrs.PlateCarree()
    default_extent = [-76, -64, -57, -16]

    def __init__(self, ax=None, projection=None, extent=None, data=None, name=None, drawlayout=False, **kwargs):
        self.projection = projection if projection else self.default_projection
        self.extent = extent if extent else self.default_extent
        self.data = data
        self.name = name

        self.ax = ax if ax else plt.axes(projection=self.projection, **kwargs)
        self.ax.set_extent(self.extent)

        if drawlayout:
            self.drawlayout(**kwargs)

    def __getattr__(self, item):
        return getattr(self.ax, item)

    def __dir__(self):
        self.ax.__dir__() + ['drawlayout']

    def drawlayout(self, country=None, locations=None, left_label=True, bottom_label=True, **kwargs):
        draw_map_layout(self.ax, name=self.name, country=country, 
                        left_label=left_label, bottom_label=bottom_label, **kwargs)

        if locations:
            draw_locations(self.ax, locations, **kwargs)
    
    def drawlocations(self, locations, **kwargs):
            draw_locations(self.ax, locations, **kwargs)


class ElevationTransect( object ):
    """ElevationTransect class
    """
    
    def __init__(self, ax=None, extent=None, data=None, name=None, drawlayout=False, **kwargs):
        self.data = data
        self.name = name

        self.ax = ax if ax else plt.axes(**kwargs)

        if drawlayout:
            self.drawlayout(**kwargs)

    def __getattr__(self, item):
        return getattr(self.ax, item)

    def __dir__(self):
        self.ax.__dir__() + ['drawlayout']

    def drawlayout(self, left_label=True, bottom_label=True, **kwargs):
        draw_elevationtransect_layout(self.ax, name=self.name, left_label=left_label,
                                      bottom_label=bottom_label, **kwargs)
