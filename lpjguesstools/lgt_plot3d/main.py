"""FILE lgt_plot3d.main.py

Plot 3D tiles for given variable.

Christian Werner, SENCKENBERG Biodiversity and Climate Research Centre (BiK-F)
email: christian.werner@senkenberg.de
2017/10/05
"""

import datetime
import glob 
import logging
from matplotlib import pylab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mayavi import mlab
import numpy as np
import os
import pandas as pd
from PIL import Image
from pylab import cm
import string
import time
from tvtk.api import tvtk
from tvtk.common import configure_input_data, configure_source_data, \
                        is_old_pipeline
from tvtk.tools import visual

import xarray as xr


log = logging.getLogger(__name__)

# import constants
from . import NODATA
from . import ENCODING


# ------------------------

def redirect_vtk_messages():
    """ Can be used to redirect VTK related error messages to a
    file."""

    import vtk
    output=vtk.vtkFileOutputWindow()
    output.SetFileName("log.txt")
    #output.SendToStdErrOn()
    vtk.vtkOutputWindow().SetInstance(output)




def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def set_text(text_property):
    """
    Set the text to sane defaults
    
    Parameters
    ----------
    text_property: mayavi TextProperty object
    """
    text_property.font_size = 8
    text_property.bold = False
    text_property.italic = False
    text_property.font_family = 'arial'
    #text_property.color = (0,0,0)

def set_cbar_text(cbar,lut_manager='scalar'):
    if lut_manager == 'scalar':
        manager = cbar.parent.scalar_lut_manager
    elif lut_manager == 'vector':
        manager = cbar.parent.vector_lut_manager
    else:
        raise Exception()
    set_text(manager.title_text_property)
    set_text(manager.label_text_property)


def Arrow_From_A_to_B(fig, x1, y1, z1, x2, y2, z2):

    visual.set_viewer(fig)
    
    ar1=visual.arrow(x=x1, y=y1, z=z1, color=(1,0,0))
    ar1.length_cone=0.4

    arrow_length=np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    ar1.actor.scale=[arrow_length, arrow_length, arrow_length]
    ar1.pos = ar1.pos/arrow_length
    ar1.axis = [x2-x1, y2-y1, z2-z1]
    return ar1

# DO NOT USE
def image_from_array(ary):
    """ Create a VTK image object that references the data in ary.
        The array is either 2D or 3D with.  The last dimension
        is always the number of channels.  It is only tested
        with 3 (RGB) or 4 (RGBA) channel images.
        Note: This works no matter what the ary type is (accept
        probably complex...).  uint8 gives results that make since
        to me.  Int32 and Float types give colors that I am not
        so sure about.  Need to look into this...
    """

    sz = ary.shape
    dims = len(sz)
    # create the vtk image data
    img = tvtk.ImageData()

    if dims == 2:
        # 1D array of pixels.
        img.whole_extent = (0, sz[0]-1, 0, 0, 0, 0)
        img.dimensions = sz[0], 1, 1
        img.point_data.scalars = ary

    elif dims == 3:
        # 2D array of pixels.
        if is_old_pipeline():
            img.whole_extent = (0, sz[0]-1, 0, sz[1]-1, 0, 0)
        else:
            img.extent = (0, sz[0]-1, 0, sz[1]-1, 0, 0)
        img.dimensions = sz[0], sz[1], 1

        # create a 2d view of the array
        ary_2d = ary[:]
        ary_2d.shape = sz[0]*sz[1],sz[2]
        img.point_data.scalars = ary_2d

    else:
        raise ValueError("ary must be 3 dimensional.")

    return img




def main(cfg):
    
    # --screen currently has no effect as the render only works offline
    mlab.options.offscreen = True
    
    redirect_vtk_messages()
    #
    ds = xr.open_dataset(cfg.TILE)
    
    # TODO: potentially allow for more tiles to be merged
    # ds = xr.merge([ds1, ds2, ds3, ds4], compat='no_conflicts')

    
    flip = False
    if ds.lat[0] < ds.lat[-1]:
        flip = True
    
    # some default limits    
    varD = dict(slope=(0, 35),
            landform_class=(None, None), 
            landform=(1, 6),
            aspect=(None, None),
            elevation=(None, None))

    # limit var
    varD = { k: varD[k] for k in cfg.VARS }

    for VAR in sorted(cfg.VARS):
        
        VAR_MIN = None
        VAR_MAX = None
        if VAR in varD.keys():
            VAR_MIN, VAR_MAX = varD[VAR]

        if VAR_MIN is None: VAR_MIN = float(ds[VAR].min())
        if VAR_MAX is None: VAR_MAX = float(ds[VAR].max())

        z = ds['elevation'].to_masked_array() #.filled(0)
        z = np.ma.masked_where((z <= 0) | (z == np.nan), z)
        w = ds[VAR].to_masked_array() #.filled(0)

        if flip:
            log.warn('Flipping latitude')
            z = np.flipud(z)
            w = np.flipud(w)

        lat_min, lat_max, lon_min, lon_max = ds.lat.min(), ds.lat.max(), ds.lon.min(), ds.lon.max()

        WIDTH = len(ds.lat)



        fig = mlab.figure(size=(1200, 900), bgcolor=(0, 0, 0))

        # Create the data source
        src = mlab.pipeline.array2d_source(z, mask=np.ma.getmaskarray(z))

        # Add the additional scalar information 'w', this is where we need to be a bit careful,
        # see
        # http://code.enthought.com/projects/mayavi/docs/development/html/mayavi/auto/example_atomic_orbital.html
        # and
        # http://code.enthought.com/projects/mayavi/docs/development/html/mayavi/data.html
        dataset = src.mlab_source.dataset
        array_id = dataset.point_data.add_array(w.T.ravel())
        dataset.point_data.get_array(array_id).name = 'color'
        dataset.point_data.update()

        # Here, we build the very exact pipeline of surf, but add a
        # set_active_attribute filter to switch the color, this is code very
        # similar to the code introduced in:
        # http://code.enthought.com/projects/mayavi/docs/development/html/mayavi/mlab.html#assembling-pipelines-with-mlab

        # calculate warp_scale
        delta_ele = np.nanmax(z) - np.nanmin(z)
        bscale = 0.25
        WS = bscale - ((delta_ele / (delta_ele+float(WIDTH)))*bscale)

        warp = mlab.pipeline.warp_scalar(src, warp_scale=WS)
        normals = mlab.pipeline.poly_data_normals(warp)
        active_attr = mlab.pipeline.set_active_attribute(normals,
                                                point_scalars='color')
        surf = mlab.pipeline.surface(active_attr, colormap='gist_earth', vmin=VAR_MIN, vmax=VAR_MAX) 

        cbar = mlab.colorbar(title=VAR, orientation='horizontal')
        cbar.scalar_bar_representation.position = [0.1, 0.03]
        cbar.scalar_bar_representation.position2 = [0.8, 0.05]

        axes = mlab.axes(xlabel='Lat.', nb_labels=3, ylabel='Long.', zlabel='',
                     ranges = [lat_max, lat_min, lon_min, lon_max, np.nanmin(z), np.nanmax(z)]) 
        axes.axes.font_factor=0.8

        set_text(axes.axes.axis_label_text_property)
        set_text(axes.axes.axis_title_text_property)
        set_cbar_text(surf)

        label_fstring='%#2.1f'
        bar = surf.parent.scalar_lut_manager.scalar_bar_widget
        bar.scalar_bar_actor.label_format = label_fstring

        # north arrow
        # TODO: fix when run for 4 tiles
        north = mlab.text(WIDTH*0.5, 0,'N',z=1050, width=0.025, color=(1,0,0)) #, font_factor=0.5)
        a = Arrow_From_A_to_B(fig, WIDTH*0.5, 0, 1000, WIDTH*0.5 - 200, 0, 1000)


        # offline rendering
        rw = tvtk.RenderWindow(size=fig.scene._renwin.size, off_screen_rendering=1)
        rw.add_renderer(fig.scene._renderer)

        w2if = tvtk.WindowToImageFilter()
        w2if.magnification = fig.scene.magnification
        w2if.input = rw
        ex = tvtk.PNGWriter()
        ex.file_name = cfg.OUTFILE[:-4] + '_' + VAR + '.png'
        configure_input_data(ex, w2if.output)
        w2if.update()
        ex.write()

        mlab.close(all=True)


def main_old(cfg):
    """Main Script."""    


    ds = xr.open_dataset(cfg.TILE)

    # set mayavi and vtk controls
    #redirect_vtk_messages()
    mlab.options.offscreen = True

    # open tile 
    # load elevation
    ele_data = xr.open_dataset( cfg.TILE)['elevation'].to_masked_array() #.filled(0)
    ele_data = ele_data.filled(0)

    # load landform class array
    lf_data = xr.open_dataset( cfg.TILE )['landform_class'].to_masked_array() #.filled(0)

    asp_comp = lf_data%1

    asp_comp   = (lf_data // 10**0 % 10).filled(0)
    slope_comp = (lf_data // 10**1 % 10).filled(0)
    ele_comp   = (lf_data // 10**2 % 10).filled(0)

    print asp_comp.mean()

    a3d = []
    for i in [asp_comp, slope_comp, ele_comp]:
        print i.min(), i.max()
        x = ((i - i.min()) / (i.max() - i.min()))*255
        a3d.append(x)
        print x.mean()
    A3D = np.array(a3d, dtype=np.uint8)
    
    A3D = np.ascontiguousarray(np.array(a3d, dtype=np.uint8).transpose(1,2,0))
    
    #arr = np.random.uniform(size=(3,256,256))*255
    #arr = np.ascontiguousarray(arr.transpose(1,2,0))
    #img = Image.fromarray(arr, 'RGB')
    #img.save('out.png')

    site_no = 1

    #im = Image.fromarray(np.uint8(new_cmap(lf_data)*255))
    im = Image.fromarray(A3D, 'RGB')

    im2 = im.transpose(Image.ROTATE_90)
    im2.save("tmp/outfile_lf_%d.png" % site_no)
    img = tvtk.PNGReader()
    img.file_name="tmp/outfile_lf_%d.png" % site_no
    texture=tvtk.Texture(input_connection=img.output_port, interpolate=1)


    # create figure
    f = mlab.figure(size=(800, 600), bgcolor=(1, 1, 1))

    surfs = []
    
    # loop over months
    surfBase = mlab.surf(ele_data - 25,
                color=(0,0.1,0.6),
                opacity=.7,
                warp_scale=0.05, vmin=0, vmax=1600)

    surf = mlab.surf(ele_data,
                color=(1,1,1), 
                warp_scale=0.05, vmin=0, vmax=1600)
        
    surf.actor.enable_texture = True
    surf.actor.tcoord_generator_mode = 'plane'
    surf.actor.actor.texture = texture

    
    mlab.view(azimuth=330, elevation=50, distance=1785.)

    no_skip = False

    if cfg.SHOW:
        mlab.options.offscreen = True

        if no_skip:
            # create dummy legend
            import pylab 
            import matplotlib as mpl 

            # Make a figure and axes with dimensions as desired. 
            fig = pylab.figure(figsize=(.5,3))
            ax = fig.add_axes([0.05, 0.1, 0.1, 0.9])

            # Set the colormap and norm to correspond to the data for which 
            # the colorbar will be used. 
            cmap = mpl.cm.cool 
            norm = mpl.colors.Normalize(vmin=_min, vmax=_max)

            # ColorbarBase derives from ScalarMappable and puts a colorbar 
            # in a specified axes, so it has everything needed for a 
            # standalone colorbar.  There are many more kwargs, but the 
            # following gives a basic continuous colorbar with ticks 
            # and labels. 
            cb = mpl.colorbar.ColorbarBase(ax, cmap=new_cmap,
                                           norm=norm, 
                                           orientation='vertical')
            cb.set_label('fractional cover [%]')


            pylab.savefig('Bplots3d/focussites_3d_legend.png', bbox_inches='tight')

            del fig

            # black background legend
            fg_color = 'white'
            bg_color = 'black'
            fig = pylab.figure(figsize=(1,3),  facecolor=bg_color)
            ax = fig.add_axes([0.05, 0.1, 0.1, 0.9])

            # Set the colormap and norm to correspond to the data for which 
            # the colorbar will be used. 
            cmap = mpl.cm.cool 
            norm = mpl.colors.Normalize(vmin=_min, vmax=_max)

            # ColorbarBase derives from ScalarMappable and puts a colorbar 
            # in a specified axes, so it has everything needed for a 
            # standalone colorbar.  There are many more kwargs, but the 
            # following gives a basic continuous colorbar with ticks 
            # and labels. 
            cb = mpl.colorbar.ColorbarBase(ax, cmap=new_cmap,
                                           norm=norm, 
                                           orientation='vertical')
            cb.set_label('vegetation cover [%]', color=fg_color)

            # set colorbar tick color
            cb.ax.yaxis.set_tick_params(color=fg_color)
            # set colorbar edgecolor 
            cb.outline.set_edgecolor(fg_color)
            # set colorbar ticklabels
            plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

            pylab.savefig('Bplots3d/focussites_3d_legend_black.png', bbox_inches='tight', dpi=150,
                    facecolor=fig.get_facecolor(), edgecolor='none')

            del fig


        # write mayavi plot to file
        fig = mlab.gcf()
        
        log.info('Plotting')
        
        rw = tvtk.RenderWindow(size=fig.scene._renwin.size, off_screen_rendering=1)
        rw.add_renderer(fig.scene._renderer)

        w2if = tvtk.WindowToImageFilter()
        w2if.magnification = fig.scene.magnification
        w2if.input = rw
        ex = tvtk.PNGWriter()
        ex.file_name = cfg.OUTFILE
        configure_input_data(ex, w2if.output)
        w2if.update()
        ex.write()
        del fig, ex, rw, w2if

    else:
        mlab.show()


    # cut ocean


    log.info("Done")


