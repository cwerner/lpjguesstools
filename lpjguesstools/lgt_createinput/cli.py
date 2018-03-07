# -*- coding: utf-8 -*-
"""lgt_createinput.cli: Commandline interpreter for lgt_createinput.py."""

import click
import logging

from .main import main

# import constants
from .. import EPILOG


log = logging.getLogger(__name__)

class Bunch(object):
    """Simple data storage class."""
    def __init__(self, adict):
        self.__dict__.update(adict)
    def overwrite(self, adict):
        self.__dict__.update(adict)

# command line arguments
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
click.Context.get_usage = click.Context.get_help

@click.command(context_settings=CONTEXT_SETTINGS, epilog=EPILOG)

@click.option('--classfication', type=click.Choice(['SIMPLE', 'WEISS']), 
                    default='SIMPLE', show_default=True,
                    help='classification scheme')

@click.option('--cutoff', default=1.0, show_default=True,
                    help='required area fraction [%]')

@click.option('--dems', metavar='PATH',
                    help='source for DEM files')

@click.option('--extent', nargs=4, type=click.FLOAT, metavar='LON1 LAT1 LON2 LAT2',
                    help='extent of output netcdf files')

@click.option('--force-overwrite', is_flag=True, default=False, 
                    help='overwrite tiles even if they already exists')

@click.option('--gridlist', default='gridlist.txt', 
                    help='name of created gridlist file')

@click.option('--masks', metavar='PATH',
                    help='source for water masks (shp)')

@click.option('--verbose', is_flag=True, 
                    help='increase logging info')

@click.version_option()

@click.argument('storage', type=click.Path(exists=True))
@click.argument('outdir', type=click.Path(exists=True)) 

def cli(cutoff, dems, masks, gridlist, extent, classfication, storage, outdir, verbose, force_overwrite):
    """LPJ-GUESS 4.0 subpixel mode input creation tool
    
    This tools creates site and landform netCDF files and a gridlist file
    from SRTM1 (or potentially other) elevation data.
     
    """
    
    # example:
    #./lgt_createinput.py processed output --dems=srtm1/*.zip --masks=srtm1_shp_mask --extent -76 -56 -66 -16

    if verbose:
        logging.getLogger(__name__).setLevel(logging.DEBUG)
    else:
        logging.getLogger(__name__).setLevel(logging.INFO)
        
    if dems is not None:
        SRTMSTORE_PATH = dems
    
    if masks is not None:
        WATERMASKSTORE_PATH = masks

    REGION = None
    if len(extent) == 4:
        REGION = list(extent)
        
    
    # the setup dictionary to convert into a bunch obj
    config_data=dict(SRTMSTORE_PATH=dems,
                     WATERMASKSTORE_PATH=masks,
                     TILESTORE_PATH=storage,
                     REGION=REGION,
                     CLASSIFICATION=classfication,
                     CUTOFF=cutoff,
                     OUTPUT_PATH=outdir,
                     GRIDLIST_TXT=gridlist,
                     OUTDIR=outdir,
                     OVERWRITE=force_overwrite)
    
    # TODO: change logging level based on verbose flag
    cfg = Bunch(config_data)

    main(cfg)
 
