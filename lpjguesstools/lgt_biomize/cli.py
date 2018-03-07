# -*- coding: utf-8 -*-
"""lgt_biomize.cli: Commandline interpreter for lgt_biomize.py."""

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

@click.option('--classification', type=click.Choice(['CHILE_ES_NEW']), 
                    default='CHILE_ES_NEW', show_default=True,
                    help='classification scheme')

@click.option('--limit', is_flag=True, default=False, 
                    help='limit site mode to 1500 yrs for debugging')

@click.option('--site-mode', '-s', is_flag=True, default=False, 
                    help='apply to single site netCDF file (also no vertical netCDF)')

@click.option('--verbose', is_flag=True, 
                    help='increase logging info')

@click.version_option()

@click.argument('infile', type=click.Path(exists=True))
@click.argument('outfile')

def cli(classification, limit, site_mode, infile, outfile, verbose):
    """LPJ-GUESS 4.0 subpixel mode biomization tool
    
    This tools creates biomization netCDF file(s)
    from lgt_convert produced output netCDFs.
    Spatial (2D) and single site (1D) files are accepted.
     
    """
    
    if verbose:
        logging.getLogger(__name__).setLevel(logging.DEBUG)
    else:
        logging.getLogger(__name__).setLevel(logging.INFO)
        
    # the setup dictionary to convert into a bunch obj
    config_data=dict(INFILE=infile,
                     OUTFILE=outfile,
                     SMODE=site_mode,
                     LIMIT=limit,
                     CLASSIFICATION=classification)
    
    # TODO: change logging level based on verbose flag
    cfg = Bunch(config_data)

    main(cfg)
 
