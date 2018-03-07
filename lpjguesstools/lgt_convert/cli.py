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

@click.option('-a', 'avg', is_flag=True, default=False,
                    help='average output years (see -y or -l)')

@click.option('-c', 'config', metavar='MYCONF', default=None,
                    help='use MYCONF file as config')

@click.option('-d', 'default', is_flag=True, default=False,
                    help='only process default stand (no landforms)')

@click.option('-l', 'last_nyears', type=click.INT, default=None, metavar='NYEARS',
                    help='use number of years at end of run (alternative to -y)')

@click.option('-m', 'use_month_dim', is_flag=True, default=False, 
                    help='create extra month dimension <time, month, lat, lon>')

@click.option('-r', 'refinfo', metavar='FILE,VAR', default=None,
                    help='refdata from (landforms,variable) netCDF file')
 
@click.option('-s', 'smode', is_flag=True, default=False,
                    help='run in site mode (only single coordinate)')

@click.option('-S', 'storeconfig', is_flag=True, default=False,
                    help='make passed config (-c) the new default')

@click.option('-y', 'years', metavar="1900-1989", default=[],
                    help='select range of calendar years (alternative to -l)')
                    
@click.option('--verbose', is_flag=True, 
                    help='increase logging info')


@click.version_option()

@click.argument('indir', type=click.Path(exists=True))
@click.argument('outname') 

def cli(avg, config, default, last_nyears, use_month_dim, refinfo, smode, storeconfig, years, indir, outname, verbose):
    """LPJ-GUESS 4.0 subpixel mode netCDF convert tool
    
    This tools nonverts default and subpixel mode output from
    LPJ-GUESS 4.0 .out (or gzipped .out.gz) files and creates
    a joined netCDF file.
     
    """
    
    # example:
    #./lgt_createinput.py processed output --dems=srtm1/*.zip --masks=srtm1_shp_mask --extent -76 -56 -66 -16

    if verbose:
        logging.getLogger(__name__).setLevel(logging.DEBUG)
    else:
        logging.getLogger(__name__).setLevel(logging.INFO)

    if years != []:
        try:
            years=[int(x) for x in years.split('-')]
        except:
            log.error("Years require format 1990-1999")
            exit(1)
    
    if refinfo != None:        
         if ',' in refinfo:
             f,v=[x.strip() for x in refinfo.split(',')[0:2]]
         else:
             f=refinfo.strip()
             v=None
         refinfo=(f,v)
    else:
         refinfo=(None,None)
 
 
    if smode and (use_month_dim == False):
         log.warn("Site mode requires the use of a monthly dim (-m). Proceeding.")
         use_month_dim = True
 
    if avg and (use_month_dim == False):
         log.warn("Average mode requires the use of a monthly dim (-m). Proceeding.")
         use_month_dim = True
 
    # the setup dictionary to convert into a bunch obj
    config_data=dict(AVG=avg,
                     CONFIG=config,
                     DEFAULT=default,
                     LAST_NYEARS=last_nyears,
                     USE_MONTH_DIM=use_month_dim,
                     REFINFO=refinfo,
                     SMODE=smode,
                     STORECONFIG=storeconfig,
                     YEARS=years,
                     INDIR=indir,
                     OUTNAME=outname)
    
    # TODO: change logging level based on verbose flag
    cfg = Bunch(config_data)


    if cfg.STORECONFIG and (cfg.CONFIG is None):
        log.critical("Option -S requires that you pass a file with -c.")
        exit(1)

    if cfg.YEARS != [] and cfg.LAST_NYEARS != -1:
        log.critical("Use either option -y or Option -l.")
        exit(1)

    main(cfg)
