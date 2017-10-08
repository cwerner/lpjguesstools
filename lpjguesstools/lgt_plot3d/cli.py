# -*- coding: utf-8 -*-
"""lgt_plot3d.cli: Commandline interpreter for lgt_plot3d.py."""

import click
import logging
import xarray as xr

from main import main

# import constants
from .. import EPILOG


log = logging.getLogger(__name__)

class Bunch(object):
    """Simple data storage class."""
    def __init__(self, adict):
        self.__dict__.update(adict)
 

def validate_variable(ctx, param, value):
    """Check that requested variables are in netCDF file """
    

    for p in ctx.command.params:

        
        #print p
        if isinstance(p, click.Argument) and p.name == 'tile':
            print '>', p.name
            tile = help(p.process_value)
            print tile
            with xr.open_dataset(tile) as ds:
                for v in value:
                    if v not in ds.data_vars:
                        raise click.BadParameter('Variable %s not in specified file %s.' % (v, ctx.args.tile))

    return value



class CommandValidateTile(click.Command):

    def make_context(self, *args, **kwargs):
        ctx = super(CommandValidateTile, self).make_context(*args, **kwargs)
        self._validate_variable(ctx)
        return ctx

    def _validate_variable(self, ctx):
        """Check that requested variable(s) are in netCDF file """

        value = ctx.params['variable']
        tile = ctx.params['tile']
        with xr.open_dataset(tile) as ds:
            missing = []
            for v in value:
                if v not in ds.data_vars.keys():
                    missing.append(v)
            if len(missing) > 0:
                #click.echo(click.Context.get_help(ctx))
                click.echo("")
                available = 'Available vars ' + ', '.join(sorted(ds.data_vars.keys()))
                raise click.UsageError('Var "%s" not in file %s.\n\n%s' % (','.join(missing), tile, available))

# command line arguments
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
click.Context.get_usage = click.Context.get_help

@click.command(cls=CommandValidateTile, context_settings=CONTEXT_SETTINGS, epilog=EPILOG)

@click.option('--screen', is_flag=True, default=False,
                    help='render on screen [NOT FUNTIONAL ATM]')

@click.option('--variable', '-v', multiple=True, metavar='VAR', 
                    help='variable to render')

@click.option('--verbose', is_flag=True, 
                    help='increase logging info')

@click.version_option()

@click.argument('tile', type=click.Path(exists=True))
@click.argument('plot') #, type=click.File()) 



def cli(tile, plot, screen, variable, verbose):
    """LPJ-GUESS 4.0 3D tile plotter
    
    This tools creates 3d plots for model output or site topography.
     
    """
    
    # example:
    #./lgt_plot3d.py -v SoilC -V mwcount_upper processed/srtm1_processed_xyz.nc myplot.png 

    if verbose:
        logging.getLogger(__name__).setLevel(logging.DEBUG)
    else:
        logging.getLogger(__name__).setLevel(logging.INFO)

    # the setup dictionary to convert into a bunch obj
    config_data=dict(TILE=tile,
                     OUTFILE=plot,
                     SHOW=screen,
                     VARS=list(variable))
    
    # TODO: change logging level based on verbose flag
    cfg = Bunch(config_data)

    print cfg.TILE

    main(cfg)
 