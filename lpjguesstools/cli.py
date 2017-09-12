# -*- coding: utf-8 -*-
#
# cli.py
# ==================
#
# Christian Werner
# christian.werner@senckenberg.de

import argparse
import datetime as dt
import logging
import os
import pkg_resources

version = pkg_resources.require("lpjguesstools")[0].version

log = logging.getLogger(__name__)


EPILOG = """
    Use this tool to create netCDF files based on standard
    LPJ-GUESS 4 txt output files
    """

DESCR = "lpt_convert :: LPJ-GUESS output converter (v%s)\n" % version


class VerbosityAction(argparse.Action):
    """ CustomAction for argparse that increases the log level """
    def __init__(self, nargs=0, **kw):
        argparse.Action.__init__(self, nargs=nargs, **kw)

    def __call__(self, parser, namespace, values=None, option_string=None):
        handlers = logging.getLogger().handlers
        for handler in handlers:
            if type(handler) is logging.StreamHandler:
                handler.setLevel(logging.DEBUG)
        setattr(namespace, self.dest, True)


class RangeAction(argparse.Action):
    """ CustomAction for argparse to be able to process value range """
    def __call__(self, parser, namespace, values, option_string=None):
        s = values.split('-')
        def is_valid_year_range(s):
            if len(s) > 2:
                return False
            for e in s:
                try:
                    _ = int(e)
                except:
                    return False
            if int(s[1]) < int(s[0]):
                return False
            return True

        if is_valid_year_range(s):
            setattr(namespace, self.dest, range(int(s[0]), int(s[-1])+1))
        else:
            log.critical('No valid range: %s' % values)
            exit(1)

class MultiArgsAction(argparse.Action):
    """ CustomAction for argparse to be able to process ,-seperated args """
    def __init__(self, **kw):
        argparse.Action.__init__(self, **kw)
        self._nsegs = self.const     # misuse of const to carry len(segs)

    def __call__(self, parser, namespace, values, option_string=None):
        s = values.split(',')
        if len(s) == self._nsegs:
            setattr(namespace, self.dest, tuple(s))
        elif len(s) == 1:
            setattr(namespace, self.dest, tuple([s[0], None]))
        else:
            log.critical("Syntax error in %s" % option_string)
            exit(1)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    if type(action.default)==list:
                        help += ' (default: %d-%d)' % ( action.default[0], action.default[-1] )
                    else:
                        help += ' (default: %(default)s)'
        return help


def cli():
    """ command line interface """

    DESCR = "lpt_convert :: LPJ-GUESS output converter (v%s)" % version

    GREETING = '\n'.join(['-'*78, DESCR, '-'*78])

    EPILOG  = "Use this tool to create netCDF files based on standard\n"
    EPILOG += "LPJ-GUESS 4 txt output files\n"


    parser = argparse.ArgumentParser(description=DESCR,
                                     epilog=EPILOG,
                                     formatter_class=CustomFormatter)

    parser.add_argument('indir', help="location of lpjguess txt or txt.gz files")
    parser.add_argument('outname', help="output netCDF file")

    # TODO: fully implement the avg. functionality
    parser.add_argument("-a",
                    dest="avg",
                    action='store_true',
                    default=False,
                    help="average output years (see -y)")

    parser.add_argument("-c",
                    dest="config",
                    metavar='MYCONF',
                    help="use MYCONF file as config (disabled)")

    # TODO: not yet implemented
    parser.add_argument("-l",
                    dest="last_nyears",
                    default=-1,
                    type=int,
                    metavar='NYEARS',
                    help="number of years to use (or specify range with -y)")

    parser.add_argument("-m",
                    dest="use_month_dim",
                    default=False,
                    action='store_true',
                    help="create extra month dimension <time, month, lat, lon>")

    # TODO: implemented, but needs more flexibility (allow any netCDF file etc.)
    parser.add_argument("-r",
                    dest="refinfo",
                    action=MultiArgsAction,
                    const=2,
                    metavar="FILE,VAR",
                    help="refdata from (landforms) netCDF file")

    parser.add_argument("-S",
                    dest="storeconfig",
                    action='store_true',
                    default=False,
                    help="make passed config (-c) the new default")

    parser.add_argument("-v",
                    dest="verbose",
                    action=VerbosityAction,
                    default=False,
                    help="increase output verbosity")

    parser.add_argument("-y",
                    dest="years",
                    default=range(1950,1989),
                    action=RangeAction,
                    help="range of years to consider")

    # extra subpixel/ landform switches
    parser.add_argument("--ns",
                    dest="north_south",
                    default=False,
                    action='store_true',
                    help="create N, S aspect data for landform variables")

    print GREETING

    args = parser.parse_args()


    log.debug('-' * 50)
    log.debug('lpt_convert called at: %s' % dt.datetime.now())

    if args.storeconfig and (args.config is None):
        log.critical("Option -S requires that you pass a file with -c.")
        exit(1)

    if args.years != range(1950,1989) and args.last_nyears != -1:
        log.critical("Use either option -y or Option -l.")
        exit(1)

    return (parser, args)
