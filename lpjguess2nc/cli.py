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

version = pkg_resources.require("lpjguess2nc")[0].version

#version = "0.0.1"

log = logging.getLogger(__name__)


EPILOG = """
    Use this tool to create netCDF files based on standard
    LPJ-GUESS 4 txt output files
    """

DESCR = "lpjguess2nc :: LPJ-GUESS output converter (v%s)\n" % version


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

    DESCR = "lpjguess2nc :: LPJ-GUESS output converter (v%s)" % version

    GREETING = '\n'.join(['-'*78, DESCR, '-'*78])

    EPILOG  = "Use this tool to create netCDF files based on standard\n"
    EPILOG += "LPJ-GUESS 4 txt output files\n"


    parser = argparse.ArgumentParser(description=DESCR,
                                     epilog=EPILOG,
                                     formatter_class=CustomFormatter)

    parser.add_argument('indir', help="location of source lpjguess txt files")
    parser.add_argument('outdir', help="destination of created netCDF files")

    parser.add_argument("-c",
                    dest="config",
                    metavar='MYCONF',
                    help="use MYCONF file as config")

    parser.add_argument("-l",
                    dest="limiter",
                    metavar='PATTERN',
                    help="limit files by PATTERN")

    parser.add_argument("-o",
                    dest="outfile",
                    default="outfile.nc",
                    help="name of the output netCDF file")

    parser.add_argument(
                    "-r",
                    dest="refinfo",
                    action=MultiArgsAction,
                    const=2,
                    metavar="FILE,VAR",
                    help="refdata from netCDF file")

    #parser.add_argument("-s",
    #                dest="split",
    #                action='store_true',
    #                default=False,
    #                help="split output to yearly netCDF files")

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
                    default=range(1990,2014),
                    action=RangeAction,
                    help="range of years to consider")

    print GREETING

    args = parser.parse_args()


    log.debug('-' * 50)
    log.debug('lpjguess2nc called at: %s' % dt.datetime.now())

    if args.storeconfig and (args.config is None):
        log.critical("Option -S requires that you pass a file with -c.")
        exit(1)

    return args
