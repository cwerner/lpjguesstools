LPJGUESSTOOLS: LPJ-GUESS4 pre- and postprocessing tools
=======================================================

version number: 0.0.1 author: Christian Werner
(<christian.werner@senckenberg.de>)

Overview
--------

This package bundles a set of helper scripts. For now, these include a script to produce
LPJ-GUESS 4 (subpixel) input files and a converter to produce netcdf files from model output.
Currently used for EarthShape project, but in general
applicable for all LPJ-GUESS simulations.

Installation / Usage
--------------------

To install use pip:

::

    $ pip install lpjguesstools

Or clone the repo:

::

    $ git clone https://gitlab.com/cw_code/lpjguesstools.git
    $ python setup.py install

Contributing
------------

TBD

Example (lgt_convert)
---------------------

lgt_convert -r REFDATA.nc,cid -y 1990-2014 -o output.nc lpjguess\_results\_dir
lpjguess\_netcdf\_dir

Usage
-----

    usage: lgt_convert [-h] [-c MYCONF] [-l PATTERN] [-o OUTFILE] [-r FILE,VAR]
                       [-S] [-v] [-y YEARS]
                       indir outdir

    positional arguments:
      indir        location of source ldndc txt files
      outdir       destination of created netCDF files

    optional arguments:
      -h, --help   show this help message and exit
      -c MYCONF    use MYCONF file as config (default: None)
      -l PATTERN   limit files by PATTERN (default: None)
      -o OUTFILE   name of the output netCDF file (default: outfile.nc)
      -r FILE,VAR  refdata from netCDF file (default: None)
      -S           make passed config (-c) the new default (default: False)
      -v           increase output verbosity (default: False)
      -y YEARS     range of years to consider (default: 1990-2014)

