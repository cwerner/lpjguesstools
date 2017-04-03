LPJ-GUESS 4 output postprocessor to netCDF format
======================================================

version number: 0.0.1 author: Christian Werner
(christian.werner@senckenberg.de)

Overview
--------

This package converts LPJ-GUESS 4 output to netCDF files for selected
variables. Currently used for EarthShape project, but in general
applicable for all LPJ-GUESS simulations.

Installation / Usage
--------------------

To install use pip:

::

    $ pip install lpjguess2nc

Or clone the repo:

::

    $ git clone https://gitlab.com/cw_code/lpjguess2nc.git
    $ python setup.py install

Contributing
------------

TBD

Example
-------

lpjguess2nc -r REFDATA.nc,cid -y 1990-2014 -o output.nc lpjguess\_results\_dir
lpjguess\_netcdf\_dir

Usage
-----

::

    usage: lpjguess2nc [-h] [-c MYCONF] [-l PATTERN] [-o OUTFILE] [-r FILE,VAR]
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
