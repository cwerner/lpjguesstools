# -*- coding: utf-8 -*-
"""lpjguesstools.extra: extra module within the lpjguesstools package."""

import logging
import os
from pkg_resources import Requirement, resource_filename
import shutil
import string

import numpy as np
import param
import yaml

log = logging.getLogger(__name__)


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def _copy_default_config():
    """ copy default conf file to user dir """

    #TODO somewhat redundand, merge with set_config code

    fname = resource_filename(
        Requirement.parse("lpjguesstools"), "lpjguesstools/data/lgt_convert.conf")
    shutil.copyfile(fname, os.path.join(
        os.path.expanduser("~"), "lgt_convert.conf"))


def _find_config():
    """ look for cfgFile in the default locations """
    cfgFile = None
    locations = [os.curdir, os.path.expanduser("~"), "/etc/lpjguesstools",
                 os.environ.get("LPJGUESSTOOLS_CONF")]
    locations = [x for x in locations if x is not None]

    for loc in locations:
        f = os.path.join(loc, "lgt_convert.conf")
        if os.path.isfile(f):
            cfgFile = f
            break

    return cfgFile


def _parse_config(cfgFile):
    """ read yaml config file and modify special properties"""

    with open(cfgFile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    def is_multipart_item(x):
        return ';' in x

    def is_assigned_item(x):
        return '=' in x

    for ix, v in enumerate(cfg['data']):
        vs_new = []
        if is_multipart_item(v):
            xs = v.split(';')
            d = {}
            for x in xs[1:]:
                if is_assigned_item(x):
                    _k, _v = x.split('=')
                    d[_k] = _v
                else:
                    d[x] = x
            # return tuple file: variables, variables(renamed)
            vs_new.append((xs[0], d))
        else:
            # return just a list of items
            vs_new.append(v)

        cfg['data'][ix] = vs_new

    return cfg


def parse_config(cfg, section=None):
    """ parse config data structure, return data of required section """

    def is_valid_section(s):
        valid_sections = ['info', 'project', 'data', 'refdata', 'insfile_keywords']
        return s in valid_sections

    cfg_data = None
    if is_valid_section(section):
        try:
            cfg_data = cfg[section]
        except KeyError:
            log.critical(cfg.keys())
            log.critical("Section <%s> not found in config" % section)
            exit(1)
    else:
        log.critical("Section <%s> not a valid name" % section)
        exit(1)
    return cfg_data


def get_config(cfgFile=None):
    """ locate and read config file """

    cfg = None
    locations = []

    def cfgfile_exists(cfgFile):
        return cfgFile != None

    if cfgfile_exists(cfgFile):
        if not os.path.isfile(cfgFile):
            log.critical("Specified configuration file not found.")
            exit(1)
    else:
        cfgFile = _find_config()
        if not cfgfile_exists(cfgFile):
            log.info("Copying config file")
            _copy_default_config()
            cfgFile = _find_config()
    
    log.debug("Using conf file: %s" % cfgFile)
    cfg = _parse_config(cfgFile)

    return cfg


def set_config(cfg):
    """ write cfg file to user dir """
    fname = os.path.join(os.path.expanduser("~"), 'lgt_convert.conf')
    with open(fname, 'w') as f:
        f.write(yaml.dump(cfg, default_flow_style=False))


class RefDataBuilder(param.Parameterized):
    lonmin = param.Number(None, bounds=(-180, 180), doc="min longitude")
    latmin = param.Number(None, bounds=(-90, 90), doc="min latitude")
    lonmax = param.Number(None, bounds=(-180, 180), doc="max longitude")
    latmax = param.Number(None, bounds=(-90, 90), doc="max latitude")
    res = param.Number(None, bounds=(0, 5), doc="cell resolution in degrees")

    def __init__(self, cfg):
        _cfg = parse_config(cfg, section='refdata')

        try:
            self.lonmin, self.latmin, self.lonmax, self.latmax = _cfg['bbox']
        except Exception as e:
            log.critical(str(e))

        try:
            self.res = _cfg['res']
        except Exception as e:
            log.critical("No <res> statement in refdata")
            exit(1)

        # we work with cell centers
        cell_half = self.res * 0.5
        self.lons = np.arange(self.lonmin + cell_half, self.lonmax, self.res)
        self.lats = np.arange(self.latmin + cell_half, self.latmax, self.res)
        self.globlons = np.arange(-180 + cell_half, 180, self.res)
        self.globlats = np.arange(-90 + cell_half, 90, self.res)

    def _find_nearest(self, array, value):
        """ locate closest value match """
        idx = (np.abs(array - value)).argmin()
        return array.flat[idx]

    def build(self):
        """ actually populate cellid array """
        cell_ids = np.empty((len(self.lats), len(self.lons)), np.int64)
        for j, lat in enumerate(self.lats):
            for i, lon in enumerate(self.lons):
                if self.formula != 'continuous':
                    cid = self._compute_formula_cid(j, i)
                else:
                    cid = self._compute_continuous_cid(j, i)
                cell_ids[j, i] = cid
        return (cell_ids, self.lats, self.lons)


def create_refdata(cfg):
    """ produce refdata using info from cfg file
        :param: yaml cfg
        :return: cell_ids, lats, lons
        :rtype: tuple
    """
    rdb = RefDataBuilder(cfg)
    cell_ids = rdb.build()
