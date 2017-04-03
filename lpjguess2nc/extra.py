# -*- coding: utf-8 -*-
"""lpjguess2nc.extra: extra module within the lpjguess2nc package."""

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
        Requirement.parse("ldndc2nc"), "ldndc2nc/data/ldndc2nc.conf")
    shutil.copyfile(fname, os.path.join(
        os.path.expanduser("~"), "ldndc2nc.conf"))


def _find_config():
    """ look for cfgFile in the default locations """
    cfgFile = None
    locations = [os.curdir, os.path.expanduser("~"), "/etc/ldndc2nc",
                 os.environ.get("LDNDC2NC_CONF")]
    locations = [x for x in locations if x is not None]

    for loc in locations:
        f = os.path.join(loc, "ldndc2nc.conf")
        if os.path.isfile(f):
            cfgFile = f
            break

    return cfgFile


def _parse_config(cfgFile):
    """ read yaml config file and modify special properties"""

    with open(cfgFile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    for k, vs in cfg['variables'].items():
        vs_new = []
        for v in vs:

            def is_multipart_item(x):
                return ';' in x

            if is_multipart_item(v):
                x = v.split(';')
                vs_new.append((x[0], x[1:]))
            else:
                vs_new.append(v)

            cfg['variables'][k] = vs_new

    return cfg


def parse_config(cfg, section=None):
    """ parse config data structure, return data of required section """

    def is_valid_section(s):
        valid_sections = ['info', 'project', 'variables', 'refdata']
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

    cfg = _parse_config(cfgFile)

    return cfg


def set_config(cfg):
    """ write cfg file to user dir """
    fname = os.path.join(os.path.expanduser("~"), 'ldndc2nc.conf')
    with open(fname, 'w') as f:
        f.write(yaml.dump(cfg, default_flow_style=False))


class RefDataBuilder(param.Parameterized):
    lonmin = param.Number(None, bounds=(-180, 180), doc="min longitude")
    latmin = param.Number(None, bounds=(-90, 90), doc="min latitude")
    lonmax = param.Number(None, bounds=(-180, 180), doc="max longitude")
    latmax = param.Number(None, bounds=(-90, 90), doc="max latitude")
    local = param.Boolean(False, doc="number ids for regional subset")
    formula = param.String(default='continuous')
    res = param.Number(None, bounds=(0, 5), doc="cell resolution in degrees")
    i_shift = 0
    j_shift = 0

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

        try:
            self.local = _cfg['local']
        except Exception as e:
            log.debug("No <local> statement in refdata: using 'global'")

        try:
            formula = _cfg['formula']
            formula = self._check_formula(formula)
            self.formula = formula
        except:
            log.debug("No <formula> statement in refdata: using 'continuous'")

        # we work with cell centers
        cell_half = self.res * 0.5
        self.lons = np.arange(self.lonmin + cell_half, self.lonmax, self.res)
        self.lats = np.arange(self.latmin + cell_half, self.latmax, self.res)
        self.globlons = np.arange(-180 + cell_half, 180, self.res)
        self.globlats = np.arange(-90 + cell_half, 90, self.res)

        if not self.local:
            # compute shift of local bbox in respect to global domain
            m_lon = self._find_nearest(self.globlons, self.lons[0])
            m_lat = self._find_nearest(self.globlats, self.lats[::-1][0])
            self.i_shift = np.where(self.globlons == m_lon)[0]
            self.j_shift = np.where(self.globlats[::-1] == m_lat)[0]

    def _check_formula(formula_str):
        """ check norefdata formula given in conf file """
        valid_chars = "xyij0123456789+*^"
        safe_method = []
        formula_str = string.lower(string.replace(formula_str, '^', '**'))
        for c in formula_str:
            if c in valid_chars:
                formula_str_validated.append(c)
        return ''.join(formula_str_validated)

    def _compute_formula_cid(self, j, i):
        """ calculate cellid based on formula """
        i += self.i_shift
        j += self.j_shift
        return eval(self.rd.formula, {'__builtins__': None}, {})

    def _compute_continuous_cid(self, j, i):
        """ calculate continuous cellid """
        if self.local:
            i_len = len(self.lons)
        else:
            i_len = len(self.globlons)
        return (j + self.j_shift) * i_len + (i + self.i_shift)

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
