import logging
from logging.handlers import RotatingFileHandler
import numpy as np
import sys

try:  # python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

# TODO: create sublogger for different scripts
logPath = '.'
fileName = 'lpjguesstools'

class MultiLineFormatter(logging.Formatter):
    """ A custom multi-line logging formatter """

    def format(self, record):
        str = logging.Formatter.format(self, record)
        header, footer = str.split(record.message)
        str = str.replace('\n', '\n' + ' ' * len(header))
        return str

# optional colored console logger (nice!)
try:
    import colorlog
    class MultiLineFormatterColor(colorlog.ColoredFormatter):
        def format(self, record):
            record.__dict__.update(colorlog.escape_codes)
            record.log_color = self.color(self.log_colors, record.levelname)

            str = logging.Formatter.format(self, record)
            header, footer = str.split(record.message)
            str = str.replace('\n', '\n' + ' ' * len(header))
            return str
    CONS_FORMAT = "[%(log_color)s%(levelname)-8s%(reset)s] %(log_color)s%(message)s%(reset)s"

except ImportError:
    # both formatters should use the default (non-color)
    MultiLineFormatterColor = MultiLineFormatter
    CONS_FORMAT = "[%(levelname)-8s] %(message)s"


FILE_FORMAT = "%(asctime)s [%(levelname)-8s] %(message)s (%(filename)s:%(lineno)s)"

lfCons = MultiLineFormatterColor(CONS_FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
lfFile = MultiLineFormatter(FILE_FORMAT, datefmt='%Y-%m-%d %H:%M:%S')

rootLogger = logging.getLogger(__name__)
rootLogger.setLevel(logging.DEBUG)

hCons = logging.StreamHandler()
hCons.setFormatter(lfCons)
hCons.setLevel(logging.DEBUG)
rootLogger.addHandler(hCons)

hFile = RotatingFileHandler("{0}/{1}.log".format(logPath, fileName), maxBytes=10000)
hFile.setFormatter(lfFile)
hFile.setLevel(logging.DEBUG)
rootLogger.addHandler(hFile)

EPILOG = """Christian Werner, SENCKENBERG Biodiversity and Climate Research Centre (BiK-F)
email: christian.werner@senkenberg.de
2017/09/26"""
