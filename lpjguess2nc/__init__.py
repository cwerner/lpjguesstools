import logging
from logging.handlers import RotatingFileHandler
import sys

try:  # python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

# TODO make these flexible (ENV var and/ or lpjguess2nc.conf) 
logPath = '.'
fileName = 'lpjguess2nc'

class MultiLineFormatter(logging.Formatter):
    """ A custom multi-line logging formatter """

    def format(self, record):
        str = logging.Formatter.format(self, record)
        header, footer = str.split(record.message)
        str = str.replace('\n', '\n' + ' ' * len(header))
        return str


CONS_FORMAT = "[%(levelname)-8s] %(message)s"
FILE_FORMAT = "%(asctime)s [%(levelname)-8s] %(message)s (%(filename)s:%(lineno)s)"

lfCons = MultiLineFormatter(CONS_FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
lfFile = MultiLineFormatter(FILE_FORMAT, datefmt='%Y-%m-%d %H:%M:%S')

rootLogger = logging.getLogger()
rootLogger.setLevel(logging.DEBUG)

hCons = logging.StreamHandler()
hCons.setFormatter(lfCons)
hCons.setLevel(logging.INFO)
rootLogger.addHandler(hCons)

hFile = RotatingFileHandler("{0}/{1}.log".format(logPath, fileName), maxBytes=10000)
hFile.setFormatter(lfFile)
hFile.setLevel(logging.DEBUG)
rootLogger.addHandler(hFile)
