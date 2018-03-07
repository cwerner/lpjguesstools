# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
from setuptools.command.sdist import sdist as SDistCommand

from codecs import open
from os import path
import re

version = '0.0.2'
#version = re.search('^__version__\s*=\s*"(.*)"',
#                    open('lpjguesstools/lgt_convert/lgt_convert.py').read(), re.M).group(1)

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
#with open("README.md", "rb") as f:
#    long_descr = f.read().decode("utf-8")

# Auto convert md to rst for PyPI
readme_file = 'README.md'
try:
    from pypandoc import convert
    long_descr = convert(readme_file, 'rst', 'md')
    with open(path.join(here, 'README.rst'), 'w', encoding='utf-8') as f:
        f.write(long_descr)
except ImportError:
    print(
        "warning: pypandoc module not found, could not convert Markdown to RST")
    long_descr = open(readme_file).read()

# get the dependencies and installs
with open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
dependency_links = [x.strip().replace('git+', '') for x in all_reqs
                    if 'git+' not in x]


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

class PyPack(SDistCommand):
    description = "Custom sdist command, runs pandoc prior to building the sdist package"
    def run(self):
        try:
            from pypandoc import convert
            long_descr = convert(readme_file, 'rst', 'md')
            with open(path.join(here, 'README.rst'), 'w', encoding='utf-8') as f:
                f.write(long_descr)
        except ImportError:
            print("warning: pypandoc not found, could not convert Markdown to RST")
            long_descr = open(readme_file).read()
        SDistCommand.run(self)
        

setup(name='lpjguesstools',
      version=version,
      description='This package contains pre- and postprocessing tools for LPJ-GUESS 4.0',
      long_description=long_descr,
      url='https://gitlab.com/cw_code/lpjguesstools',
      author='Christian Werner',
      author_email='christian.werner@senckenberg.de',
      license='ND',
      download_url='https://gitlab.com/cw_code/lpjguesstools/tarball/' + version,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: LPJ-GUESS scientists',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 3',
      ],
      keywords='LPJ-GUESS preprocessing postprocessing netcdf',
      zip_safe= False,
      packages=find_packages(exclude=['docs', 'tests']),
      install_requires=install_requires,
      extras_require={'test': ['pytest'], },
      package_data={'lpjguesstools': ['data/lgt_convert.conf', 
                                      'data/GLOBAL_WISESOIL_DOM_05deg.nc',
                                      'data/GLOBAL_ELEVATION_05deg.nc']},
      include_package_data=True,
      entry_points={'console_scripts': [
          'lgt_biomize=lpjguesstools.lgt_biomize.cli:cli',
          'lgt_convert=lpjguesstools.lgt_convert.cli:cli',
          'lgt_createinput=lpjguesstools.lgt_createinput.cli:cli']},
      test_suite='lpjguess2nc.test.test_lpjguess2nc',
      cmdclass={'test': PyTest, 'sdist': PyPack},
      dependency_links=dependency_links)
