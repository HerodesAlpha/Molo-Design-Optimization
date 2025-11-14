import sys
import subprocess
import os
import glob
import shutil

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

##
import importlib
pyNemoh_spec = importlib.util.find_spec("pyNemoh")
pyNemoh_found = pyNemoh_spec is not None

if pyNemoh_found:
    print(pyNemoh_spec)
##
__author__ = "Eivind Sonju"
__copyright__ = "Copyright (C) 2017-2019 Verbun AS. All rights reserved."
__version__ = "2.0"

args = sys.argv[1:]

# Make a `cleanall` rule to get rid of intermediate and library files
if "cleanall" in args:
    print("Deleting cython files...")
    # Just in case the build directory was created by accident,
    for name in glob.glob('*.so'):
        os.remove(name)
    for name in glob.glob('*.pyd'):
        os.remove(name)
    for name in glob.glob('*.pyc'):
        os.remove(name)
    shutil.rmtree('build', True)
    for name in glob.glob('*.c'):
        os.remove(name)

    # Now do a normal clean
    sys.argv[1] = "clean"




# Get the path to the build directory where libnemoh.dll is located
build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'build'))
pyNemoh_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'pyNemoh'))

ext_modules = [Extension(language='c',
                         name='pyNemoh.solver_fortran',
                         sources=['pyNemoh/solver_fortran.pyx'],
                         libraries=['libnemoh'],
                         library_dirs=[build_dir, pyNemoh_dir],
                         include_dirs=[np.get_include()]
                         ),
               ]

setup(
    name='pyNemoh',
    version='2.0',
    description='An adatation of OpenWARP CLI to MOLO',
    packages=['pyNemoh'],
    url='',
    license='Apache License 2.0',
    author='Eivind Sonju',
    author_email='es@verbun.com',
    data_files=[(r'Lib\site-packages\pyNemoh', [r'pyNemoh\libnemoh.dll', r'pyNemoh\libgfortran-3.dll'])],
    #zip_safe=False,
    ext_modules=cythonize(ext_modules),
)
