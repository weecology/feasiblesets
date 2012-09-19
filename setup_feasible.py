from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("feasible_functions", ["feasible_functions.pyx"])]

setup(
  name = 'feasible_functions',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

""" This is a Cython setup file for the Python script (.pyx extension) feasible_functions.pyx """