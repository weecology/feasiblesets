from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("combine_macrostates", ["combine_macrostates.pyx"])]

setup(
  name = 'combine_macrostates',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

""" This is a Cython setup file for the Python script (.pyx extension) combine_macrostates.pyx """
