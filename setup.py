from distutils.core import setup
from Cython.Build import cythonize

setup(
    name='itm',
    version='0.1.0',
    description='ITM <-> WGS84 conversions',
    ext_modules=cythonize('itm.pyx', language='c++'),
)
