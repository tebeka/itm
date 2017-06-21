from distutils.core import setup, Extension


setup(
    name='itm',
    version='0.1.0',
    description='Lat/lng to ITM',
    ext_modules=[
        Extension(
            '_itm',
            sources=['itm.i', 'isr84lib.cc'],
            swig_opts=['-c++'],
        ),
    ],
)
