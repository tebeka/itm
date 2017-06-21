from setuptools import setup
from Cython.Build import cythonize


def find_version():
    with open('itm.pyx') as fp:
        for line in fp:
            if '__version__' not in line:
                continue
            _, version = line.split(' = ')
            return version.replace("'", '').strip()
    assert False, 'cannot find version'


setup(
    name='itm',
    version=find_version(),
    description='ITM <-> WGS84 conversions',
    long_description=open('README.md').read(),
    author='Miki Tebeka',
    author_email='miki.tebeka@gmail.com',
    license='MIT',
    url='https://github.com/tebeka/fastavro',
    ext_modules=cythonize('itm.pyx', language='c++'),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
