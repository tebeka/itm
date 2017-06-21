# distutils: language=c++
# distutils: sources = isr84lib.cc
"""Convert Lat/lng to/from ITM (Israeli Transverse Mercator)"""

cimport decl

__version__ = '0.1.1'

def wgs842itm(double lat, double lng):
    """Convert WGS84 lat/lng to ITM N/E"""
    cdef int N, E

    N = E = 0  # Make the compilter happy

    decl.wgs842itm(lat, lng, N, E)
    return N, E


def itm2wgs84(int N, int E):
    """Convert ITM N/E to lat/lng"""
    cdef double lat, lng

    lat = lng = 0.0  # Make the compilter happy

    decl.itm2wgs84(N, E, lat, lng)

    return lat, lng
