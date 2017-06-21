# distutils: language=c++

cdef extern from "isr84lib.h":
    void wgs842itm(double lat,double lon,int& N,int& E);
    void itm2wgs84(int N,int E,double& lat,double &lon);
