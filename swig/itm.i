%module itm
%include <std_string.i>

%{
extern void wgs842itm(double lat,double lon,int& N,int& E);
extern void itm2wgs84(int N,int E,double& lat,double &lon);
%}


%inline {

struct ITM {
    int N;
    int E;
};

typedef struct {
    double lat;
    double lng;
} WGS84;

ITM
WGS84toITM(WGS84 wgs) {
    ITM itm;

    wgs842itm(wgs.lat, wgs.lng, itm.N, itm.E);
    return itm;
}

WGS84
ITMtoWGS84(ITM itm) {
    WGS84 wgs;

    itm2wgs84(itm.N, itm.E, wgs.lat, wgs.lng);
    return wgs;
}

}
