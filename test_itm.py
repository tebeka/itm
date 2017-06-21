import itm


def test_wgs842itm():
    assert itm.wgs842itm(32.3, 35.7) == (689686, 266124)


def test_itm2wgs84():
    lat, lng = itm.itm2wgs84(689686, 266124)
    assert round(lat, 5) == 32.3
    assert round(lng, 5) == 35.7
