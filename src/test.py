from pyglow.pyglow import Point
from datetime import datetime
from multiprocessing import Process
import numpy as np


def main():
    s1 = datetime.utcnow()
    dn = datetime(2011, 3, 23, 9, 30)
    lat = 0.
    lon = -80.
    alt = 250.00, 300.00

    pt = Point(dn, lat, lon, alt)
    pt2 = Point(dn, lat, lon, alt)
    pt3 = Point(dn, lat, lon, alt)

    pt.run_hwm(version=1993)
    pt2.run_hwm(version=2007)
    pt3.run_hwm(version=2014)

    print("v1993\nu: {}     v: {}\n".format(pt.u, pt.v))
    print("v2007\nu: {}     v: {}\n".format(pt2.u, pt2.v))
    print("v2014\nu: {}     v: {}\n".format(pt3.u, pt3.v))
    print('\n\n{}'.format(datetime.utcnow() - s1))

def main2():
    dn = datetime(2011, 3, 23, 9, 30)
    lat = 0.
    lon = -80.
    alt = 250.

    pt = Point(dn, lat, lon, alt)
    pt.run_igrf()
    pt.run_hwm93()
    pt.run_msis()
    pt.run_iri()

    print pt
    print pt.nn
    print pt.Tn_msis


def create_point(dn, lat, lon, alt, ver):
    pt = Point(dn, lat, lon, alt)

    pt.run_hwm(version=ver)

    print("v{}\nu: {}     v: {}\n".format(ver, pt.u, pt.v))


def run_parallel():
    s1 = datetime.utcnow()

    dn = datetime(2011, 3, 23, 9, 30)
    lat = 0.
    lon = -80.
    alt = 250.00, 300.00

    processes = []
    processes.append(Process(target=create_point, args=(dn, lat, lon, alt, 1993)))
    processes.append(Process(target=create_point, args=(dn, lat, lon, alt, 2007)))
    processes.append(Process(target=create_point, args=(dn, lat, lon, alt, 2014)))

    for p in processes:
        p.start()
    print('\n\n{}'.format(datetime.utcnow() - s1))


def get_ne_alt(point):
    DEN_I = []
    i = 200

    while i <= 800:
        point.alt = i
        point.run_iri()
        DEN_I.append(point.ne)

        i += 5  # passo de 5km

    return DEN_I

if __name__ == '__main__':
    # main()
    # main2()
    # run_parallel()

    JMAX = 121
    DEN_I = np.zeros((JMAX))

    dn = datetime(2015, 10, 23, 21)
    msis_alt = 335
    msis_lat = 0
    msis_lon = -40

    # Criando o Point do PyGlow
    ponto = Point(dn, msis_lat, msis_lon, msis_alt)
    DEN_I2 = get_ne_alt(ponto)

    print len(DEN_I)
    print len(DEN_I2)
