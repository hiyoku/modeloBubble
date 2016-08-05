import numpy as np
from scipy.io import netcdf

class Den():
    def __init__(self, xlat, den, vy, vz, dz, dy, co, co2, cn2, e0):
        self.xlat = xlat
        self.den = den
        self.vy = vy
        self.vz = vz
        self.dz = dz
        self.dy = dy
        self.co = co
        self.co2 = co2
        self.cn2 = cn2
        self.e0 = e0

    def get_info(self):
        return self.xlat, self.e0


# Function to search in a ListDen
def search(denList, xlat, e0):
    for den in denList:
        if den.xlat == xlat and den.e0 == np.float64(0.19E-04 * e0):
            return den

def create_netcdf(denList):
    ds = netcdf.netcdf_file('dataset_den.nc', 'w')
    ds.createDimension('dataListSize', len(denList))
    ds.createDimension('x', 81)
    ds.createDimension('y', 121)
    # ds.createDimension('xlat', 21)
    # ds.createDimension('e0', 540)

    arrayLat = ds.createVariable('xlat', 'f', ('dataListSize', ))
    arrayDen = ds.createVariable('den', 'f', ('dataListSize', 'x', 'y'))
    arrayE0 = ds.createVariable('e0', 'f', ('dataListSize', ))

    i = 0
    for den in denList:
        arrayLat[i] = den.xlat
        arrayE0[i] = den.e0
        arrayDen[i, ::] = den.den
        i += 1

    ds.close()

