# !/usr/bin/python2
# Coding=UTF-8
import numpy as np
from scipy.io import netcdf


class BubbleModelStorage():
    def __init__(self, xlat, e0):
        # Indices
        self.xlat = xlat
        self.e0 = e0

        self.den = None
        self.vy = None
        self.vz = None
        self.dz = None
        self.dy = None
        self.co = None
        self.co2 = None
        self.cn2 = None

    def get_indices(self):
        """
            This function returns Lat and Time (e0)
        """
        return self.xlat, self.e0

    def set_info(self, den, vy, vz, dy, dz, co, co2, cn2):
        self.den = den
        self.vy = vy
        self.vz = vz
        self.dy = dy
        self.dz = dz
        self.co = co
        self.co2 = co2
        self.cn2 = cn2

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

