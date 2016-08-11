# !/usr/bin/python2
# -*- coding: utf-8 -*-
import numpy as np
from scipy.io import netcdf


class BubbleModelStorage():
    def __init__(self, xlat, e0):
        # Indices
        self.xlat = xlat
        self.e0 = e0

        # Iniciando todas as variáveis
        # 2D Array
        self.den = None
        self.vy = None
        self.vz = None
        self.t = None

        # 1D Array
        self.co = None
        self.co2 = None
        self.cn2 = None
        self.beta = None
        self.cfo = None
        self.den_i = None
        self.ux = None
        self.uy = None
        self.eoz = None

        # Constants
        self.dy = None
        self.dz = None
        self.bo = None
        self.jmax = None
        self.imax = None
        self.dtt = None
        self.ut0 = None
        self.hb = None
        self.rk1 = None
        self.rk2 = None
        self.tn = None
        self.te = None
        self.omega1 = None
        self.zlamda = None
        self.onda = None
        self.amp = None
        self.pyglow_o = None
        self.pyglow_o2 = None
        self.pyglow_n2 = None

    def get_indices(self):
        """
            This function returns Lat and Time (e0)
        """
        return self.xlat, self.e0

    def set_2d_array(self, den, vy, vz, t):
        self.den = den
        self.vy = vy
        self.vz = vz
        self.t = t

    def set_1d_array(self, co, co2, cn2, beta, cfo, den_i, ux, uy, eoz):
        self.co = co
        self.co2 = co2
        self.cn2 = cn2
        self.beta = beta
        self.cfo = cfo
        self.den_i = den_i
        self.ux = ux
        self.uy = uy
        self.eoz = eoz

    def set_escalar(self, dy, dz, bo, jmax, imax, dtt, ut0, hb, rk1, rk2, tn, te, omega1, zlamda, onda, amp, o, o2, n2):
        self.dy = dy
        self.dz = dz
        self.bo = bo
        self.jmax = jmax
        self.imax = imax
        self.dtt = dtt
        self.ut0 = ut0
        self.hb = hb
        self.rk1 = rk1
        self.rk2 = rk2
        self.tn = tn
        self.te = te
        self.omega1 = omega1
        self.zlamda = zlamda
        self.onda = onda
        self.amp = amp
        self.pyglow_o = o
        self.pyglow_o2 = o2
        self.pyglow_n2 = n2


# Static Functionctions
# Function to search in a List of BMS
def find_one(denList, xlat, e0):
    for den in denList:
        if den.xlat == xlat and den.e0 == np.float64(0.19E-04 * e0):
            return den


# Find a Storage with the same latitude
def find_by_lat(denList, xlat):
    resultList = []
    for den in denList:
        if den.xlat == xlat:
            resultList.append(den)

    return resultList

def find_by_tempo(denList, e0):
    resultList = []
    for den in denList:
        if den.e0 == np.float64(0.19E-04 * e0):
            resultList.append(den)
    return resultList


# Create a NETCDF file, Name is a optional input
def create_netcdf(denList, name=None):
    test = denList[1]

    if (name is None):
        ds = netcdf.netcdf_file('dataset_modelo_bubble_v6.nc', 'w')
    else:
        ds = netcdf.netcdf_file(str(name + '.nc'), 'w')
    ds.createDimension('dataListSize', len(denList))
    ds.createDimension('x', test.imax)
    ds.createDimension('y', test.jmax)
    ds.createDimension('constants', 1)

    arrayLat = ds.createVariable('xlat', 'f', ('dataListSize', ))
    arrayE0 = ds.createVariable('e0', 'f', ('dataListSize', ))

    arrayDen = ds.createVariable('den', 'f', ('dataListSize', 'x', 'y'))
    arrayVy = ds.createVariable('vy', 'f', ('dataListSize', 'x', 'y'))
    arrayVz = ds.createVariable('vz', 'f', ('dataListSize', 'x', 'y'))
    arrayT = ds.createVariable('t', 'f', ('dataListSize', 'x', 'y'))

    arrayCo = ds.createVariable('co', 'f', ('y', ))
    arrayCo2 = ds.createVariable('co2', 'f', ('y', ))
    arrayCn2 = ds.createVariable('cn2', 'f', ('y', ))
    arrayBeta = ds.createVariable('beta', 'f', ('y', ))
    arrayCfo = ds.createVariable('cfo', 'f', ('y', ))
    arrayDenI = ds.createVariable('den_i', 'f', ('y', ))
    ux = ds.createVariable('ux', 'f', ('y', ))
    uy = ds.createVariable('uy', 'f', ('y', ))
    eoz = ds.createVariable('eoz', 'f', ('y', ))

    Dy = ds.createVariable('dy', 'f', ('constants',))
    Dz = ds.createVariable('dz', 'f', ('constants',))
    Bo = ds.createVariable('bo', 'f', ('constants',))
    jmax = ds.createVariable('jmax', 'f', ('constants',))
    imax = ds.createVariable('imax', 'f', ('constants',))
    dtt = ds.createVariable('dtt', 'f', ('constants',))
    ut0 = ds.createVariable('ut0', 'f', ('constants',))
    hb = ds.createVariable('hb', 'f', ('constants',))
    rk1 = ds.createVariable('rk1', 'f', ('constants',))
    rk2 = ds.createVariable('rk2', 'f', ('constants',))
    tn = ds.createVariable('tn', 'f', ('constants',))
    te = ds.createVariable('te', 'f', ('constants',))
    omega1 = ds.createVariable('omega1', 'f', ('constants',))
    zlamda = ds.createVariable('zlamda', 'f', ('constants',))
    onda = ds.createVariable('onda', 'f', ('constants',))
    amp = ds.createVariable('amp', 'f', ('constants',))
    pyglow_o = ds.createVariable('pyglow_o', 'f', ('constants',))
    pyglow_o2 = ds.createVariable('pyglow_o2', 'f', ('constants',))
    pyglow_n2 = ds.createVariable('pyglow_n2', 'f', ('constants',))

    i = 0
    for den in denList:
        arrayLat[i] = den.xlat
        arrayE0[i] = den.e0
        arrayDen[i, ::] = den.den
        arrayVy[i, ::] = den.vy
        arrayVz[i, ::] = den.vz
        arrayT[i, ::] = den.t

        i += 1

    arrayCo[:] = den.co
    arrayCo2[:] = den.co2
    arrayCn2[:] = den.cn2
    arrayBeta[:] = den.beta
    arrayCfo[:] = den.cfo
    arrayDenI[:] = den.den_i
    ux[:] = den.ux
    uy[:] = den.uy
    eoz[:] = den.eoz

    Dy[:] = den.dy
    Dz[:] = den.dz
    Bo[:] = den.bo
    jmax[:] = den.jmax
    imax[:] = den.imax
    dtt[:] = den.dtt
    ut0[:] = den.ut0
    hb[:] = den.hb
    rk1[:] = den.rk1
    rk2[:] = den.rk2
    tn[:] = den.tn
    te[:] = den.te
    omega1[:] = den.omega1
    zlamda[:] = den.zlamda
    onda[:] = den.onda
    amp[:] = den.amp
    pyglow_o[:] = den.pyglow_o
    pyglow_o2[:] = den.pyglow_o2
    pyglow_n2[:] = den.pyglow_n2

    ds.close()
