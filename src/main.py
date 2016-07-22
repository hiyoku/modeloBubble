# Making a compartible code
from __future__ import print_function

# Python Native lib imports
import numpy as np
import os
from datetime import datetime
from multiprocessing import Process

# Import PyGlow
from pyglow.pyglow import Point

# Program imports
import leitores
from rotinas_fortran import (potencial, coef_pot, veloc_ion, density, line1356, line6300, tecn)


def func_line1356(IMAX, JMAX, DZ, DEN, DY, XLAT):
    R1356 = np.zeros((IMAX))
    R1356 = line1356(IMAX, JMAX, DZ, DEN, R1356)

    return R1356


def func_line6300(IMAX, JMAX, DZ, DEN, CO, CO2, CN2):
    R6300 = np.zeros((IMAX))
    R6300 = line6300(IMAX, JMAX, DZ, DEN, R6300, CO, CO2, CN2)

    return R6300


def func_tecn(IMAX, JMAX, DZ, DEN):
    TEC = np.zeros((IMAX))
    TEC = tecn(IMAX, JMAX, DZ, DEN, TEC)

    return TEC


def write_pot_log(OMEGA, K, EPS1, ANORM):
    with open('logs/log-pot.txt', 'a') as log:
        log.write("OMEGA: {}\nInteracoes: {}\nEpsilon: {}\nError: {}\n".format(OMEGA, K, EPS1, ANORM * EPS1))
        log.write("============================================\n")


def write_all(DY, DZ, HB, NY, XLAT, IMAX, JMAX, DEN, VY, VZ, CO, CO2, CN2):
    R1356 = func_line1356(IMAX, JMAX, DZ, DEN, DY, XLAT)
    R6300 = func_line6300(IMAX, JMAX, DZ, DEN, CO, CO2, CN2)
    TEC = func_tecn(IMAX, JMAX, DZ, DEN)

    r1 = open('output/output-line1356.dat', 'a')
    r2 = open('output/output-line6300.dat', 'a')
    te = open('output/output-tec.dat', 'a')
    d1 = open('output/output-den.dat', 'a')

    for I in range(0, IMAX):
        X1 = - DY * (NY - 1) / 2.0 + DY * I

        r1.write("{}   {}   {}\n".format(X1, XLAT, np.log10(R1356[I] / 1E6)))
        r2.write("{}   {}   {}\n".format(X1, XLAT, np.log10(R6300[I] / 1E6)))
        te.write("{}   {}   {}\n".format(X1, XLAT, TEC[I] / 1.0E16))

        for J in range(0, JMAX):
            Z1 = HB + DZ * J
            d1.write("{}   {}   {}   {}   {}\n".format(X1, Z1, np.log10(DEN[I, J]), VY[I, J], VZ[I, J]))

    r1.close()
    r2.close()
    te.close()
    d1.close()


def start_norm(DY, DZ, HB, NY, XLAT, IMAX, JMAX, DEN, VY, VZ, CO, CO2, CN2):
    # t = Process(target=write_all, args=(DY, DZ, HB, NY, XLAT, IMAX, JMAX, DEN, VY, VZ, CO, CO2, CN2))
    # t.start()
    write_all(DY, DZ, HB, NY, XLAT, IMAX, JMAX, DEN, VY, VZ, CO, CO2, CN2)


def start_logging(OMEGA, K, EPS1, ANORM):
    l = Process(target=write_pot_log, args=(OMEGA, K, EPS1, ANORM))
    l.start()


def get_ne_alt(point):
    DEN_I = []
    i = 200

    while i <= 800:
        point.alt = i
        point.run_iri()
        DEN_I.append(point.ne)

        i += 5  # passo de 5km

    return DEN_I


def get_all_uy_ux(point, ver):
    UY = []
    UX = []
    i = 200.0

    while i <= 800.0:
        point.alt = i
        point.run_hwm(version=ver)
        UY.append(point.u)
        UX.append(point.v)

        i += 5  # passo em 5km

    return UY, UX


def loop_lat(M, DY, DZ, ponto, tempo_inicio):
    print('Processo {} iniciado!'.format(M))
    XLAT = -30.0 + 3.0 * M
    PHI = XLAT * PI / 180.0
    COSLAT = np.cos(PHI)
    SENLAT = np.sin(PHI)
    COSDIP = COSLAT / np.sqrt(1.0 + 3.0 * SENLAT ** 2.0)
    SENDIP = 2. * SENLAT / np.sqrt(1. + 3. * SENLAT ** 2.0)
    # F1 = COSLAT ** 3.0 / np.sqrt(1.0 + 3. * SENLAT ** 2.0)   # to  zonal electric
    # F2 = COSLAT ** 3.0                                 # to vertical electric
    # F3 = np.sqrt(1. + 3. * SENLAT ** 2.)
    # FU = 1.0 / (0.7 + 0.4 * np.exp(-(XLAT) ** 2.0 / 200.0) - 0.035)

    # Parametros do Point
    ponto.lat = XLAT
    UY, UX = get_all_uy_ux(ponto, versao_hwm)
    # np.savetxt('uy.txt', UY)
    # np.savetxt('ux.txt', UX)

    # Malha inicial + perturbacao zonal em t = 0
    DEN = np.zeros((IMAX, JMAX))
    T = np.zeros((IMAX, JMAX))
    AA = np.zeros((IMAX, JMAX))
    CC = np.zeros((IMAX, JMAX))
    SOURCE = np.zeros((IMAX, JMAX))  # Iniciando matrizes

    for J in range(0, JMAX):
        for I in range(0, IMAX):
            X1 = -DY * (NY - 1) / 2.0 + DY * I
            PERT = (1.0 - AMP * np.cos(ONDA * X1))
            DEN[I, J] = DEN_I[J] * PERT

    # Programa Principal
    # Criando matrizes necessarias
    EOZ = np.zeros((JMAX))
    VY = np.zeros((IMAX, JMAX))
    VZ = np.zeros((IMAX, JMAX))
    DTDZ = np.zeros((IMAX, JMAX))
    DTDY = np.zeros((IMAX, JMAX))
    UY = np.zeros((JMAX))
    OMEGA = 0
    K = 0
    EPS1 = 0
    ANORM = 0

    for E0 in EE0:  # Campo Electrico = Tempo = 540 interactions
        # for J in range(0, JMAX):
        #     EOZ[J] = -E0
        #     UY = UYZ
        #     UX[J] = 0.0
        EOZ = [(i - i - E0) for i in EOZ]
        # Tempo - Vento
        # UY = UYZ
        # UX = [(i - i) for i in UX]

        # Vetor Tempo - Vento (PONTO)

        DEN, AA, CC, SOURCE = coef_pot(IMAX, JMAX, DY, DZ, CFO, DEN, AA, CC, SOURCE, E0, OMEGA1, EOZ, UY, UX, COSDIP, SENDIP, BO)

        T, OMEGA, K, EPS1, ANORM = potencial(T, DY, DZ, IMAX, JMAX, AA, CC, SOURCE, OMEGA, K, EPS1, ANORM)

        start_logging(OMEGA, K, EPS1, ANORM)

        VY, VZ, DTDZ, DTDY = veloc_ion(IMAX, JMAX, DY, DZ, T, CFO, OMEGA1, E0, UY, UX, EOZ, COSDIP, SENDIP, BO, VY, VZ, DTDZ, DTDY)

        DEN, DY, DZ, VY, VZ = density(DEN, IMAX, JMAX, DY, DZ, BETA, DTT, VY, VZ)

    # Normalizacoes
    start_norm(DY, DZ, HB, NY, XLAT, IMAX, JMAX, DEN, VY, VZ, CO, CO2, CN2)

    print('Processo {} terminado em: {}'.format(M, str(datetime.utcnow() - tempo_inicio)))


if __name__ == '__main__':
    # Pegando o tempo de inicio do programa.
    tempo_inicio = datetime.utcnow()

    # Definindo constantes
    NY = 81
    NZ = 121
    BO = 0.285E-04
    PI = np.pi

    # Parametros PyGlow
    versao_hwm = 1993   # Ano da versao do HWM
    dn = datetime(2015, 10, 23, 21)
    lon = 0

    msis_alt = 335.0
    msis_lat = 0
    msis_lon = -40.0

    F0F2 = leitores.leitor_frequencia_fof2(os.getcwd() + "/inputs/input-foF2.dat")
    EE0 = leitores.leitor_campo_eletrico_zonal("inputs/input-vz.dat")
    UX, UYZ = leitores.leitor_vento_termosferico("inputs/input-wind.dat")

    # Criando o Point do PyGlow
    ponto = Point(dn, msis_lat, msis_lon, msis_alt)
    ponto.run_msis()

    # Definindo parametros da malha
    IMAX = NY           # Contador en la horizontal
    JMAX = NZ           # Contador en la vertical
    DY = 5.0            # Paso en Longitud  [km]
    DZ = 5.0            # Paso en altura    [km]
    DTT = 10.0          # Paso en tiempo    [s]
    UT0 = 21.666        # Tiempo inicial en UT = 21:40

    # Parametros da regiao F
    HB = 200.0              # Altura base da regiao F [km]
    RK1 = 4.0E-11           # Recombinacion [O+] con [O2]
    RK2 = 1.3E-12           # Recombinacion [O+] con [N2]
    # TN = 1200.4             # Temperatura neutra exosferica
    TN = np.float64(ponto.Tn_msis)
    TE = 1800.0             # Temperatura electronica
    OMEGA1 = 150.19         # Girofrecuencia O+ [s]
    OMEGA2 = -4.396E06      # Girofrecuencia e- [s] q = -e
    ZMASAO = 2.66E-26       # Masa O+ [kg]
    ZMASAE = 9.11E-31       # Masa e- [kg]

    # Parametros da onda de perturbacao zonal
    ZLAMDA = 400.0                  # Longitud de onda en km
    ONDA = 2.0 * PI / ZLAMDA        # Numero de onda
    AMP = 0.0470                    # Amplitud de la perturbacion

    # Atmosfera Neutra
    # Criando vetores vazios
    CO = np.zeros((JMAX))
    CO2 = np.zeros((JMAX))
    CN2 = np.zeros((JMAX))
    BETA = np.zeros((JMAX))
    CFO = np.zeros((JMAX))

    for J in range(0, JMAX):
        Z1 = HB + DZ * (J)
        GR = 1. / (1.0 + Z1 / 6370.0) ** 2.0
        HO = 0.0528 * TN / GR                           # Escala de altura O  [km]
        HO2 = 0.0264 * TN / GR                          # Escala de altura O2 [km]
        HN2 = 0.0302 * TN / GR                          # Escala de altura N2 [km]
        # CO[J] = 8.557E+08 * np.exp(-(Z1 - 335.0) / HO)  # Oxigeno atomico    [cm-3]
        # CO2[J] = 4.44E+06 * np.exp(-(Z1 - 335.0) / HO2)  # Oxigeno molecular [cm3]
        # CN2[J] = 2.264E+08 * np.exp(-(Z1 - 335.0) / HN2)  # Nitrogeno molecular[cm3]
        CO[J] = np.float64(ponto.nn['O']) * np.exp(-(Z1 - 335.0) / HO)  # Oxigeno atomico    [cm-3]
        CO2[J] = np.float64(ponto.nn['O2']) * np.exp(-(Z1 - 335.0) / HO2)  # Oxigeno molecular [cm3]
        CN2[J] = np.float64(ponto.nn['N2']) * np.exp(-(Z1 - 335.0) / HN2)  # Nitrogeno molecular[cm3]
        BETA[J] = (RK1 * CO2[J]) + (RK2 * CN2[J])  # Recombinacion [1/s]

        # Colisoes
        CFO[J] = 4.45E-11 * CO[J] * np.sqrt(TN) * (1.04 - 0.067 * np.log10(TN)) ** 2. + 6.64E-10 * CO2[J] + 6.82E-10 * CN2[J]

    # Perfil Vertical Electronico inicial
    # DEN_I = np.zeros((JMAX))  # Criando vetor vazio

    DEN_I = get_ne_alt(ponto)

    # for J in range(0, JMAX):
    #     DEN_I[J] = 1.24E04 * F0F2[J] ** 2.0

    # Parametros plano meridional magnetico
    for M in range(0, 21):
        Process(target=loop_lat, args=(M, DY, DZ, ponto, tempo_inicio)).start()

    # for p in prs:
    #     p.join()

    # print("Tempo de execucao -> " + str(datetime.utcnow() - tempo_inicio))
