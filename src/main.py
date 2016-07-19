import leitores
# import potencial  # Função antiga de potencial
# from pot import potencial
from rotinas_fortran import (potencial, coef_pot, veloc_ion, density, line1356, line6300, tecn)
# import coef_pot
import numpy as np
import os
from datetime import datetime
from multiprocessing import Process


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
    t = Process(target=write_all, args=(DY, DZ, HB, NY, XLAT, IMAX, JMAX, DEN, VY, VZ, CO, CO2, CN2))
    t.start()


if __name__ == '__main__':
    # Pegando o tempo de inicio do programa.
    tempo_inicio = datetime.utcnow()

    # Definindo constantes
    NY = 81
    NZ = 121
    BO = 0.285E-04
    PI = np.pi

    F0F2 = leitores.leitor_frequencia_fof2(os.getcwd() + "/inputs/input-foF2.dat")
    EE0 = leitores.leitor_campo_eletrico_zonal("inputs/input-vz.dat")
    UX, UYZ = leitores.leitor_vento_termosferico("inputs/input-wind.dat")


    # Definindo parâmetros da malha
    IMAX = NY           # Contador en la horizontal
    JMAX = NZ           # Contador en la vertical
    DY = 5.0            # Paso en Longitud  [km]
    DZ = 5.0            # Paso en altura    [km]
    DTT = 10.0          # Paso en tiempo    [s]
    UT0 = 21.666        # Tiempo inicial en UT = 21:40

    # Parametros da região F
    HB = 200.0              # Altura base da região F [km]
    RK1 = 4.0E-11           # Recombinacion [O+] con [O2]
    RK2 = 1.3E-12           # Recombinacion [O+] con [N2]
    TN = 1200.4             # Temperatura neutra exosferica
    TE = 1800.0             # Temperatura electronica
    OMEGA1 = 150.19         # Girofrecuencia O+ [s]
    OMEGA2 = -4.396E06      # Girofrecuencia e- [s] q = -e
    ZMASAO = 2.66E-26       # Masa O+ [kg]
    ZMASAE = 9.11E-31       # Masa e- [kg]

    # Parametros da onda de perturbação zonal
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
        CO[J] = 8.557E+08 * np.exp(-(Z1 - 335.0) / HO)  # Oxigeno atomico    [cm-3]
        CO2[J] = 4.44E+06 * np.exp(-(Z1 - 335.0) / HO2)  # Oxigeno molecular [cm3]
        CN2[J] = 2.264E+08 * np.exp(-(Z1 - 335.0) / HN2)  # Nitrogeno molecular[cm3]
        BETA[J] = (RK1 * CO2[J]) + (RK2 * CN2[J])  # Recombinacion [1/s]

        # Colisões
        CFO[J] = 4.45E-11 * CO[J] * np.sqrt(TN) * (1.04 - 0.067 * np.log10(TN)) ** 2. + 6.64E-10 * CO2[J] + 6.82E-10 * CN2[J]

    # Perfil Vertical Electronico inicial
    DEN_I = np.zeros((JMAX))  # Criando vetor vazio
    for J in range(0, JMAX):
        DEN_I[J] = 1.24E04 * F0F2[J] ** 2.0

    # Parametros plano meridional magnetico
    for M in range(0, 21):
        XLAT = -30.0 + 3.0 * M
        PHI = XLAT * PI / 180.0
        COSLAT = np.cos(PHI)
        SENLAT = np.sin(PHI)
        COSDIP = COSLAT / np.sqrt(1.0 + 3.0 * SENLAT ** 2.0)
        SENDIP = 2. * SENLAT / np.sqrt(1. + 3. * SENLAT ** 2.0)
        F1 = COSLAT ** 3.0 / np.sqrt(1.0 + 3. * SENLAT ** 2.0)   # to  zonal electric
        F2 = COSLAT ** 3.0                                 # to vertical electric
        F3 = np.sqrt(1. + 3. * SENLAT ** 2.)
        FU = 1.0 / (0.7 + 0.4 * np.exp(-(XLAT) ** 2.0 / 200.0) - 0.035)

        # Malha inicial + perturbação zonal em t = 0
        DEN = np.zeros((IMAX, JMAX))
        DEN2 = np.zeros((IMAX, JMAX))
        T = np.zeros((IMAX, JMAX))
        AA = np.zeros((IMAX, JMAX))
        CC = np.zeros((IMAX, JMAX))
        SOURCE = np.zeros((IMAX, JMAX))  # Iniciando matrizes

        for J in range(0, JMAX):
            for I in range(0, IMAX):
                X1 = -DY * (NY - 1) / 2.0 + DY * I
                PERT = (1.0 - AMP * np.cos(ONDA * X1))
                DEN[I, J] = DEN_I[J] * PERT
                # print(DEN[I, J], DEN_I[J] * PERT)

        # print("PRINTANDO DEN")
        # print(DEN)
        # Programa Principal
        # Criando matrizes necessárias
        EOZ = np.zeros((JMAX))
        VY = np.zeros((IMAX, JMAX))
        VZ = np.zeros((IMAX, JMAX))
        DTDZ = np.zeros((IMAX, JMAX))
        DTDY = np.zeros((IMAX, JMAX))
        UY = np.zeros((JMAX))

        for K in range(0, 360):
            E0 = EE0[K]
            # for J in range(0, JMAX):
            #     EOZ[J] = -E0
            #     UY = UYZ
            #     UX[J] = 0.0
            EOZ = [(i - i - E0) for i in EOZ]
            UY = UYZ
            UX = [(i - i) for i in UX]

            DEN, AA, CC, SOURCE = coef_pot(IMAX, JMAX, DY, DZ, CFO, DEN, AA, CC, SOURCE, E0, OMEGA1, EOZ, UY, UX, COSDIP, SENDIP, BO)

            T = potencial(T, DY, DZ, IMAX, JMAX, AA, CC, SOURCE)

            VY, VZ, DTDZ, DTDY = veloc_ion(IMAX, JMAX, DY, DZ, T, CFO, OMEGA1, E0, UY, UX, EOZ, COSDIP, SENDIP, BO, VY, VZ, DTDZ, DTDY)

            DEN, DY, DZ, VY, VZ = density(DEN, IMAX, JMAX, DY, DZ, BETA, DTT, VY, VZ)

        # Normalizações
        start_norm(DY, DZ, HB, NY, XLAT, IMAX, JMAX, DEN, VY, VZ, CO, CO2, CN2)

    print("Tempo de execução -> " + str(datetime.utcnow() - tempo_inicio))
