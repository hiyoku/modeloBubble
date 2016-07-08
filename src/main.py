import leitores
from math import (pi as PI, exp as EXP,
    sqrt as SQRT, log10 as LOG10)
from datetime import datetime

# Pegando o tempo de inicio do programa.
tempo_inicio = datetime.utcnow()

# Definindo constantes
NY = 81
NZ = 121
BO = 0.285E-04

F0F2 = leitores.leitor_frequencia_fof2("inputs\\input-fof2.dat")
EE0 = leitores.leitor_campo_eletrico_zonal("inputs\\input-vz.dat")
UX, UYZ = leitores.leitor_vento_termosferico("inputs\\input-wind.dat")


# Definindo parâmetros da malha
IMAX = NY			# Contador en la horizontal
JMAX = NZ			# Contador en la vertical
DY = 5.0			# Paso en Longitud  [km]
DZ = 5.0			# Paso en altura    [km]
DTT = 10.0			# Paso en tiempo    [s]
UT0 = 21.666		# Tiempo inicial en UT = 21:40

# Parametros da região F
HB = 200.0				# Altura base da região F [km]
RK1 = 4.0E-11			# Recombinacion [O+] con [O2]
RK2 = 1.3E-12			# Recombinacion [O+] con [N2]
TN = 1035.66			# Temperatura neutra exosferica
TE = 1800.0				# Temperatura electronica
OMEGA1 = 150.19			# Girofrecuencia O+ [s]
OMEGA2 = -4.396E06		# Girofrecuencia e- [s] q = -e
ZMASAO = 2.66E-26		# Masa O+ [kg]
ZMASAE = 9.11E-31		# Masa e- [kg]

# Parametros da onda de perturbação zonal
ZLAMDA = 400.0					# Longitud de onda en km
ONDA = 2.0 * PI / ZLAMDA		# Numero de onda
AMP = 0.0470					# Amplitud de la perturbacion

# Atmosfera Neutra
for J in range(0, JMAX):
    print(J)
    # Criando vetores vazios
    CO = []
    CO2 = CN2 = BETA = CFO = [NZ]

    Z1 = HB + DZ * (J - 1)
    GR = 1. / (1.0 + Z1 / 6370.0) ** 2.
    HO = 0.0528 * TN / GR                           # Escala de altura O  [km]
    HO2 = 0.0264 * TN / GR                          # Escala de altura O2 [km]
    HN2 = 0.0302 * TN / GR                          # Escala de altura N2 [km]
    CO.append(5.29E+08 * EXP(-(Z1 - 335.0) / HO))   # Oxigeno atomico    [cm-3]
    CO2.append(2.84E+06 * EXP(-(Z1 - 335.0) / HO2)) # Oxigeno molecular  [cm-3]
    CN2.append(1.05E+08 * EXP(-(Z1 - 335.0) / HN2)) # Nitrogeno molecular[cm-3]
    BETA.append((RK1 * CO2[J]) + (RK2 * CN2[J]))    # Recombinacion       [1/s]

    # Colisões
    CFO.append(4.45E-11 * CO[J] * SQRT(TN) * (1.04 - 0.067 * LOG10(TN)) ** 2. + 6.64E-10 * CO2[J] + 6.82E-10 * CN2[J])


print("Tempo de execução -> " + str(datetime.utcnow() - tempo_inicio))
