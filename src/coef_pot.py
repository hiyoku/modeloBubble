import numpy as np


def coef_pot(IMAX, JMAX, DY, DZ, CFO, DEN, AA, CC, SOURCE, E0, OMEGA1, EOZ, UY, UX, COSDIP, SENDIP, BO):
    EOX = E0
    AJUS = 1.0E3

    # AA, CC, SOURCE - Parametros inseridos
    DYLN_DEN = DZLN_DEN = np.zeros((IMAX, JMAX), dtype=np.float64)

    for I in range(1, JMAX - 1):
        BF = BO * np.power((6370.0 / ((I - 1) * DZ + 200.0 + 6370.0)), 3.0)
        GRAV = 9.81 * np.power((6370.0 / ((I - 1) * DZ + 200.0 + 6370.0)), 2.0)
        for J in range(1, IMAX - 2):
            AA[J, I] = (1.0 / DY) * (DEN[J + 1, I] - DEN[J - 1, I]) / (DEN[J + 1, I] + DEN[J - 1, I])
            CC[J, I] = (1.0 / DZ) * (CFO[I + 1] * DEN[J, I + 1] - CFO[I - 1] * DEN[J, I - 1]) / (CFO[I + 1] * DEN[J, I + 1] + CFO[I - 1] * DEN[J, I - 1])
            DYLN_DEN[J, I] = (1.0 / DY) * (DEN[J + 1, I] - DEN[J - 1, I]) / (DEN[J + 1, I] + DEN[J - 1, I])
            DZLN_DEN[J, I] = (1.0 / DZ) * (DEN[J, I + 1] - DEN[J, I - 1]) / (DEN[J, I + 1] + DEN[J, I - 1])
            DZCFO = (CFO[I + 1] - CFO[I - 1]) / (2.0 * DZ)

            SOURCE[J, I] = (EOX + COSDIP * GRAV * BF / CFO[I] + UX[I] * BF * SENDIP + BF * CFO[I] * UY[I] / OMEGA1) * DYLN_DEN[J, I] * AJUS + (EOZ[I] - GRAV * BF / OMEGA1 + COSDIP * UY[I] * BF) * CC[J, I] * AJUS + (EOZ[I + 1] - EOZ[I - 1]) * AJUS / (2.0 * DZ) + COSDIP * BF * (UY[I + 1] - UY[I - 1]) * AJUS / (2.0 * DZ)

    # Boundary Top - Bottom
    for K in range(0, IMAX):
        I1 = I2 = 0

        if K == 0:
            I1 = I2 = 1

        if K == IMAX - 1:
            I1 = I2 = -1

        print(DY, DEN[K + 1 + I2, JMAX - 1], DEN[K - 1 + I1, JMAX - 1], DEN[K + 1 + I2, JMAX - 1], DEN[K - 1 + I1, JMAX - 1])
        AA[K, JMAX - 1] = (1.0 / DY) * (DEN[K + 1 + I2, JMAX - 1] - DEN[K - 1 + I1, JMAX - 1]) / (DEN[K + 1 + I2, JMAX - 1] + DEN[K - 1 + I1, JMAX - 1])
        AA[K, 1] = (1.0 / DY) * (DEN[K + 1 + I2, 1] - DEN[K - 1 + I1, 1]) / (DEN[K + 1 + I2, 1] + DEN[K - 1 + I1, 1])
        CC[K, JMAX - 1] = (1.0 / DZ) * (CFO[JMAX - 1] - CFO[JMAX - 3]) / (CFO[JMAX - 1] + CFO[JMAX - 3])
        CC[K, 1] = (1.0 / DZ) * (CFO[3] - CFO[1]) / (CFO[3] + CFO[1])

        # Top for SOURCE
        DYLN_DEN[K, JMAX - 1] = AA[K, JMAX - 1]
        DZLN_DEN[K, JMAX - 1] = 0.0
        DZCFO = (CFO[JMAX - 1] - CFO[JMAX - 3]) / (2. * DZ)
        SOURCE[K, JMAX - 1] = (EOX + COSDIP * GRAV * BF / CFO[JMAX - 1]) * DYLN_DEN[K, JMAX - 1] * AJUS - (GRAV * BF / OMEGA1) * CC[K, JMAX - 1] * AJUS - 1

        # Bottom for SOURCE
        DYLN_DEN[K, 1] = AA[K, 1]
        DZLN_DEN[K, 1] = 0.0
        DZCFO = (CFO[3] - CFO[1]) / (2. * DZ)
        SOURCE[K, 1] = (EOX + COSDIP * GRAV * BF / CFO[1]) * DZLN_DEN[K, 1] * AJUS - (GRAV * BF / OMEGA1) * CC[K, 1] * AJUS

    # Boundary Left - Right
    for KX in range(1, JMAX - 1):
        AA[1, KX] = (1.0 / DY) * (DEN[3, KX] - DEN[1, KX]) / (DEN[3, KX] + DEN[1, KX])
        AA[IMAX - 1, KX] = (1.0 / DY) * (DEN[IMAX - 1, KX] - DEN[IMAX - 3, KX]) / (DEN[IMAX - 1, KX] + DEN[IMAX - 3, KX])
        CC[1, KX] = (1.0 / DZ) * (CFO[KX + 1] * DEN[1, KX + 1] - CFO[KX - 1] * DEN[1, KX - 1]) / (CFO[KX + 1] * DEN[1, KX + 1] + CFO[KX - 1] * DEN[1, KX - 1])
        CC[IMAX - 1, KX] = (1.0 / DZ) * (CFO[KX + 1] * DEN[IMAX - 1, KX + 1] - CFO[KX - 1] * DEN[IMAX - 1, KX - 1]) / (CFO[KX + 1] * DEN[IMAX - 1, KX + 1] + CFO[KX - 1] * DEN[IMAX - 1, KX - 1])

        # LEFT FOR SOURCE
        DYLN_DEN[1, KX] = (1.0 / DY) * (DEN[3, KX] - DEN[1, KX]) / (DEN[3, KX] + DEN[1, KX])
        DZLN_DEN[1, KX] = (1.0 / DY) * (DEN[1, KX + 1] - DEN[1, KX - 1]) / (DEN[1, KX + 1] + DEN[1, KX - 1])
        DZCFO = (CFO[KX + 1] - CFO[KX - 1]) / (2.0 * DZ)
        SOURCE[1, KX] = (EOX + COSDIP * GRAV * BF / CFO[KX]) * DYLN_DEN[1, KX] * AJUS - (GRAV * BF / OMEGA1) * CC[1, KX] * AJUS

        # RIGHT FOR SOURCE
        DYLN_DEN[IMAX - 1, KX] = (1.0 / DZ) * (DEN[IMAX - 1, KX] - DEN[IMAX - 3, KX]) / (DEN[IMAX - 1, KX] + DEN[IMAX - 3, KX])
        DZLN_DEN[IMAX - 1, KX] = (1.0 / DZ) * (DEN[IMAX - 1, KX + 1] - DEN[IMAX - 1, KX - 1]) / (DEN[IMAX - 1, KX + 1] + DEN[IMAX - 1, KX - 1])
        DZCFO = (CFO[KX + 1] - CFO[KX - 1]) / (2.0 * DZ)
        SOURCE[IMAX - 1, KX] = (EOX + COSDIP * GRAV * BF / CFO[KX]) * DYLN_DEN[IMAX - 1, KX] * AJUS - (GRAV * BF / OMEGA1) * CC[IMAX - 1, KX] * AJUS

    return AA, CC, SOURCE
