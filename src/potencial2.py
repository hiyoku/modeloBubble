import numpy as np
from datetime import datetime


def potencial(t, den, dy, dz, imax, jmax, kn, aa, cc, source):
    dy2 = dy ** 2
    dz2 = dz ** 2
    s1 = 2.0 * (1.0 / dy2 + 1.0 / dz2)

    a = np.zeros((imax, jmax))
    b = np.zeros((imax, jmax))
    c = np.zeros((imax, jmax))
    d = np.zeros((imax, jmax))
    f = np.zeros((imax, jmax))

    inicio = datetime.utcnow()
    for j in range(0, jmax):
        for i in range(0, imax):
            a[i, j] = (1.0 / s1) * (1.0 / dy2 + aa[i, j] / (2.0 * dy))
            b[i, j] = (1.0 / s1) * (1.0 / dy2 - aa[i, j] / (2.0 * dy))
            c[i, j] = (1.0 / s1) * (1.0 / dz2 + cc[i, j] / (2.0 * dz))
            d[i, j] = (1.0 / s1) * (1.0 / dz2 - cc[i, j] / (2.0 * dz))
            f[i, j] = source[i, j] * (1.0 / s1)

    fim = datetime.utcnow()
    print(fim - inicio)
    # print(a, b, c, d, f)
    # SOR parameters
    itera = 9001    # Relaxation Counter
    omega = 1.7     # 1 < omega < 2
    eps1 = 1.0E-05   # threshold

    # SOR start

    # Calculus: GRID INTERNAL POINTS
    anormf = 0.0

    for ii in range(1, jmax - 1):
        for jj in range(1, imax - 1):
            anormf = anormf + abs(f[jj, ii])

    # SOR loop
    to = 0.0
    anorm = 1 + (eps1 * anormf)
    K = 0

    while (anorm > (eps1 * anormf)):
        inicio = datetime.utcnow()
        K = K + 1
    # for K in range(1, itera):
        anorm = 0.0
        for J in range(0, jmax):
            t[0, J] = t[imax - 1, J]    # periodic condition

        for kx in range(0, jmax):
            if kx == 0:
                aux = [1, 1]
            elif kx == jmax - 1:
                aux = [-1, -1]
            else:
                aux = [1, -1]

            for ky in range(1, imax - 1):
                to = a[ky, kx] * t[ky + 1, kx] + b[ky, kx] * t[ky - 1, kx] + c[ky, kx] * t[ky, kx + aux[0]] + d[ky, kx] * t[ky, kx + aux[1]] - f[ky, kx]

                if kx == jmax - 1:
                    resd = omega * (to - t[ky, kx])
                    anorm = anorm + abs(resd / omega)

                t[ky, kx] = t[ky, kx] + omega * (to - t[ky, kx])

        fim = datetime.utcnow()
        print(fim - inicio)

    print('Omega = ' + str(omega))
    print('Interações = ' + str(K))
    print('Epsilon = ' + str(eps1))
    error = eps1 * anorm
    print('Error = ' + str(error))
    print('')
    return t
