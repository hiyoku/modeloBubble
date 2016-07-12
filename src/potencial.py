import numpy as np
import os


def potencial(t, den, dy, dz, hb, imax, jmax, kn, aa, cc, source):
    try:
        dy2 = np.float64(np.power(dy, 2))
        dz2 = np.float64(np.power(dz, 2))
        s1 = np.float64(2.0 * (1.0 / dy2 + 1.0 / dz2))

        a = b = c = d = f = np.zeros((imax, jmax), dtype=np.float64)

        for j in range(0, jmax):
            for i in range(0, imax):
                a[i, j] = np.float64((1.0 / s1) * (1.0 / dy2 + aa[i, j] / (2.0 * dy)))
                b[i, j] = np.float64((1.0 / s1) * (1.0 / dy2 - aa[i, j] / (2.0 * dy)))
                c[i, j] = np.float64((1.0 / s1) * (1.0 / dz2 + cc[i, j] / (2.0 * dz)))
                d[i, j] = np.float64((1.0 / s1) * (1.0 / dz2 - cc[i, j] / (2.0 * dz)))
                f[i, j] = np.float64(source[i, j] * (1.0 / s1))

        print(a, b, c, d, f)
        # SOR parameters
        itera = 9001    # Relaxation Counter
        omega = np.float64(1.7)     # 1 < omega < 2
        eps1 = np.float64(1.0E-05)   # threshold

        # SOR start

        # Calculus: GRID INTERNAL POINTS
        anormf = np.float64(0.0)

        for ii in range(1, jmax - 1):
            for jj in range(1, imax - 1):
                anormf = np.float64(anormf + np.abs(f[jj, ii]))

        # SOR loop
        to = np.float64(0.0)

        for K in range(1, itera):
            anorm = np.float64(0.0)
            # print("Interação -> {}".format(K))
            for J in range(0, jmax):
                t[0, J] = t[imax - 1, J]    # periodic condition

            for kx in range(0, jmax):
                if kx == 0:
                    aux = [1, 1]
                if kx > 0 and kx < jmax - 1:
                    aux = [1, -1]
                if kx == jmax - 1:
                    aux = [-1, -1]

                for ky in range(1, imax - 1):
                    a1 = a[ky, kx] * t[ky + 1, kx]
                    print(a1)
                    a2 = b[ky, kx] * t[ky - 1, kx]
                    print(a2)
                    a3 = c[ky, kx] * t[ky, kx + aux[0]]
                    print(a3)
                    a4 = d[ky, kx] * t[ky, kx + aux[1]]
                    print(a4)
                    to = a1 + a2 + a3 + a4 - f[ky, kx]
                    print("to: {}".format(to))

                    if kx == jmax - 1:
                        resd = omega * (to - t[ky, kx])
                        anorm = anorm + np.abs(resd / omega)

                    t[ky, kx] = t[ky, kx] + omega * (to - t[ky, kx])

            print("anorm: {} eps1: {} anormf: {}".format(anorm, eps1, anormf))
            if anorm < (eps1 * anormf):
                print('Omega = ' + str(omega))
                print('Interações = ' + str(K))
                print('Epsilon = ' + str(eps1))
                error = eps1 * anorm
                print('Error = ' + str(error))
                print('')
                return t
    except Exception as e:
        print(e)
