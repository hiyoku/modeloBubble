C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 20/07/2012
C     http://onlinelibrary.wiley.com/doi/10.1029/1999JA900136/full
C     http://iri.gsfc.nasa.gov/
C
C            1ra  version Nov.  2010
C            2da  version Jun.  2011
C            3ra  version Dic.  2012
C            4ta  version Abril 2014
C            5ta  version Marzo 2015
C            6ta  version Mayo  2016 Tesis Diego
C
C          INSTABILIDADE DE RAYLEIGH-TAYLOR
C                EQUADOR MAGNETICO
C              Altura-vs-Leste-Oeste
C
C               ALEXANDER CARRASCO, PhD
C               Departamento de Fisica
C              Universidad de Los Andes
C                  Merida-Venezuela
C       Post-Doutorado em DAE/INPE-BRASIL 2010-2011
C
C                    +z ^  Up
C                       |
C                       |
C                       |
C                       |
C             <---------o---------->
C           West      +x (out)     +y  East
C
C
C
C     ==================================================================
C                       ROTINAS FORTRAN PARA PYTHON
C                            UTILIZANDO F2PY
C                hideki.hiyoku@gmail.com (Marcos Hideki)
C     ==================================================================
C     ==================================================================

C     =============== INICIO DA ROTINA COEF_POT ========================
      SUBROUTINE COEF_POT(IMAX, JMAX, DY, DZ, CFO, DEN, AA, CC, SOURCE,
     *                    E0,OMEGA1,EOZ,UY,UX,COSDIP,SENDIP,BO)
C     IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NY=81, NZ=121, PI=3.141592654)
      DIMENSION  DEN(NY, NZ), AA(NY, NZ), CC(NY, NZ), SOURCE(NY, NZ),
     *           CFO(NZ), EOZ(NZ), DYLN_DEN(NY, NZ), UY(NZ), UX(NZ),
     *           DZLN_DEN(NY, NZ)

C     ===================== Parametros do F2PY =========================
Cf2py intent(in) IMAX, JMAX, DY, DZ, CFO, DEN, AA, CC, SOURCE, E0, OMEGA1, EOZ, UY, UX, COSDIP, SENDIP, BO
Cf2py intent(out) DEN, AA, CC, SOURCE

C     Inicio da rotina
      DY2 = DY * DY
      DZ2 = DZ * DZ
      EOX = E0          !-- Campo eletrico orden cero
      AJUS = 1.0E03     !-- Ajuste de unidades para MKS
      BF = 0.0
      GRAV = 0.0

C     =================== COEFICIENTES AA, CC, SOURCE
      DO J = 2, JMAX - 1
        BF = BO * (6370.0 / ((J-1)*DZ + 200.0 + 6370.0)) ** 3.
        GRAV = 9.81 * (6370.0/((J-1)*DZ + 200.0 + 6370.0)) ** 2.
        DO I = 2, IMAX - 1
          AA(I, J) = (1.0/DY)*(DEN(I+1, J) - DEN(I-1, J))
     *             / (DEN(I+1, J) + DEN(I-1, J))

          CC(I, J) = (1.0/DZ)*(CFO(J+1) * DEN(I, J+1)
     *               - CFO(J-1) * DEN(I, J-1)) / (CFO(J + 1)
     *               * DEN(I, J+1) + CFO(J-1)*DEN(I, J-1))

          DYLN_DEN(I,J) = (1./DY)*(DEN(I+1, J) - DEN(I-1, J))     !--DERIVADA ZONAL LN
     *                    /(DEN(I+1, J) + DEN(I-1, J))

          DZLN_DEN(I,J) = (1./DZ)*(DEN(I, J+1) - DEN(I, J-1))     !--DERIVADA VERT LN
     *                    /(DEN(I,J+1)+DEN(I,J-1))

          DZCFO = (CFO(J+1)-CFO(J-1))/(2.*DZ)              !--DERIVADA COLISION
C
          SOURCE(I,J) = (EOX + COSDIP*GRAV*BF/CFO(J) + UX(J)*BF*SENDIP
     *                  + BF*CFO(J)*UY(J)/OMEGA1  )*DYLN_DEN(I,J)*AJUS
     *                  + (EOZ(J)-GRAV*BF/OMEGA1 + COSDIP*UY(J)*BF)
     *                  * CC(I,J)*AJUS + (EOZ(J+1)-EOZ(J-1))
     *                  * AJUS/(2.0*DZ) + COSDIP*BF*(UY(J+1) - UY(J-1))
     *                  * AJUS/(2.0*DZ)
C    *                  - (GRAV*OMEGA1*BF/CFO(J)**2.0) * SENDIP**2.0
C    *                  * DZLN_DEN(I, J) * AJUS
        ENDDO
      ENDDO
C     ==================================================================
C                      Boundary TOP - BOTTOM
      DO I = 1, IMAX
        I1 = 0
        I2 = 0

        IF(I == 1) THEN
          I1 = 1
          I2 = 1
        ELSE IF(I == IMAX) THEN
          I1 = -1
          I2 = -1
        ENDIF

        AA(I, JMAX) = (1.0/DY) * (DEN(I + 1 + I2, JMAX)
     *                - DEN(I - 1 + I1, JMAX))
     *                / (DEN(I + 1 + I2, JMAX) + DEN(I - 1 + I1, JMAX))

        AA(I,1) = (1.0/DY) * (DEN(I + 1 + I2, 1) - DEN(I - 1 + I1, 1))
     *            /(DEN(I + 1 + I2, 1) + DEN(I - 1 + I1, 1))

        CC(I,JMAX) = (1.0 / DZ) * (CFO(JMAX) - CFO(JMAX - 2))
     *               /(CFO(JMAX) + CFO(JMAX - 2))

        CC(I,1) = (1.0 / DZ) * (CFO(3) - CFO(1)) / (CFO(3) + CFO(1))

C     =================== Top for SOURCE
        DYLN_DEN(I, JMAX) = AA(I, JMAX)
        DZLN_DEN(I, JMAX) = 0.0
        DZCFO = (CFO(JMAX) - CFO(JMAX-2)) / (2.0*DZ)

        SOURCE(I,JMAX) = (EOX + COSDIP * GRAV * BF / CFO(JMAX))
     *                   * DYLN_DEN(I, JMAX) * AJUS
     *                   - (GRAV*BF/OMEGA1) * CC(I, JMAX) * AJUS

C     =================== Bottom for SOURCE
        DYLN_DEN(I, 1) = AA(I, 1)
        DZLN_DEN(I, 1) = 0.0
        DZCFO = (CFO(3) - CFO(1)) / (2.*DZ)
        SOURCE(I, 1) = (EOX + COSDIP * GRAV * BF / CFO(1))
     *                 * DYLN_DEN(I, 1) * AJUS
     *                 - (GRAV * BF / OMEGA1) * CC(I, 1) * AJUS
      ENDDO

C     ==================================================================
C                      Boundary LEFT - RIGHT
      DO J = 2, JMAX - 1
        AA(1, J) = (1.0 / DY) * (DEN(3, J) - DEN(1, J)) / (DEN(3, J)
     *            + DEN(1, J))

        AA(IMAX, J) = (1.0 / DY) * (DEN(IMAX, J) - DEN(IMAX-2, J))
     *               / (DEN(IMAX, J) + DEN(IMAX-2, J))

        CC(1, J) = (1.0 / DZ) * (CFO(J+1) * DEN(1, J+1) - CFO(J-1)
     *             * DEN(1, J-1)) / (CFO(J+1) * DEN(1, J+1) + CFO(J-1)
     *             * DEN(1, J-1))

        CC(IMAX, J) = (1.0/DZ) * (CFO(J+1) * DEN(IMAX, J+1) - CFO(J-1)
     *                * DEN(IMAX, J-1)) / (CFO(J+1) * DEN(IMAX, J+1)
     *                + CFO(J-1) * DEN(IMAX, J-1))

C     =============== LEFT for SOURCE
        DYLN_DEN(1, J) = (1.0 / DY) * (DEN(3, J) - DEN(1, J))
     *                   / (DEN(3, J) + DEN(1, J))

        DZLN_DEN(1, J) = (1.0 / DZ) * (DEN(1, J+1) - DEN(1, J-1))
     *                  / (DEN(1, J+1) + DEN(1, J-1))

        DZCFO = (CFO(J+1) - CFO(J-1)) / (2.0 * DZ)

        SOURCE(1, J) = (EOX + COSDIP*GRAV*BF/CFO(J)) * DYLN_DEN(1, J)
     *                 * AJUS - (GRAV*BF/OMEGA1) * CC(1, J) * AJUS

C     =============== RIGHT for SOURCE
        DYLN_DEN(IMAX, J) = (1./DY) * (DEN(IMAX, J) - DEN(IMAX-2, J))
     *                      / (DEN(IMAX, J) + DEN(IMAX-2, J))

        DZLN_DEN(IMAX, J) = (1./DZ) * (DEN(IMAX, J+1) - DEN(IMAX, J-1))
     *                      /(DEN(IMAX,J+1)+DEN(IMAX,J-1))

        DZCFO = (CFO(J+1) - CFO(J-1)) / (2.0 * DZ)

        SOURCE(IMAX, J) = (EOX + COSDIP*GRAV*BF/CFO(J))
     *                    * DYLN_DEN(IMAX, J) * AJUS
     *                    - (GRAV*BF/OMEGA1) * CC(IMAX, J) * AJUS

      ENDDO
      END
C     ==================================================================
C     ===================== END COEF_POT FUNCTION ======================


C     ==================================================================
C     Função do pontencial eletrico
C     ==================================================================
      SUBROUTINE POTENCIAL(T,DY,DZ,IMAX,JMAX,AA,CC,SOURCE)
      PARAMETER (NY = 81, NZ = 121, PI = 3.141592654)
      DIMENSION A(NY, NZ),   B(NY, NZ),  C(NY, NZ),
     *          D(NY, NZ),   F(NY, NZ),  SOURCE(NY, NZ),
     *          T(NY, NZ),  AA(NY, NZ),  CC(NY, NZ)
      INTEGER AUX(2)
      REAL*4 TO

C     Despues de la discretizacion
C     T(i,j)= A*T(i+1,j) + B*T(i-1,j) + C*T(i,j+1) + D*T(i,j-1) - F

C     ===================== Parametros do F2PY =========================
Cf2py intent(in) T, DY, DZ, IMAX, JMAX, AA, CC, SOURCE
Cf2py intent(out) T
C     ==================================================================
      DY2 = DY * DY
      DZ2 = DZ * DZ
      S1 = 2.0 * (1.0/DY2 + 1.0/DZ2)

      DO J = 1, JMAX
        DO I = 1, IMAX
          A(I, J) = (1.0/S1) * (1.0/DY2 + AA(I, J) / (2.0*DY))
          B(I, J) = (1.0/S1) * (1.0/DY2 - AA(I, J) / (2.0*DY))
          C(I, J) = (1.0/S1) * (1.0/DZ2 + CC(I, J) / (2.0*DZ))
          D(I, J) = (1.0/S1) * (1.0/DZ2 - CC(I, J) / (2.0*DZ))
          F(I, J) = SOURCE(I, J) * (1.0/S1)
        ENDDO
      ENDDO

C     =================== Parametros do SOR ============================
      ITERA = 9000        !-- Contador de relajacion
      OMEGA = 1.7         !--  1 < OMEGA < 2
      EPS1 = 1.0E-05     !-- Tolerancia

C     ====================== Calculo do SOR ============================
C     CALCULO: PUNTOS INTERNOS DE LA MALLA
C     ==================================================================

      ANORMF = 0.0
      DO J = 2, JMAX - 1
        DO I = 2, IMAX - 1
          ANORMF = ANORMF + ABS(F(I, J))
        ENDDO
      ENDDO

C     ======================= Loop do SOR
      DO K = 1, ITERA
        ANORM = 0.0

        DO J = 1, JMAX
          T(1, J) = T(IMAX, J)  !--- Condicion periodica
C         T(1, J) = T(2, J) / 2.0
C         T(IMAX, J) = T(1, J)
C         T(I, 1) = T(I, 2) / 2
C         T(I, JMAX) = T(I, JMAX-1) / 2.0
          IF (J == 1) THEN
            AUX = [1, 1]
          ELSE IF(J == JMAX) THEN
            AUX = [-1, -1]
          ELSE
            AUX = [1, -1]
          ENDIF

          DO I = 2, IMAX-1
            TO = A(I, J) * T(I + 1, J) + B(I, J) * T(I - 1, J)
     *         + C(I, J) * T(I, J + AUX(1)) + D(I, J) * T(I, J + AUX(2))
     *         - F(I, J)

            IF (J == JMAX) THEN
              RESD = OMEGA*( TO-T(I,J) )
              ANORM = ANORM+ABS(RESD/OMEGA)
            ENDIF

            T(I,J) = T(I,J) + OMEGA*( TO-T(I,J) )
          ENDDO
        ENDDO

        IF(ANORM.LT.EPS1*ANORMF) THEN
          EXIT
        ENDIF
      ENDDO
C     ============================= End Loop
C       WRITE(*,550)OMEGA
C   550 FORMAT('Omega:',F12.7)

C       WRITE(*,600)K
C   600 FORMAT('Interações:',I6)

C       WRITE(*,610)EPS1
C   610 FORMAT('Epsilon:', 1PE12.4)

C       WRITE(*,620)EPS1*ANORM
C   620 FORMAT('Error:', 1PE12.4)
      END
C     ==================================================================
C     ===================== END POTENCIAL FUNCTION =====================


C     ==================================================================
C     Função de velocidade do Plasma O+
C     ==================================================================
      SUBROUTINE VELOC_ION(IMAX,JMAX,DY,DZ,T,CFO,OMEGA1,
     *              E0,UY,UX,EOZ,COSDIP,SENDIP,BO, VY, VZ, DTDZ, DTDY)
      PARAMETER (NY = 81, NZ = 121, PI = 3.141592654)
      DIMENSION T(NY,NZ),  DTDY(NY,NZ),  DTDZ(NY,NZ),     UX(NZ),
     *         VY(NY,NZ),    VZ(NY,NZ),      EOZ(NZ),     UY(NZ),
     *           CFO(NZ)

C     ===================== Parametros do F2PY =========================
Cf2py intent(in) IMAX,JMAX,DY,DZ,T,CFO,OMEGA1,E0,UY,UX,EOZ,COSDIP,SENDIP,BO, VY, VZ, DTDZ, DTDY
Cf2py intent(out) VY, VZ, DTDZ, DTDY
C     ==================================================================

      RT = 6370.0           !-- Radio de la Tierra en Km
      EOX = E0              !-- Campo electrico orden cero    [V/m]
      AJUS = 1.0E-03         !-- Ajuste de unidades para MKS

C     ===================== Pontos Internos
      DO J = 2,JMAX-1
        DO I = 2,IMAX-1
          DTDY(I, J) = (T(I+1, J) - T(I-1, J)) / (2.0*DY)
          DTDZ(I, J) = (T(I, J+1) - T(I, J-1)) / (2.0*DZ)
        ENDDO
      ENDDO

C     ===================== Boundary TOP - BOTTOM ======================
      DO I = 1,IMAX
        IF(I == 1) THEN
          I1 = 1
          I2 = 1
        ELSE IF(I == IMAX) THEN
          I1 = -1
          I2 = -1
        ELSE
          I1 = 0
          I2 = 0
        ENDIF

        DTDZ(I, 1) = 0.0
        DTDZ(I, JMAX) = 0.0
        DTDY(I, 1) = (T(I+1+I2, 1) - T(I-1+I1, 1)) / (2.0*DY)
        DTDY(I, JMAX) = (T(I+1+I2, JMAX) - T(I-1+I1, JMAX)) / (2.0*DY)
      ENDDO

C     ==================== Boundary LEFT - RIGHT =======================
      DO J = 2,JMAX-1
        DTDY(1, J) = (T(3, J) - T(1, J)) / (2.0*DY)
        DTDY(IMAX, J) = (T(IMAX, J) - T(IMAX-2, J)) / (2.0*DY)
        DTDZ(1, J) = (T(1, J+1) - T(1, J-1)) / (2.0*DZ)
        DTDZ(IMAX, J) = (T(IMAX, J+1) - T(IMAX, J-1)) / (2.0*DZ)
      ENDDO

C     ================ V = (Vy) i + (Vz) k   [m/s] =====================
C     Vy = COMPONENTE HORIZONTAL
C     Vz = COMPONENTE VERTICAL

      DO J = 1, JMAX
        Z1 = DZ*(J-1) + 200.0
        BF = BO !/(1.+Z1/RT)**3.
        GRAV = 9.81/(1.0 + Z1/RT)**2.
        DO I = 1, IMAX
          VY(I, J) =  CFO(J) * EOX / (BF*OMEGA1) - COSDIP * EOZ(J) / BF
     *             - (CFO(J)/OMEGA1*DTDY(I, J)/BF) * AJUS
     *             + (COSDIP*DTDZ(I, J)/BF) * AJUS
     *             + COSDIP*GRAV/OMEGA1

          VZ(I, J) = COSDIP*EOX/BF + COSDIP*UY(J)*CFO(J)/OMEGA1
     *             - ( COSDIP*DTDY(I,J)/BF )*AJUS
     *             - ( CFO(J)/OMEGA1*DTDZ(I,J)/BF )*AJUS
     *             + CFO(J)/OMEGA1*EOZ(J)/BF
     *             + UX(J)*COSDIP*SENDIP
C     *            -GRAV*SENDIP**2./CFO(J)
        ENDDO
      ENDDO
      END

C     ==================================================================
C     Função da Densidade de Plasma: Height vs East-West
C     ==================================================================
      SUBROUTINE DENSITY(DEN,IMAX,JMAX,DY,DZ,BETA,DTT,VY,VZ)
      PARAMETER (NY = 81, NZ = 121)
      DIMENSION    BETA(NZ),    DEN(NY, NZ),
     *       DDFO(NY, NZ),   DDGO(NY, NZ),
     *         VY(NY, NZ),     VZ(NY, NZ),   SOURCE(NY, NZ),
     * DYDZDEN_ES(NY, NZ), DEN_ES(NY, NZ),DYDZDEN_T(NY, NZ),
     *       DDSF(NY, NZ),   DDSG(NY, NZ),SOURCE_ES(NY, NZ),
     *      DEN_T(NY, NZ),    DDF(NY, NZ),      DDG(NY, NZ),
     * DYDZDEN_TD(NY, NZ), DEN_TD(NY, NZ),
     *       DDUF(NY, NZ),   DDUG(NY, NZ),
     *        ZFA(NY, NZ),    ZFB(NY, NZ),
C    *        G1C(NY, NZ),    G2C(NY, NZ),
C    *        F1C(NY, NZ),    F2C(NY, NZ),
C    *       DTDY(NY, NZ),   DTDZ(NY, NZ),
C    *          T(NY, NZ),
     *        ZGA(NY, NZ),    ZGB(NY, NZ)

C     ===================== Parametros do F2PY =========================
Cf2py intent(in) DEN, IMAX, JMAX, DY, DZ, BETA, DTT, VY, VZ
Cf2py intent(out) DEN, DY, DZ, VY, VZ
C     ==================================================================

C
C     ***************************************************************
C     Ecuacion de Continuidad
C
C        d[O+]    d [O+]Vy    d [O+]Vz
C       -----  +  --------- +  --------- = S(i,j)
C        dt         dY            dZ
C
C     El metodo Predictor-Corrector es usado para la integracion
C     en tiempo para obtener la solucion en un paso t*= to + dt/2.
C     Asi, la ecuacion discretizada queda en la forma:
C     ------------------------------------------------------------
C               *          o           o            o
C      dYdZ [O+] = dYdZ[O+](i,j)  - ( F(i+1/2,j) - F(i-1/2,j) )/2
C                                      o            o
C                                 - ( G(i,j+1/2) - G(i,j-1/2) )/2
C                                 + dYdZdt S(i,j)/2
C     -----------------------------------------------------------
C     La solucion anterior es usada para obtener [O+] en to + dt.
C     ----------------------------------------------------------
C               t          o           *            *
C      dYdZ [O+] = dYdZ[O+](i,j)  - ( F(i+1/2,j) - F(i-1/2,j) )
C                                      *            *
C                                 - ( G(i,j+1/2) - G(i,j-1/2) )
C                                 + dYdZdt S*(i,j)
C
C
C     Suma del termino DIFUSION FUERTE para reducir oscilaciones
C     inducidas por el metodo numerico (ver Boris-Book, 1976):
C      -----------------------------------------------------
C              td          t       d            d
C      dYdZ [O+]  = dYdZ[O+]  + ( F(i+1/2,j) - F(i-1/2,j) )
C                                  d            d
C                             + ( G(i,j+1/2) - G(i,j-1/2) )
C      -----------------------------------------------------
c
C                   td
C     La solucion [O+]   contiene un exceso de difusion numerica
C     que debe ser compensada para obtener una solucion adecuada.
C     Esto se resuelve con el uso de un termino de ANTI-DIFUSION
C     de flujo. Sin  embargo, la  solucion  obtenidad nuevamente,
C     contiene  un flujo  no  corregido. Asi que,  para corregir
C     estos flujos se debe  utilizar un corrector  de  flujo  de
C     transporte del tipo  Boris-Book. En  resumen, la  ecuacion
C     de continuidad es resolvida 4 veces.
C
C
C     ==================================================================
      DY = DY * 1.0E05
      DZ = DZ * 1.0E05
      DO I = 1, IMAX
        DO J = 1, JMAX
          VY(I, J) = VY(I, J) * 1.0E02
          VZ(I, J) = VZ(I, J) * 1.0E02
        ENDDO
      ENDDO

C     ==================================================================
C     ===================== STEP [1] - to + dt/2 =======================
C         Fo(I+1/2,J)  Y  Fo(I-1/2,J)
C     -------------------------------------  PUNTOS INTERNOS

      DO J = 2, JMAX - 1
        DO I = 2, IMAX - 1

          FO1 = (DEN(I+1, J) + DEN(I, J))
     *          * (VY(I+1, J) + VY(I, J)) * DZ * DTT / 4.0

          FO2 = (DEN(I, J) + DEN(I-1, J))
     *          * (VY(I, J) + VY(I-1, J)) * DZ * DTT / 4.0

          DDFO(I, J) = FO1 - FO2
        ENDDO
      ENDDO

C     -------------------------------------  PUNTOS INTERNOS
C         Go(I,J+1/2)  Y  Go(I,J-1/2)

      DO I = 2, IMAX - 1
        DO J = 2, JMAX - 1

          GO1 = (DEN(I, J+1) + DEN(I, J))
     *          * (VZ(I, J+1) + VZ(I, J)) * DY * DTT/4.0
          GO2 = (DEN(I, J) + DEN(I, J-1))
     *          * (VZ(I, J) + VZ(I, J-1)) * DY * DTT/4.0
          DDGO(I, J) = GO1 - GO2
        ENDDO
      ENDDO

C     ================= Boundary TOP - BOTTOM ==========================
      DO I = 1, IMAX
        IF(I == 1) THEN
          I1 = 1
          I2 = 1
        ELSE IF(I == IMAX) THEN
          I1 = -1
          I2 = -1
        ELSE
          I1 = 0
          I2 = 0
        ENDIF

        DDGO(I, JMAX) = DEN(I, JMAX) * (VZ(I, JMAX) - VZ(I, JMAX-1))/DZ
        DDGO(I, 1) = DEN(I, 1) * (VZ(I, 2) - VZ(I, 1)) / DZ

        DDFO(I, JMAX) = VY(I, JMAX)
     *                  * (DEN(I+1+I2, JMAX) - DEN(I-1+I1, JMAX))
     *                  /(2.0 * DY)
     *                  +DEN(I, JMAX)* (VY(I+1+I2,JMAX)-VY(I-1+I1,JMAX))
     *                 /(2.*DY)

        DDFO(I,1) = VY(I, 1) * (DEN(I+1+I2, 1) - DEN(I-1+I1, 1))/(2.*DY)
     *             + DEN(I, 1) * (VY(I+1+I2, 1) - VY(I-1+I1, 1))/(2.*DY)
      ENDDO

C     ================== Boundary LEFT - RIGHT =========================
      DO J = 2, JMAX-1
        DDFO(1, J) = (DEN(3, J)*VY(3, J) - DEN(1, J)*VY(1, J)) / (2.*DY)

        DDFO(IMAX, J) = (DEN(IMAX, J) * VY(IMAX, J) - DEN(IMAX-2, J)
     *                  * VY(IMAX-2, J)) / (2.0 * DY)

        DDGO(1, J) = (DEN(1, J+1) * VZ(1, J+1) - DEN(1, J-1)*VZ(1, J-1))
     *               / (2.0 * DZ)

        DDGO(IMAX, J) = (DEN(IMAX, J+1) * VZ(IMAX, J+1) - DEN(IMAX, J-1)
     *                  * VZ(IMAX, J-1)) / (2.0 * DZ)
      ENDDO

C     ============= TERMINO FUENTE: S(i,j)

      DO I = 1,IMAX
        DO J = 1,JMAX
          SOURCE(I, J) = - BETA(J) * DEN(I, J)
        ENDDO
      ENDDO

C     ============= SOLUCION  t* = to + dt/2

      DO I = 1, IMAX
        DO J = 1, JMAX
          DYDZDEN_ES(I,J) = DY * DZ * DEN(I, J)
     *                      - (DDFO(I, J) + DDGO(I, J)) / 2.0
     *                      + (1.0/2.0) * DY * DZ * DTT * SOURCE(I, J)

          DEN_ES(I, J) = DYDZDEN_ES(I, J) / (DY * DZ)
        ENDDO
      ENDDO

C     ==================================================================
C     ======================= STEP [2] - to + dt =======================
C         F(I+1/2,J)  Y  F(I-1/2,J)
C     ------------------------------------- PUNTOS INTERNOS

      DO J = 2, JMAX-1
        DO I = 2, IMAX-1
          FO1 = (DEN_ES(I+1, J) + DEN_ES(I, J)) * (VY(I+1, J)+VY(I, J))
     *          * DZ * DTT / 4.0

          FO2 = (DEN_ES(I, J) + DEN_ES(I-1, J)) * (VY(I, J) + VY(I-1,J))
     *          * DZ * DTT / 4.0

          DDSF(I, J) = FO1 - FO2
        ENDDO
      ENDDO

C         G(I,J+1/2)  Y  G(I,J-1/2)
C     -------------------------------------  PUNTOS INTERNOS
      DO I = 2, IMAX-1
        DO J = 2, JMAX-1

          GO1 = (DEN_ES(I, J+1) + DEN_ES(I, J)) * (VZ(I, J+1) + VZ(I,J))
     *         * DY * DTT / 4.0

          GO2 = (DEN_ES(I, J) + DEN_ES(I, J-1)) * (VZ(I, J) + VZ(I,J-1))
     *         * DY * DTT / 4.0

          DDSG(I, J) = GO1 - GO2
        ENDDO
      ENDDO

C     =================== Boundary TOP - BOTTOM ========================

      DO I = 1, IMAX
        IF(I == 1) THEN
          I1 = 1
          I2 = 1
        ELSE IF(I == IMAX) THEN
          I1 = -1
          I2 = -1
        ELSE
          I1 = 0
          I2 = 0
        ENDIF

        DDSG(I, JMAX) = DEN_ES(I, JMAX)
     *                  * (VZ(I, JMAX) - VZ(I, JMAX-1)) / DZ

        DDSG(I, 1) = DEN_ES(I, 1) * (VZ(I, 2) - VZ(I, 1)) / DZ
        DDSF(I,JMAX) = VY(I, JMAX) * (DEN_ES(I+1+I2, JMAX)
     *                 - DEN_ES(I-1+I1, JMAX)) / (2.0 * DY)
     *                 + DEN_ES(I, JMAX) * (VY(I+1+I2, JMAX)
     *                 - VY(I-1+I1, JMAX)) / (2.*DY)

        DDSF(I, 1) = VY(I, 1) * (DEN_ES(I+1+I2, 1) - DEN_ES(I-1+I1, 1))
     *              / (2.*DY) + DEN_ES(I, 1) * (VY(I+1+I2, 1)
     *              - VY(I-1+I1, 1)) / (2.0 * DY)
      ENDDO

C     ================== Boundary LEFT - RIGHT =========================
      DO J = 2, JMAX-1
        DDSF(1, J) = (DEN_ES(3, J) * VY(3, J) - DEN_ES(1, J) * VY(1, J))
     *              / (2.0 * DY)

        DDSF(IMAX, J) = (DEN_ES(IMAX, J) * VY(IMAX, J)
     *                - DEN_ES(IMAX-2, J)*VY(IMAX-2, J)) / (2.0 * DY)

        DDSG(1, J) = (DEN_ES(1, J+1) * VZ(1, J+1)
     *               - DEN_ES(1, J-1) * VZ(1, J-1)) / (2.0 * DZ)

        DDSG(IMAX, J) = (DEN_ES(IMAX, J+1) * VZ(IMAX, J+1)
     *                - DEN_ES(IMAX, J-1) * VZ(IMAX, J-1)) / (2.0 * DZ)
      ENDDO


C     ================= TERMINO FUENTE: S(i,j)*

      DO I = 1, IMAX
        DO J = 1, JMAX
          SOURCE_ES(I, J) = - BETA(J) * DEN_ES(I, J)
        ENDDO
      ENDDO

C     *****************************
C     SOLUCION t = to + dt

      DO I = 1,IMAX
        DO J = 1,JMAX
          DYDZDEN_T(I,J) = DY*DZ*DEN(I,J)
     *                  -( DDSF(I,J)+DDSG(I,J) )
     *                     + DY*DZ*DTT*SOURCE_ES(I,J)
          DEN_T(I,J) = DYDZDEN_T(I,J)/(DY*DZ)
        ENDDO
      ENDDO

C     ==================================================================
C     ======================= STEP [3] =================================
C                    SUMANDO EL FLUJO DISPERSIVO
C                     PARA REDUCIR OSCILACIONES
C                      EN LA SOLUCION NUMERICA

C         F(I+1/2,J)  Y  F(I-1/2,J)
C     ------------------------------------- PUNTOS INTERNOS
      DO J = 2, JMAX-1
        DO I = 2, IMAX-1
          EPSF1 = (1.0/2.0) * (VY(I+1, J) + VY(I, J)) * DTT / DY
          UF1 = 1.0/6.0 + (1.0/3.0) * (EPSF1**2.)
          EPSF2 = (1.0/2.0) * (VY(I, J) + VY(I-1, J)) * DTT / DY
          UF2 = 1.0/6.0 + (1.0/3.0) * (EPSF2**2.)
          FD1 = UF1 * DY * DZ * (DEN(I+1, J) - DEN(I, J))
          FD2 = UF2 * DY * DZ * (DEN(I, J) - DEN(I-1, J))
          DDF(I, J) = FD1 - FD2
        ENDDO
      ENDDO

C         G(I,J+1/2)  Y  G(I,J-1/2)
C   ------------------------------------- PUNTOS INTERNOS
      DO I = 2,IMAX-1
        DO J = 2,JMAX-1
          EPSG1 = (1.0/2.0) * (VZ(I, J+1) + VZ(I, J)) * DTT / DZ
          UG1 = 1.0/6.0 + (1.0/3.0) * (EPSG1**2.)
          EPSG2 = (1.0/2.0) * (VZ(I, J) + VZ(I, J-1)) * DTT / DZ
          UG2 = 1.0/6.0 + (1.0/3.0) * (EPSG2**2.)
          GD1 = UG1 * DY * DZ * (DEN(I, J+1) - DEN(I, J))
          GD2 = UG2 * DY * DZ * (DEN(I, J) - DEN(I, J-1))
          DDG(I, J) = GD1 - GD2
        ENDDO
      ENDDO

C     =================== Boundary TOP - BOTTOM ========================
      DO I = 1,IMAX
        IF(I == 1) THEN
          I1 = 0
          I2 = 1
        ELSE IF(I == IMAX) THEN
          I1 = -1
          I2 = 0
        ELSE
          I1 = 0
          I2 = 0
        ENDIF

        DDG(I,JMAX) = 0.0
        DDG(I,1) = 0.0

        EPSF1 = (1.0/2.0) * (VY(I+1+I1, JMAX) + VY(I+I1, JMAX)) * DTT/DY
        UF1 = 1.0/6.0 + (1.0/3.0) * (EPSF1**2.)

        EPSF2 = (1.0/2.0) * (VY(I+I2, JMAX) + VY(I-1+I2, JMAX)) * DTT/DY
        UF2 = 1.0/6.0 + (1.0/3.0) * (EPSF2**2.)

        FD1 = UF1 * DY * DZ * (DEN(I+1+I1, JMAX) - DEN(I+I1, JMAX))
        FD2 = UF2 * DY * DZ * (DEN(I+I2, JMAX) - DEN(I-1+I2, JMAX))

        DDF(I, JMAX) = FD1 - FD2

        EPSF1 = (1.0/2.0) * (VY(I+1+I1, 1) + VY(I+I1, 1)) * DTT / DY
        UF1 = 1.0/6.0 + (1.0/3.0) * (EPSF1**2.0)

        EPSF2 = (1.0/2.0) * (VY(I+I2, 1) + VY(I-1+I2, 1)) * DTT / DY
        UF2 = 1.0/6.0 + (1.0/3.0) * (EPSF2**2.)

        FD1 = UF1 * DY * DZ * (DEN(I+1+I1, 1) - DEN(I+I1, 1))
        FD2 = UF2 * DY * DZ * (DEN(I+I2, 1) - DEN(I-1+I2, 1))

        DDF(I, 1) = FD1 - FD2
      ENDDO

C     ================= Boundary LEFT - RIGHT ==========================
      DO J = 2, JMAX-1
        DDF(1, J) = 0.0
        DDF(IMAX, J) = 0.0

        EPSG1 = (1.0 / 2.0) * (VZ(1, J+1) + VZ(1, J)) * DTT/DZ
        UG1 = 1.0/6.0 + (1.0/3.0) * (EPSG1**2.)

        EPSG2 = (1.0/2.0) * (VZ(1, J) + VZ(1, J-1)) * DTT/DZ
        UG2 = 1.0/6.0 + (1.0/3.0) * (EPSG2**2.)

        GD1 = UG1 * DY * DZ * (DEN(1, J+1) - DEN(1, J))
        GD2 = UG2 * DY * DZ * (DEN(1, J) - DEN(1, J-1))

        DDG(1, J) = GD1 - GD2

        EPSG1 = (1.0/2.0) * (VZ(IMAX, J+1) + VZ(IMAX, J)) * DTT / DZ
        UG1 = 1.0/6.0 + (1.0/3.0) * (EPSG1**2.)

        EPSG2 = (1.0/2.0) * (VZ(IMAX, J) + VZ(IMAX, J-1)) * DTT / DZ
        UG2 = 1.0/6.0 + (1.0/3.0) * (EPSG2**2.)

        GD1 = UG1 * DY * DZ * (DEN(IMAX, J+1) - DEN(IMAX, J))
        GD2 = UG2 * DY * DZ * (DEN(IMAX, J) - DEN(IMAX, J-1))

        DDG(IMAX, J) = GD1 - GD2
      ENDDO

C     =============== SOLUCION CON FLUJO DISPERSIVO
      DO I = 1, IMAX
        DO J = 1, JMAX
          DYDZDEN_TD(I, J) = DYDZDEN_T(I, J) + (DDF(I, J) + DDG(I, J))
          DEN_TD(I, J) = DYDZDEN_TD(I, J) / (DY * DZ)
        ENDDO
      ENDDO
C     ==================================================================
C     ====================== STEP [4]
C       NOTA: La solucion anterior contiene un exceso  de  difusion
C         numerica que debe ser compensada  para que  se  pueda
C         obtener la solucion adecuada. El problema es resuelto
C         con la incorporacion  de  FLUJOS  ANTI-DIFUSIVOS, que
C         luego son adicionalmente corregidos.

C         F(I+1/2,J)  Y  F(I-1/2,J)
C   --------------------------------- PUNTOS INTERNOS
      DO J = 2, JMAX-1
        DO I = 2, IMAX-1
          EPSF1 = (1.0/2.0) * (VY(I+1, J) + VY(I, J)) * DTT/DY
          UBF1 = (1.0/6.0) * (1.0 - EPSF1**2.0)

          EPSF2 = (1.0/2.0) * (VY(I, J) + VY(I-1, J)) * DTT/DY
          UBF2 = (1.0/6.0) * (1.0 - EPSF2**2.0)

          FU1 = UBF1 * (DEN_T(I+1, J) - DEN_T(I, J)) * DY * DZ
          FU2 = UBF2 * (DEN_T(I, J) - DEN_T(I-1, J)) * DY * DZ

          ZFA(I, J) = FU1 / (DY * DZ)
          ZFB(I, J) = FU2 / (DY * DZ)

          DDUF(I, J) = FU1 - FU2

          EPSG1 = (1.0/2.0) * (VZ(I, J+1) + VZ(I, J)) * DTT/DZ
          UBG1 = (1.0/6.0) * (1.0 - EPSG1**2.0)

          EPSG2 = (1.0/2.0) * (VZ(I, J) + VZ(I, J-1)) * DTT/DZ
          UBG2 = (1.0/6.0) * (1.0 - EPSG2**2.0)

          GU1 = UBG1 * (DEN_T(I, J+1) - DEN_T(I, J)) * DY * DZ
          GU2 = UBG2 * (DEN_T(I, J) - DEN_T(I, J-1)) * DY * DZ

          ZGA(I, J) = GU1 / (DY * DZ)
          ZGB(I, J) = GU2 / (DY * DZ)

          DDUG(I, J) = GU1 - GU2
        ENDDO
      ENDDO

C     ================== Boundary TOP - BOTTOM
      DO I = 1, IMAX
        IF(I == 1) THEN
          I1 = 0
          I2 = 1
        ELSE IF(I == IMAX) THEN
          I1 = -1
          I2 = 0
        ELSE
          I1 = 0
          I2 = 0
        ENDIF

        DDUG(I, JMAX) = 0.0
        DDUG(I, 1) = 0.0

        EPSF1 = (1.0/2.0) * (VY(I+1+I1, JMAX) + VY(I+I1, JMAX)) * DTT/DY
        UBF1 = (1.0/6.0) * (1.0 - EPSF1**2.0)

        EPSF2 = (1.0/2.0) * (VY(I+I2, JMAX) + VY(I-1+I2, JMAX)) * DTT/DY
        UBF2 = (1.0/6.0) * (1.0 - EPSF2**2.0)

        FU1 = UBF1 * (DEN_T(I+1+I1, JMAX) - DEN_T(I+I1, JMAX)) * DY * DZ
        FU2 = UBF2 * (DEN_T(I+I2, JMAX) - DEN_T(I-1+I2, JMAX)) * DY * DZ
        DDUF(I, JMAX) = FU1 - FU2

        EPSF1 = (1.0/2.0) * (VY(I+1+I1, 1) + VY(I+I1, 1)) * DTT/DY
        UBF1 = (1.0/6.0) * (1.0 - EPSF1**2.0)

        EPSF2 = (1.0/2.0) * (VY(I+I2, 1) + VY(I-1+I2, 1)) * DTT/DY
        UBF2 = (1.0/6.0) * (1.0 - EPSF2**2.0)

        FU1 = UBF1 * (DEN_T(I+1+I1, 1) - DEN_T(I+I1, 1)) * DY * DZ
        FU2 = UBF2 * (DEN_T(I+I2, 1) - DEN_T(I-1+I2, 1)) * DY * DZ

        DDUF(I, 1) = FU1 - FU2
      ENDDO

C     ================== Boundary LEFT - RIGHT
      DO J = 2, JMAX-1
        DDUF(1, J) = 0.0
        DDUF(IMAX, J) = 0.0

        EPSG1 = (1.0/2.0) * (VZ(1, J+1) + VZ(1, J)) * DTT/DZ
        UBG1 = (1.0/6.0) * (1.0 - EPSG1**2.0)

        EPSG2 = (1.0/2.0) * (VZ(1, J) + VZ(1, J-1)) * DTT/DZ
        UBG2 = (1.0/6.0) * (1.0 - EPSG2**2.0)

        GU1 = UBG1 * (DEN_T(1, J+1) - DEN_T(1, J)) * DY * DZ
        GU2 = UBG2 * (DEN_T(1, J) - DEN_T(1, J-1)) * DY * DZ

        DDUG(1, J) = GU1 - GU2

        EPSG1 = (1.0/2.0) * (VZ(IMAX, J+1) + VZ(IMAX, J)) * DTT/DZ
        UBG1 = (1.0/6.0) * (1.0 - EPSG1**2.0)

        EPSG2 = (1.0/2.0) * (VZ(IMAX, J) + VZ(IMAX, J-1)) * DTT/DZ
        UBG2 = (1.0/6.0) * (1.0 - EPSG2**2.0)

        GU1 = UBG1 * (DEN_T(IMAX, J+1) - DEN_T(IMAX, J)) * DY * DZ
        GU2 = UBG2 * (DEN_T(IMAX, J) - DEN_T(IMAX, J-1)) * DY * DZ

        DDUG(IMAX, J) = GU1 - GU2
      ENDDO

C     NOTA: Los flujos [Fu] y [Gu]  deben ser corregidos  con el
C             limitador de Boris-Book si deseamos obtener soluciones
C             fisicas reales. El limitador es  un coeficiente  que
C             puede tomar valores entre 0 y 1.

C     =================== SOLUCION FINAL ===============================
C        LIMITADOR de Boris-Book (1976)

      ZLIM = 1.0/8.0

      DO I = 1, IMAX
        DO J = 1, JMAX
          DDUFC = ZLIM * DDUF(I, J)
          DDUGC = ZLIM * DDUG(I, J)

          DEN(I, J) = (DYDZDEN_TD(I, J) - DDUFC - DDUGC) / (DY*DZ)

          IF(DEN(I, J) .LT. 0.) DEN(I, J) = ABS(DEN(I, J))
        ENDDO

         DEN(I, 1) = DEN(I, 2) / 3.0
      ENDDO

      DO J=1, JMAX
        DEN(IMAX, J) = DEN(1, J)
      ENDDO

C     Retornando  VY, VZ, a las unidades MKS
      DY = DY / 1.0E05
      DZ = DZ / 1.0E05

      DO I=1, IMAX
        DO J=1, JMAX
          VY(I, J) = VY(I, J) / 1.0E02
          VZ(I, J) = VZ(I, J) / 1.0E02
        ENDDO
      ENDDO

      END

C     ==================================================================
C     Integracion metodo de Simpson 1/3  Line 1356

      SUBROUTINE LINE1356(IMAX, JMAX, DZ, DEN, R1356)
      PARAMETER(NY = 81, NZ = 121)
      DIMENSION DEN(NY, NZ),  R1356(NY)

C     ===================== Parametros do F2PY =========================
Cf2py intent(in) IMAX, JMAX, DZ, DEN, R1356
Cf2py intent(out) R1356
C     ==================================================================
C     Linea espectral 1356: Modelo de Melendez-Alvira, 1999 JGR

      ALFA = 7.3E-13  !-- cm^3 s^-1

C     =============== Integração
      DO I=1, IMAX
        SUMA = 0.0
        DO J=1, JMAX-52   !--52 MEJOR
          SUMA = SUMA + (DZ/3.0) * (DEN(I, J)**2.0
     *           + 4.0 * DEN(I, J+1)**2. + DEN(I, J+2)**2.)
        ENDDO

      R1356(I) = SUMA * ALFA * 1.0E5
      ENDDO
      END

C     ==================================================================
C     Integracion metodo de Simpson 1/3  Line 6300 [O+] = [e]

      SUBROUTINE LINE6300(IMAX, JMAX, DZ, DEN, R6300, CO, CO2, CN2)
      PARAMETER(NY = 81, NZ = 121)
      DIMENSION DEN(NY, NZ), V630(NY, NZ), R6300(NY),
     *               CO(NZ),      CO2(NZ),   CN2(NZ)

C     ===================== Parametros do F2PY =========================
Cf2py intent(in) IMAX, JMAX, DZ, DEN, R6300, CO, CO2, CN2
Cf2py intent(out) R6300
C     ==================================================================
      ZK2 = 2.30E-11
      ZK3 = 1.06E-11
      ZK5 = 3.20E-11
      ZK6 = 6.60E-10
      ZK7 = 9.20E-13
      A1D = 7.45E-03
      F1D = 1.1

      DO I=1, IMAX
        DO J=1, JMAX
          AX1 = 0.756 * F1D * ZK3 * CO2(J) * DEN(I, J)
          AX2 = 1.0 + (ZK2*CN2(J) + ZK5*CO2(J) + ZK6*DEN(I, J)
     *          + ZK7 * CO(J)) / A1D
          V630(I, J) = AX1 / AX2
        ENDDO
      ENDDO

      DO I=1,IMAX
        SUMA = 0.0
        DO J=1,JMAX-52
          SUMA = SUMA + (DZ/3.)*(V630(I,J)+4.*V630(I,J+1)+V630(I,J+2))
        ENDDO
        R6300(I) = SUMA*1.0E05
      ENDDO
      END

C     ==================================================================
C                     TEC Contenido electronico

      SUBROUTINE TECN(IMAX, JMAX, DZ, DEN, TEC)
      PARAMETER(NY = 81, NZ = 121)
      DIMENSION DEN(NY, NZ), TEC(NY)
C     ===================== Parametros do F2PY =========================
Cf2py intent(in) IMAX, JMAX, DZ, DEN, TEC
Cf2py intent(out) TEC
C     ==================================================================

C     TEC franja 250-500 km    J=11,61
      DO I=1,IMAX
        SUMA = 0.0

        DO J=1,JMAX-52
          SUMA = SUMA + (DZ/3.)*(DEN(I, J) + 4.*DEN(I, J+1)+DEN(I, J+2))
        ENDDO

        TEC(I) = SUMA*1.0E5*1.0E04
      ENDDO
      END