C   ================================
C     Potencial Electrico:  En MKS
C   ================================
      SUBROUTINE POTENCIAL(T,DY,DZ,IMAX,JMAX,AA,CC,SOURCE)
      PARAMETER (NY=81,NZ=121,PI=3.141592654)
      DIMENSION A(NY, NZ),   B(NY, NZ),  C(NY, NZ),
     *          D(NY, NZ),   F(NY, NZ),  SOURCE(NY, NZ),
     *          T(NY, NZ),  AA(NY, NZ),  CC(NY, NZ)

C   --------------------------------------------------------------
C   Despues de la discretizacion
C   T(i,j)= A*T(i+1,j) + B*T(i-1,j) + C*T(i,j+1) + D*T(i,j-1) - F
C   --------------------------------------------------------------
Cf2py intent(in) T, DY, DZ, IMAX, JMAX, AA, CC, SOURCE
Cf2py intent(out) T
C   *****************************
C    COEFICIENTES A, B, C, D, F
C
             DY2 = DY*DY
             DZ2 = DZ*DZ
              S1 = 2.0*(1.0/DY2 + 1.0/DZ2)

      DO 40 J=1,JMAX
      DO 39 I=1,IMAX
         A(I,J) = (1.0/S1)*(1.0/DY2 + AA(I,J)/(2.0*DY))
         B(I,J) = (1.0/S1)*(1.0/DY2 - AA(I,J)/(2.0*DY))
         C(I,J) = (1.0/S1)*(1.0/DZ2 + CC(I,J)/(2.0*DZ))
         D(I,J) = (1.0/S1)*(1.0/DZ2 - CC(I,J)/(2.0*DZ))
         F(I,J) = SOURCE(I,J)*(1.0/S1)
   39 CONTINUE
   40 CONTINUE
C -----------------------------
C         PARAMETROS DE SOR
C -----------------------------
C
      ITERA = 9000        !-- Contador de relajacion
      OMEGA = 1.7         !--  1 < OMEGA < 2
      EPS1  = 1.0E-05     !-- Tolerancia
C
C    -------------------------------------
C     Inicio de SOR
C     CALCULO: PUNTOS INTERNOS DE LA MALLA
C    -------------------------------------
      ANORMF = 0.0
      DO J=2,JMAX-1
      DO I=2,IMAX-1
      ANORMF = ANORMF + ABS(F(I,J))
      ENDDO
      ENDDO
C
      DO 300 K = 1,ITERA    !--- LOOP de SOR
      ANORM = 0.0
      CHGMAX = 0.0
      DO 110 J = 1,JMAX
        T(1,J) = T(IMAX,J)  !--- Condicion periodica
C           T(1,J) = T(2,J)/2.
C        T(IMAX,J) = T(1,J)
C           T(I,1) = T(I,2)/2
C        T(I,JMAX) = T(I,JMAX-1)/2.
      DO 105 I = 2,IMAX-1
        IF(J.EQ.1) GO TO 85
        IF(J.GT.1.AND.J.LT.JMAX) GO TO 86
        IF(J.EQ.JMAX) GO TO 87
C--------o
   85   BOTT = T(I,J+1)              !-- Bottom boundary dT/dZ=0
          TO = A(I,J)*T(I+1,J)+B(I,J)*T(I-1,J)+C(I,J)*T(I,J+1)
     *          +D(I,J)*BOTT-F(I,J)
      GO TO 99
C--------o
   86     TO =  A(I,J)*T(I+1,J)+B(I,J)*T(I-1,J)+C(I,J)*T(I,J+1)
     *          +D(I,J)*T(I,J-1)-F(I,J)
      GO TO 99
C--------o
   87    TOP = T(I,J-1)              !-- Top boundary  dT/dZ=0
          TO =  A(I,J)*T(I+1,J)+B(I,J)*T(I-1,J)+C(I,J)*TOP
     *          +D(I,J)*T(I,J-1)-F(I,J)
C--------o
          RESD = OMEGA*( TO-T(I,J) )
          ANORM = ANORM+ABS(RESD/OMEGA)
          IF(CHGMAX.LT.EPS1) CHGMAX = ABS(RESD)
   99     T(I,J) = T(I,J) + OMEGA*( TO-T(I,J) )
  105 CONTINUE
  110 CONTINUE
      IF(ANORM.LT.EPS1*ANORMF) GO TO 335
  300 CONTINUE                       !-- FIN DE ITERA
      WRITE(*,*)(' NOT CONVERGED!')
      STOP
  335 CONTINUE
      WRITE(*,550)OMEGA
  550 FORMAT('Omega:',F12.7)
      WRITE(*,600)K
  600 FORMAT('Interações:',I6)
      WRITE(*,610)EPS1
  610 FORMAT('Epsilon:', 1PE12.4)
      ERROR = EPS1*ANORM
      WRITE(*,620)ERROR
  620 FORMAT('Error:', 1PE12.4)
      WRITE(*,*)' '
      END
C
C