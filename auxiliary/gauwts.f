C*********************************************************************
C subroutine GAUWTS **************************************************
C____________________________________________________________________C
C> Calculates points and weights for gauss-legendre quadrature
C>____________________________________________________________________
      SUBROUTINE GAUWTS ( X1, X2, GAUX, GAUW, NTHPTS )
      IMPLICIT NONE
C> (IN) Number of theta points
      INTEGER NTHPTS
C> (IN) Starting value for integration
      DOUBLE PRECISION X1
C> (IN) Ending value for integration
      DOUBLE PRECISION X2
C> (OUT) Array containing abscissae of the Gauss-Legendre 
C> NTHPTS-points quadrature formula.
      DOUBLE PRECISION GAUX( NTHPTS )
C> (OUT) Array containing the weights for the above points.
      DOUBLE PRECISION GAUW( NTHPTS )
C Variable declarations - Working variables 
      INTEGER M,I,J
      DOUBLE PRECISION XM,XL,P1,P2,P3,EPS,PP,Z,Z1
      PARAMETER (EPS=1.0d-13)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
      IF ( ABS( X1 - X2 ).LT.EPS ) THEN
        PRINT *,' Subroutine GAUWTS,'
        PRINT *,' X1 = ', X1
        PRINT *,' X2 = ', X2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C ....... the roots are symmetric in the interval so need only find
C half of them.
      M = ( NTHPTS + 1 )/2
      XM = 0.5d0 * ( X2 + X1 )
      XL = 0.5d0 * ( X2 - X1 )
C ........start looping over the desired roots .......................
      DO I = 1, M
         Z = DCOS( PI*( DBLE(I) - 0.25d0)/( DBLE(NTHPTS) + 0.5d0 ))
C           ..... starting with this approximation to the Ith root, we
C                enter the main loop of refinement by Newton's method.
 100     CONTINUE
            P1 = 1.0D0
            P2 = 0.0D0
C           ........... Loop up the recurrence relation to get the
C                      legendre Polynomial evaluated at Z.
            DO J = 1, NTHPTS
               P3 = P2
               P2 = P1
               P1 = ((2.0d0*J-1.0d0)*Z*P2 - (J-1.0d0)*P3)/DBLE( J )
            ENDDO
C           ..................... finish recurrence relation loop ...
C ... P1 is now the desired Legendre Polynomial. We now compute PP,
C    its derivative by a standard relation involving also P2, the 
C    polynomial of one order lower.
            PP = NTHPTS*(Z*P1-P2)/(Z*Z-1.0d0)
            Z1 = Z
            Z = Z1 - P1/PP
         IF ( ABS(Z-Z1).GT.EPS ) GOTO 100
C ...........scale the root to the desired interval .................
         GAUX( I ) = XM - XL*Z
C ...........and add its symmetric counterpart ......................
         GAUX( NTHPTS+1-I ) = XM + XL*Z
C ...........calculate the weight ...................................
         GAUW( I ) = 2.0d0*XL/((1.0d0-Z*Z)*PP*PP)
C ...........and add its symmetric counterpart ......................
         GAUW( NTHPTS + 1 - I ) = GAUW( I )
      ENDDO
C ......... end looping over the desired roots .......................
      RETURN
      END
C*********************************************************************

