C*********************************************************************
C function DPMM1 *****************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C> Calculates d(P_(m+1)^m/d(theta) by equation 185 in my notes        
C>____________________________________________________________________
      FUNCTION DPMM1 ( M, X, S, PMM0, DPMM0)
      IMPLICIT NONE
C> d(P_(m+1)^m/d(theta)
      DOUBLE PRECISION DPMM1
C> Polynomial order
      INTEGER M
C> Cos(theta)
      DOUBLE PRECISION X
C> Sin(theta)
      DOUBLE PRECISION S
C> P_m^m(X) as evaluated by Function PMM. 
      DOUBLE PRECISION PMM0
C> d(P_m^m(X))/d(theta) as evaluated by Function DPMM
      DOUBLE PRECISION DPMM0
C Variable declarations - Working variables
      DOUBLE PRECISION RM
C START OF FUNCTION
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      DPMM1 = DSQRT( 2.0d0*RM+1.0d0 )*( X*DPMM0 - S*PMM0 )
      RETURN
      END
C*********************************************************************
