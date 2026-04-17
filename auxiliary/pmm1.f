C*********************************************************************
C function PMM1 ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C> Evaluates the Schmidt Normalised Legendre Function P_(m+1)^m (X)   
C> according to equation 179 in my notes ; i.e.                       
C>                                                                    
C> P_(m+1)^m (X) = SQRT( 2m+1 ).X.P^m_m(X)                         
C>                                                                    
C>____________________________________________________________________
      FUNCTION PMM1 ( M, X, PMM0 )
      IMPLICIT NONE
C> P_(m+1)^m (X) = SQRT( 2m+1 ).X.P^m_m(X)
      DOUBLE PRECISION PMM1
C> Polynomial Order
      INTEGER M
C> Cos(theta)
      DOUBLE PRECISION X
C> P_m^m(X) as evaluated by Function PMM. 
      DOUBLE PRECISION PMM0
C Variable declarations - Working variables
      DOUBLE PRECISION RM
C START OF FUNCTION 
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      PMM1 = X*PMM0*DSQRT( 2.0d0*RM+1.0d0 )
      RETURN
      END
C*********************************************************************
