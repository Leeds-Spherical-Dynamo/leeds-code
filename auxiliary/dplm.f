C*********************************************************************
C function DPLM *******************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C                                                                    
C> Calculates general P_l^m derivative from eqn 187 in my notes       
C>____________________________________________________________________
      FUNCTION DPLM ( L, M, X, S, PLMIN1, DPLMN1, DPLMN2 )
      IMPLICIT NONE
C> Derivative of Plm(x)
      DOUBLE PRECISION DPLM
C> Associated legendre order
      INTEGER M
C> Associated legendre degree
      INTEGER L
C> Cos(theta)
      DOUBLE PRECISION X
C> Sin(theta)
      DOUBLE PRECISION S
C> P_(l-1)^m ( X )  
      DOUBLE PRECISION PLMIN1
C> P_(l-1)^m ( X ) derivative 
      DOUBLE PRECISION DPLMN1
C> P_(l-2)^m ( X ) derivative   
      DOUBLE PRECISION DPLMN2
C Variable declarations - Working variables
      DOUBLE PRECISION RM,RL
C START OF FUNCTION
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function DPLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function DPLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' DPLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      DPLM = ( 2.0d0*RL - 1.0d0 )*(X*DPLMN1-S*PLMIN1)
      DPLM = DPLM - DPLMN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      DPLM = DPLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
