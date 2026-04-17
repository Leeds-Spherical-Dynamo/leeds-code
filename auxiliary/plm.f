C*********************************************************************
C function PLM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C> Calculates the Schmidt Normalised Legendre Function P_l^m (x)      
C> given P_(l-1)^m and P_(l-2)^m according to equation 183 in my notes
C> i.e.                                                               
C>   P_l^m( X ) = { A * P_(l-1)^m - B * P_(l-2)^m }/C                 
C>                                                                    
C> where A = (2*l - 1)*X ,                                            
C>                                                                    
C> B = SQRT( (L+M-1)*(L-M-1) ) and C = SQRT( (L+M)*(L-M) )
C>__________________________________________________________________
      FUNCTION PLM ( L, M, X, PLMIN1, PLMIN2 )
      IMPLICIT NONE
C> Plm(x)
      DOUBLE PRECISION PLM
C> Order of Associated legendre polynomial
      INTEGER M
C> Degree of Associated legendre polynomial
      INTEGER L
C> Cos(theta)
      DOUBLE PRECISION X 
C> P_(l-1)^m ( X )
      DOUBLE PRECISION PLMIN1
C> P_(l-2)^m ( X ) 
      DOUBLE PRECISION PLMIN2
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C START OF FUNCTION *************************************************C

      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function PLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function PLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' PLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      PLM = ( 2.0d0*RL - 1.0d0 )*X*PLMIN1
      PLM = PLM - PLMIN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      PLM = PLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
