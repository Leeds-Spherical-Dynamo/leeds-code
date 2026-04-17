C*********************************************************************
C function PMM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C> Gives the Schmidt Normalised Legendre Function P_m^m ( X )         
C> from eqn. 175 in my notes ie
C> \f[
C> P_m^m(x) = \frac{ ( 2m - 1)!! (1- XX)^{m/2} * \sqrt{2}}
C>   {\sqrt{(2m)!}}
C> \f]
C>                                                                    
C> for m non-zero and                                           
C>                                                                    
C> \f[  P_0^0( X ) = 1.0d0 \f]                                               
C>                                                                    
C> N.B. The double factorial sign means the product of all ODD        
C> integers between 1 and ( 2m - 1 ).                                 
C> Best to use the form in eq 184 to calculate PMM.                   
C>____________________________________________________________________
      FUNCTION PMM ( M, S )
      IMPLICIT NONE
C> Pmm(x)
      DOUBLE PRECISION PMM
C> Order of polynomial
      INTEGER M
C> Sin(theta)
      DOUBLE PRECISION S
      INTEGER I
      DOUBLE PRECISION TWOI
C START OF FUNCTION
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         PMM = 1.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      PMM = DSQRT ( 2.0d0 )
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         PMM = PMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
