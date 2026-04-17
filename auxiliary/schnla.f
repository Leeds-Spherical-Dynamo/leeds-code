C*********************************************************************
C subroutine SCHmidt Normalised Legendre function Array **************
C            ---     -          -                 -     **************
C Steve Gibbons 22.4.97                                              C
C____________________________________________________________________C
C> Does the same as SCHNLF except that instead of a single valued X   
C> for one theta point, it fills arrays PA and DPA with the           
C> legendre Functions etc. for each of the NTHPTS values of cos(theta)
C> in the array GAUX.                                                 
C>____________________________________________________________________
C> Functions ... Calling Proceedures:                              
C> Double Precision                                                  
C>_________________                                       
C> PMM ( M, S )					                     
C> DPMM ( M, C, S )				                     
C> PMM1 ( M, X, PMM0 )                                                
C> PLM ( L, M, X, PLMIN1, PLMIN2 )				     
C> DPMM1 ( M , X , S, PMM , DPMM)                                     
C> DPLM ( L, M , X , S, PMM1 , DPMM1, DPMM2 )                         
C>____________________________________________________________________
      SUBROUTINE SCHNLA ( PA, DPA, GAUX, LH, NTHPTS )
      IMPLICIT NONE
C> Highest degree, l, of spherical harmonic.
      INTEGER LH
C> Number of theta points.
      INTEGER NTHPTS
C> Array of cosines to the NTHPTS angles. 
      DOUBLE PRECISION GAUX( NTHPTS )
C> Schmidt Normalised Legendre Functions Dimension.
C> {  ( LH + 1 )*( LH + 2 )/2 , NTHPTS }
C> P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C> Derivatives of the above.
      DOUBLE PRECISION DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C Variable declarations - Working variables
      INTEGER INDEX,L,M,IOLD1,IOLD2,NTHETA
      DOUBLE PRECISION SINE,PMIN1,PMIN2,TOL,DPMIN1,
     1                 DPMIN2,X
      PARAMETER (TOL=1.0d-6)
C Variable declarations - Functions called 
      DOUBLE PRECISION PMM,PMM1,PLM,DPMM,DPMM1,DPLM
C START OF PROGRAM 
C ........................... now loop around theta points
      DO NTHETA = 1, NTHPTS
C
        X = GAUX ( NTHETA )
        SINE = X*X
        IF ( SINE.GT.1.0d0 ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' Illegal Cos(theta) has been entered.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C
C Set SINE (theta) in terms of X
        SINE = DSQRT ( (1.0d0 + X)*(1.0d0 - X) )
        IF ( SINE.LT.TOL ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' SINE is too small. Division by zero imminent.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C..................... first calculate the P_l^m ..........
        DO M = 0, LH - 2
C                        ............. Calculate P_M^M .....
           L = M
           INDEX = L*(L+1)/2+M+1
           PA ( INDEX , NTHETA) = PMM ( M , SINE )
           DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C                         ............. Calculate P_(M+1)^M .
           PMIN1 = PA ( INDEX , NTHETA)
           DPMIN1 = DPA ( INDEX , NTHETA)
           IOLD2 = INDEX
           L = L + 1
           INDEX = L*(L+1)/2+M+1
           PA (INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
           DPA (INDEX , NTHETA) = DPMM1 (M,X,SINE,PMIN1, DPMIN1)
           IOLD1 = INDEX
C                         ......... Calculate P_L^M general .
           DO L = M + 2, LH
              PMIN2 = PA ( IOLD2 , NTHETA)
              PMIN1 = PA ( IOLD1 , NTHETA)
              DPMIN2 = DPA ( IOLD2 , NTHETA)
              DPMIN1 = DPA ( IOLD1 , NTHETA)
              INDEX = L*(L+1)/2+M+1
              PA ( INDEX , NTHETA) = PLM ( L,M,X,PMIN1,PMIN2 )
              DPA ( INDEX , NTHETA) = DPLM (L,M,X,SINE , PMIN1, 
     1                              DPMIN1, DPMIN2 )
              IOLD2 = IOLD1
              IOLD1 = INDEX
           ENDDO
        ENDDO
        M = LH - 1
        L = M
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
        PMIN1 = PA( INDEX , NTHETA)
        DPMIN1 = DPA( INDEX , NTHETA)
        L = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
        DPA(INDEX ,NTHETA) = DPMM1 (M ,X ,SINE,PMIN1,DPMIN1 )
        M = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C......................finished calculating P_l^m .........
      ENDDO
C......................finished looping around theta points

      RETURN
      END
C*********************************************************************
