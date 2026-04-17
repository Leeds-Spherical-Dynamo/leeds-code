#include "../parallel.h"
!**************************************************************************
!> Definition of mesh-types.
!!
!! Public variables:
!!~~~
!!   r(n,k)              radial points r_j^k, j in [1,N]
!!
!!   dr(p)%M             p-th derivative of f = dr(p)%M . f
!!
!!   radLap%M            radial part of laplacian  drr + (2/r) dr
!!~~~
!!
!! (mesh)-type matrices are stored s.t.   M(   KL+1+n-j, j) = A(n,j)
!! (lumesh)-type matrices are stored s.t. M( 2*KL+1+n-j, j) = A(n,j)
!!
!! Note that the index n is usually used to denote r_n and the operator on
!! the n-th point corresponds to the n-th row of the matrix.  The index
!! j usually denotes the column and is therefore in the outer loop below.
!! Matrix multiplication  y := M . x  can be achieved accessing M only
!! once along the memory by the loop form
!! ~~~
!!   y = 0
!!   do j = 1, N
!!      do n = max(1,j-KL), min(j+KL,N)
!!         y(n) = y(n) + M(KL+1+n-j,j) * x(j)
!! ~~~
!!**************************************************************************
 module meshs
!**************************************************************************
   use parameters, only : i_Kl, i_N, i_Nic
   implicit none
   save
                                ! M(KL+1+n-j, j) = A(n,j)
   !> (mesh)-type matrices are stored as\n
   !! M(   resolution::i_Kl+1+n-j, j) = A(n,j)\n
   type mesh
      double precision :: M(2*i_KL+1, i_N)
   end type mesh
   
   !> (lumesh)-type for outer core 
   !!
   !! matrices are stored s.t. M( 2*KL+1+n-j, j) = A(n,j)\n
   !! M(2*KL+1+n-j, j) = A(n,j)\n
   !! see lapack dgbtrf 
   type lumesh
      integer          :: ipiv(i_N)
      double precision :: M(3*i_KL+1, i_N)
   end type lumesh

   !> lumesh type for inner core.
   !!
   !! matrices are stored s.t. M( 2*KL+1+n-j, j) = A(n,j)\n
   !! M(2*KL+1+n-j, j) = A(n,j)\n
   !! see lapack dgbtrf 
   type iclumesh
      integer          :: ipiv(i_Nic+i_N-1)
      double precision :: M(3*i_KL+1, i_Nic+i_N-1)
   end type iclumesh

   integer, parameter  :: i_rpowmax =  3 !< 3, max power of radial points
   integer, parameter  :: i_rpowmin = -3 !< -3, min power for radial points

   !> The "rdom" type defines finite difference radial mesh.\n
   !! It defines, not only the positions, but also the differentiation and integration weights in the radial mesh
   type rdom
      integer          :: N !< Total number of radial points
      integer          :: pNi !< Location of the first (innermost) radial point treated by the current CPU.
      integer          :: pN !< Number of radial points treated by the present CPU. 
      integer          :: pNi_(0:_Nr-1) !< array, indices of first radial point for each process (rank)
      integer        :: pN_(0:_Nr-1) !< array, number of radial points for each process (rank)
      !> Radial values and their powers: r(n,p) = \f$ (r_n)^p\f$.
      double precision :: r(i_N,i_rpowmin:i_rpowmax) 
      !> Flayer enhanced diffusion based on radial position
      double precision :: f_factor(i_N)
! Modified to add in compressible attributes      
      double precision :: xi(i_N,0:1),zeta(i_N,0:1)
      double precision :: rho(i_N,0:2),Temp0(i_N,0:1)
      double precision :: eta(i_N,0:1)
      double precision :: intr2dr(i_N)  !< Weights for volume integrals.
! weights for integrating with rho (for velocity integrals)
      double precision :: intr2dr_comp(i_N)
      double precision :: intr2dr_compS(i_N), intr2dr_magdiss(i_N)
      double precision :: intr2dr_veldiss(i_N)
      double precision :: intr2drh(i_N,0:2*i_KL)
      !> Weights for taking the \f$ p^{th} \f$ derivative.
      type (mesh)      :: dr(i_KL) 
      type (mesh)      :: radLap !< Weights for the radial part of the laplacian.
      type (mesh)      :: radLap_flayer !< radial laplacian with enhanced diffusion in flayer
!
!  Introduce the compressible Laplacians needed
!
      type (mesh)      :: CompLap
      type (mesh)      :: GreLap
      type (mesh)      :: r1dr !< Weights for 1/r dX/dr
      type (mesh)      :: compr1dr

!
!  Introduce the variable conductivity Laplacians needed
!
      type (mesh)      :: MagTorLap
      type (mesh)      :: MagPolLap
   end type rdom

 contains


!------------------------------------------------------------------------
!>  Precompute mesh variables
!!
!!  Sets radial collocation points for outer and inner core meshs, 
!!  and initialises @link mes_oc @endlink and @link mes_ic @endlink
!!  using @link mes_rdom_init @endlink
!!
!------------------------------------------------------------------------
   subroutine mes_precompute(mes_oc, mes_ic)
      use parameters, only : d_alpha, d_rratio, d_flayer_radius, d_Pi, b_boussinesq,&
                             & b_diff_profile
      use dyn_mpi, only : mpi_rnk
      type(rdom), intent(inout) :: mes_oc !< outer core domain
      type(rdom), intent(inout) :: mes_ic !< inner core domain
      double precision :: ri,min_diff
      integer :: n, N_, sta_input_file,n_closest


      ! Initialize
      min_diff = 1.0d99
      n_closest = -1
      ! get collocation pts
      ! N extrema of T_{N-1}(x)
                                
      ! outer core
      if(_Nr>i_N) stop '_Nr>i_N'
      mes_oc%N   = i_N
      do n = 1, i_N
         if (d_alpha>0) then
            mes_oc%r(n,1) = 0.5d0 * (1d0 + dasin(d_alpha * dcos(d_Pi*(i_N-n)/dble(i_N-1)))/dasin(d_alpha) )
         else
            mes_oc%r(n,1) = 0.5d0 * ( 1d0 + dcos(d_PI*(i_N-n)/dble(i_N-1)) )
         end if
      end do      
      ri  = d_rratio/(1d0-d_rratio)
      mes_oc%r(:,1) = mes_oc%r(:,1) + ri
      
     ! if (mpi_rnk==0) print *, mes_oc%r(:,1)
     ! if (mpi_rnk==0) print *, d_flayer_radius
                                 
      if (d_flayer_radius > 0) then
         do n = 1, i_N
            if (abs(mes_oc%r(n,1) - (d_flayer_radius )) < min_diff) then
               min_diff = abs(mes_oc%r(n,1) - d_flayer_radius )
               n_closest = n
            end if
         end do
         !if (mpi_rnk==0) then
         !    print *, "Closest index to rs is:", n_closest
         !   print *, "Original value at that index:", mes_oc%r(n_closest,1)
         !end if
      ! Replace the closest value with rs
          if (n_closest > 0) then
            mes_oc%r(n_closest,1) = d_flayer_radius
         end if   
      end if  
      !if (mpi_rnk==0) print *, mes_oc%r(:,1)                         
      mes_ic%N = i_Nic
      N_ = 2*i_Nic
      do n = i_Nic+1, N_
         mes_ic%r(n-i_Nic,1) =  ri * dcos(d_PI*(N_-n)/dble(N_-1))
      end do
!
!  Write out the anelastic reference state to file
!

      if (.not.b_boussinesq) then
        if (mpi_rnk==0) then
          open(15,file='Reference_state.txt',status='Replace')
          write(15,'(8A14)') '      r      ','      rho     ','      xi      ','   d xi / dr  ', &
            '     temp     ','  d temp / dr  ','      eta     ','   d eta / dr '
          do n=1,i_N
            write(15,'(8e14.6)') mes_oc%r(n,1), mes_oc%rho(n,0), mes_oc%rho(n,1), &
             mes_oc%rho(n,2),  mes_oc%Temp0(n,0), mes_oc%Temp0(n,1), &
             mes_oc%eta(n,0),  mes_oc%eta(n,1)
          end do
          if (i_Nic/=1) then
              write(15,'(A20)') ' Inner core present '
              write(15,'(3A16)') '         r       ','    Inner eta   ' , 'Inner d eta / dr'
              do n=1,i_Nic
                write(15,'(3e16.6)') mes_ic%r(n,1),mes_ic%eta(n,0), mes_ic%eta(n,1)
              end do
          end if
        end if
      end if

!
!  Write out the Boussinesq diffusivity profile if it was read in from
!  diffusivity.txt
!

      if (b_boussinesq.and.b_diff_profile) then
        if (mpi_rnk==0) then
          open(16,file='Diffusivity_profile.dat',status='Replace')
          write(16,'(3A14)') '      r      ','      eta     ','   d eta / dr '
          do n=1,i_N
            write(16,'(3e14.6)') mes_oc%r(n,1), mes_oc%eta(n,0),  mes_oc%eta(n,1)
          end do
          if (i_Nic/=1) then
              write(16,'(A20)') ' Inner core present '
              write(16,'(3A16)') '         r       ','    Inner eta   ' , 'Inner d eta / dr'
              do n=1,i_Nic
                write(16,'(3e16.6)') mes_ic%r(n,1),mes_ic%eta(n,0), mes_ic%eta(n,1)
              end do
          end if
        end if
      end if

      call mes_rdom_init(mes_oc)
      call mes_rdom_init(mes_ic)

   end subroutine mes_precompute


!------------------------------------------------------------------------
!>  Set up radial grid for mesh object D.
!!
!! for the radial domain D, first set the number of points D%N and their
!! values D%r, then mes_rdom_init() will fill in the rest...
!------------------------------------------------------------------------
   subroutine mes_rdom_init(D)
      use dyn_mpi, only : r_rnk
      use parameters, only : d_enhance_diff, d_flayer_radius, d_transition_width 
      type (rdom), intent(inout) :: D !< inner/outer core domain
      !> num filled rows in finite difference matrix
      integer, parameter :: lda = 2*i_KL+1 
      double precision   :: A(lda,lda), c,e,w, cc, ss
      integer :: n,i,j,nn, l,r,p


      ! WD 03/09/2013, as for initialisation of harm type and 
      !the distribution over the cores, also the radial domain needs some update

      D%pN_ = 0
      ! find number of radial points on each process
      do n = 1, D%N
         r = _Nr - modulo(n-1,_Nr) - 1
         D%pN_(r) = D%pN_(r) + 1
      end do
      ! find index of first radial point for each process
      D%pNi_(0) = 1
      do r = 1, _Nr-1
         D%pNi_(r) = D%pNi_(r-1) + D%pN_(r-1)
      end do

      ! set this processor's num radial points and first point
      p = r_rnk
      D%pN  = D%pN_ (p)
      D%pNi = D%pNi_(p)


         			! get ad_r(:,j) = r^j
      do j = i_rpowmin, i_rpowmax
         if(j==1) cycle
         D%r(1:D%N,j) = D%r(1:D%N,1)**j
      end do

      ! get finite difference coeffs    
      if(i_KL<2) stop 'KL<2'
      do i = 1, i_KL
         D%dr(i)%M = 0d0
      end do
      D%radLap%M = 0d0
      D%radLap_flayer%M = 0d0
!
! Extra laplacians for compressibility here
! 
      D%CompLap%M = 0d0
      D%GreLap%M = 0d0
      D%r1dr%M   = 0d0   
      D%compr1dr%M   = 0d0
      D%MagTorLap%M = 0d0
      D%MagPolLap%M = 0d0
      D%intr2dr  = 0d0
      D%intr2dr_comp  = 0d0
      D%intr2dr_compS  = 0d0
      D%intr2dr_magdiss  = 0d0
      D%intr2dr_veldiss  = 0d0
      D%intr2drh = 0d0
!
      ! if mes_ic, with Nic=1, finish here
      if(D%N==1) return

      ! for all radial points
      do n = 1, D%N
         l = max(1,n-i_KL)
         r = min(n+i_KL,D%N)
         nn = r-l+1
         do i = 1, nn
            A(i,1) = 1d0
         end do
         do j = 2, nn
            do i = 1, nn
               A(i,j) = A(i,j-1) * (D%r(l+i-1,1)-D%r(n,1)) / dble(j-1)
            end do
         end do
         call mes_mat_invert(nn,A,lda)
         do i = 1, i_KL
            do j = l, r
               D%dr(i)%M(i_KL+1+n-j,j) = A(i+1,j-l+1)
            end do
         end do
            			! weights for integration r^2 dr
         c = 1d0
         e = 1d0
         do i = 1, nn
            c = c * (D%r(min(n+1,D%N),1)-D%r(n,1)) / dble(i)
            e = e * (D%r(max(  1,n-1),1)-D%r(n,1)) / dble(i)
            do j = l, r
               w = 0.5d0*(c-e)*A(i,j-l+1)*D%r(n,2)
               D%intr2dr (j)     = D%intr2dr (j)     + w
! Extra weights for compressibility here
!      
!              Int r^2 rho dr weights         
               D%intr2dr_comp(j)=D%intr2dr_comp(j) + &
                w*D%rho(n,0)
!              Int r^2 rho T dr weights
               D%intr2dr_compS(j)=D%intr2dr_compS(j) + &
                w*(D%rho(n,0)*D%Temp0(n,0))
!               Int r^2 eta dr  weights
               D%intr2dr_magdiss(j)=D%intr2dr_magdiss(j) + &
                w*D%eta(n,0)
!              Int r^2 2 rho (dxi/dr  - xi^2 / 3 ) dr  wweights
               D%intr2dr_veldiss(j)=D%intr2dr_veldiss(j) + &
                w*2d0*D%rho(n,0)*(D%rho(n,2)-D%rho(n,1)*D%rho(n,1)/3D0)
!
               D%intr2drh(j,i-1) = D%intr2drh(j,i-1) + w
            end do
         end do
      end do
                                ! get  drr + (2/r) dr
      do j = 1, D%N
         l = max(1,j-i_KL)
         r = min(j+i_KL,D%N)
         do n = l, r
            D%radLap%M(i_KL+1+n-j,j)  &
               = D%dr(2)%M(i_KL+1+n-j,j)  &
               + D%dr(1)%M(i_KL+1+n-j,j) * 2d0 * D%r(n,-1)
         end do
      end do
                          ! get  drr + (2/r) dr flayer version
      do j = 1, D%N
         !cc=1d0
         !if (j<i_flayer_radius) cc=d_enhance_diff
         ss = 0.5 + 0.5 * tanh( (D%r(j,1)-d_flayer_radius)/d_transition_width )
         cc = ss + (1-ss)*d_enhance_diff
                        ! get flayer factor for each radial point
         D%f_factor(j) = cc 
         l = max(1,j-i_KL)
         r = min(j+i_KL,D%N)
         do n = l, r
            D%radLap_flayer%M(i_KL+1+n-j,j)  &
               =cc*( D%dr(2)%M(i_KL+1+n-j,j)  &
               + D%dr(1)%M(i_KL+1+n-j,j) * 2d0 * D%r(n,-1))
         end do
      end do

!
! Extra Laplacians for compressibility added here
! 
           ! Also do MagTorLap and MagPolLap here
      do j = 1, D%N
         l = max(1,j-i_KL)
         r = min(j+i_KL,D%N)
         do n = l, r
           ! CompLap: Derivatives for radial part of 1/rho del^2 rho
           ! =  d^2 / dr^2  + (2 xi + 2/r) d / dr
            D%CompLap%M(i_KL+1+n-j,j)  = D%dr(2)%M(i_KL+1+n-j,j)  &
            + 2.0d0*(D%rho(n,1)+D%r(n,-1))*D%dr(1)%M(i_KL+1+n-j,j)
! 
!  MagTorLap  eta*d2 / dr2 + 2eta/r d /dr + d eta/dr * d / dr
! 
            D%MagTorLap%M(i_KL+1+n-j,j)  &
              = D%eta(n,0)*(D%dr(2)%M(i_KL+1+n-j,j) + &
                 2.0d0*D%r(n,-1)*D%dr(1)%M(i_KL+1+n-j,j)) &
                 + D%eta(n,1)*D%dr(1)%M(i_KL+1+n-j,j)
!
!  MagPolLap  eta*d2 / dr2 + 2eta/r d /dr 
!                
            D%MagPolLap%M(i_KL+1+n-j,j)  &
               = D%eta(n,0)*(D%dr(2)%M(i_KL+1+n-j,j) + &
                 2.0d0*D%r(n,-1)*D%dr(1)%M(i_KL+1+n-j,j))
         end do
      end do
!     if (mpi_rnk==0) print*,' TEST MESHS 1',mes_oc%MagTorLap

                ! GreLap: Derivatives d^/dr^2 + (2/r +xi) d/dr
                !
      do j = 1, D%N
         l = max(1,j-i_KL)
         r = min(j+i_KL,D%N)
         do n = l, r
            D%GreLap%M(i_KL+1+n-j,j)  &
             = D%dr(2)%M(i_KL+1+n-j,j)  &
            + (D%rho(n,1)+2.0d0*D%r(n,-1))*D%dr(1)%M(i_KL+1+n-j,j)
         end do
      end do


                               ! get  (1/r + dr)
      D%r1dr%M = D%dr(1)%M
      D%r1dr%M(i_KL+1,:) = D%r1dr%M(i_KL+1,:) + D%r(:,-1)

! get  (1/r + xi+ dr), needed for the compressible P,T to qst transformation
!
      D%compr1dr%M = D%dr(1)%M
      D%compr1dr%M(i_KL+1,:) = D%compr1dr%M(i_KL+1,:) + D%r(:,-1) &
      + D%rho(:,1)



   end subroutine mes_rdom_init



!------------------------------------------------------------------------
!>  Replace nxn matrix A by its inverse
!! @param A is a matrix of dimensions @param lda x @param n
!------------------------------------------------------------------------
   subroutine mes_mat_invert(n,A,lda)
      integer, intent(in) :: n !< matrix is of size lda x n
      integer, intent(in) :: lda !< matrix is of size lda x n
      double precision, intent(inout) :: A(lda,n) !< The matrix to be inverted
      double precision :: work(n)
      integer :: info, ipiv(n)
   
      call dgetrf(n, n, A, lda, ipiv, info)
      if(info /= 0) stop 'matrix inversion error 1'

      call dgetri(n, A, lda, ipiv, work, n, info)
      if(info /= 0) stop 'matrix inversion error 2'

   end subroutine mes_mat_invert


!------------------------------------------------------------------------
!>  Replace A with its LU factorisation
!!  uses lapack routine dgbtrf()
!------------------------------------------------------------------------
   subroutine mes_lu_find(D,A)
      type (rdom),   intent(in)    :: D !< inner/outer core domain
      type (lumesh), intent(inout) :: A !< LU mesh
      integer :: info

      call dgbtrf(D%N, D%N, i_KL, i_KL, A%M, 3*i_KL+1, A%ipiv, info)
!     print *, info
      if(info /= 0) stop 'mes_lu_find'
    
   end subroutine mes_lu_find


!------------------------------------------------------------------------
!>  Replace A with its LU factorisation
!!  uses lapack routine dgbtrf()
!------------------------------------------------------------------------
   subroutine mes_iclu_find(mes_oc,mes_ic,A)
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      type (iclumesh), intent(inout) :: A !< mesh for full core
      integer :: info, N_
      
      N_ = mes_ic%N+mes_oc%N-1
      call dgbtrf(N_, N_, i_KL, i_KL, A%M, 3*i_KL+1, A%ipiv, info)
      if(info /= 0) stop 'mes_iclu_find'

   end subroutine mes_iclu_find


!**************************************************************************
 end module meshs
!**************************************************************************
