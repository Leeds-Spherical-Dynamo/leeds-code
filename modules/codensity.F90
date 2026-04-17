!> Constains all subroutines applicable to entropy/codensity.
!!
!! codensity becomes entropy in the anelastic approximation, but the
!! variable is still called codensity
 module codensity
   use variables, only : coll, var_H
   use meshs, only : lumesh, mesh, rdom
   use parameters, only : i_L1
   implicit none
   save

   type (coll), private :: C !< temporary (coll) for computations 
   type (coll), private :: NC !< temporary (coll) for computations
   
   type (lumesh), private :: XC(0:i_L1) !< cod lu mesh
   type (mesh),   private :: YC(0:i_L1) !< cod Y mesh

   type (coll), private :: cC !< temporary (coll) for
   
 contains

!------------------------------------------------------------------------
!>  Initialise codensity
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine cod_precompute(OC_dom, codensity, source, BC_Re, BC_Im)
      use variables, only : var_coll_init
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(coll), intent(out) :: codensity !< codensity in spectral space
      double precision, intent(out) :: source(:) !< codensity source
      double precision, intent(out) :: BC_Re(:,0:) !< Real part of boundary conditions
      double precision, intent(out) :: BC_Im(:,0:) !< imaginary part of boundary conditions
      call var_coll_init(OC_dom,var_H, codensity)
      source    = 0d0
      BC_Re = 0d0
      BC_Im = 0d0
   end subroutine cod_precompute


!------------------------------------------------------------------------
!>  Precomputation of codensity timestepping matrices
!!
!! Sets up codensity private variables: XC, YC
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine cod_matrices(OC_dom)
      use parameters, only : d_Pr, d_Pm, b_boussinesq, smart_hd_cod, qh_cod, ldiff_cod
      use timestep, only : tim_lumesh_x, tim_mesh_Y
      use meshs, only : mes_lu_find
      type(rdom), intent(in) :: OC_dom !< outer core domain
      integer :: l
      double precision :: c1, c2, hyp_factor, c2_default
      c1 = d_Pr/d_Pm
      c2_default = 1d0

      do l = 0, i_L1
         ! smart hyperdiffusion for codensity
         if (b_boussinesq .and. smart_hd_cod) then
            if (l>=ldiff_cod) then
               hyp_factor = qh_cod**(l-ldiff_cod)
               c2 = hyp_factor * c2_default
            else
               c2 = c2_default
            end if
         else
            c2 = c2_default
         end if

         call tim_lumesh_X      ( OC_dom, c1, c2, l, XC(l) )
         call cod_bc_C          ( OC_dom, XC(l), l)
         call mes_lu_find       ( OC_dom, XC(l) )
         call tim_mesh_Y        ( OC_dom, c1, c2, l, YC(l) )
      end do

   end subroutine cod_matrices


!-------------------------------------------------------------------------
!> Set boundary conditions on Codensity
!!
!! Boundary conditions set on codensity lumesh
!! depending on parameter i_cod_bc
!!________________________________________________________________________
!-------------------------------------------------------------------------
   subroutine cod_bc_C(D,A,l)
      use parameters, only : i_cod_bc, i_Kl
      type (rdom),   intent(in)  :: D !< outer/inner core domain
      type (lumesh), intent(out) :: A !< lu mesh (for codensity)
      integer,       intent(in)  :: l
      integer :: j

      if(i_cod_bc==1 .or. i_cod_bc==2) then	
         A%M(2*i_KL+1,1)   = 1d0
      else
         do j = 1, 1+i_KL
            A%M(2*i_KL+1+1-j,j) = D%dr(1)%M(i_KL+1+1-j,j)
         end do
      end if
      
      if(i_cod_bc==1 .or. i_cod_bc==3) then	
         A%M(2*i_KL+1,D%N) = 1d0
      else
         do j = D%N-i_KL, D%N
            A%M(2*i_KL+1+D%N-j,j) = D%dr(1)%M(i_KL+1+D%N-j,j)
         end do
      end if
! # Rob (14/02/12): Fix temperature so that the code knows what the temperature is when both boundaries are fixed flux
! This is done by setting the average (i.e. l=0) mode as a BC on T rather than on dT/dr
! WD 30/09/2013 Fixing the upper boundary temperature in l=0-mode
!      if((i_cod_bc==4) .and. (l==0)) then
!         A%M(2*i_KL+1,D%N)   = 1d0
!         do j = D%N-i_KL, D%N
!            A%M(2*i_KL+1+D%N-j,j) = 0d0
!         end do
!         A%M(2*i_KL+1,D%N)   = 1d0
!      endif
! ###
      
   end subroutine cod_bc_C


!-------------------------------------------------------------------------
!>  set the RHS for the boundary condition = 0
!!________________________________________________________________________
!-------------------------------------------------------------------------
   subroutine cod_setbc(BC_Re, BC_Im, a)
      use parameters, only : i_inhombc_C, i_cod_bc, d_qi, d_qo
      use dyn_mpi, only : mpi_rnk
      double precision, intent(in) :: BC_Re(:,0:) !< Real part of boundary conditions (harmonics)
      double precision, intent(in) :: BC_Im(:,0:) !< imaginary part of boundary conditions (harmonics)
      type (coll), intent(inout) :: a !< collocated variable to apply BCs to
      integer :: nh, n

      n = a%D%N

      ! If i_inhombc_C == 2 or 3 then the boundary values of the codensity (entropy) or its
      ! gradient are read in from either codBC.cdf.in or codBC.txt.in, which must be present in the run
      ! directory. Loaded in at io_load_codBC

      if ( (i_inhombc_C==2) .or. (i_inhombc_C==3) ) then
         do nh = 0, a%H%pH1
            a%Re( 1, nh ) = BC_Re(1,nh)
            a%Re( n, nh ) = BC_Re(2,nh)
            a%Im( 1, nh ) = BC_Im(1,nh)
            a%Im( n, nh ) = BC_Im(2,nh)
         end do

      else if (i_inhombc_C==1) then
         do nh = 0, a%H%pH1
            ! Initialise the right-hand sideto zero
            a%Re( 1, nh ) = 0.d0
            a%Re( n, nh ) = 0.d0
            a%Im( 1, nh ) = 0.d0
            a%Im( n, nh ) = 0.d0
         end do
          ! If i_inhomBC == 1 then codensity (entropy) or its gradient must be
          ! given a uniform value on a boundary.


          ! Fixed flux is controlled by parameters d_qi and d_qo
          ! given in parameters.F90. io_cod_bc =1 fixed codensity (entropy) 
          ! on r_i, r_o. io_cod_bc=2, fixed codensity 1 on r_i fixed flux (gradient) d_qo on r_o
          ! io_cod_bc=3, fixed codensity = 0 on r_o, fixed flux d_qi on r_i, 
          ! io_cod_bc=4 fixed gradient d_qi at r_i, d_qo at r_o. 
          if (mpi_rnk==0) then
            ! Set bottom boundary
            if (i_cod_bc==1 .or. i_cod_bc==2) then ! inner boundary fixed codensity
               a%Re(1, 0) = 1.0d0
            else
               if (i_cod_bc==3 .or. i_cod_bc==4) then ! inner boundary FF
                 a%Re(1, 0) = d_qi
               endif
            endif
            ! Set top boundary
            ! Now setting top boundary to be d_qo if FFFF boundary conditions.
            if (i_cod_bc==2 .or. i_cod_bc==4 ) then ! outer boundary FF
               a%Re(n, 0) = d_qo
            else
               a%Re(n, 0) = 0.d0 ! outer boundary Fixed codensity
            endif
          endif
      end if

   end subroutine cod_setbc


!------------------------------------------------------------------------
!>  Get r, theta, phi components of gradient of codensity
!!
!! Get radial gradient from finite difference, then transform to phys
!! Change codensity coll to spec, then calculate th, phi gradients using
!! tra_grad().
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine cod_transform(OC_dom, codensity, grad_C_r, grad_C_th, grad_C_phi)
      use variables, only : phys, spec, coll, var_coll2spec, var_coll_meshmult
      use transform, only : tra_spec2phys, tra_grad
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(coll), intent(in) :: codensity !< codensity in spectral space
      type(phys), intent(out) :: grad_C_r !< grad(C), radial component
      type(phys), intent(out) :: grad_C_th !< grad(C), theta component
      type(phys), intent(out) :: grad_C_phi !< grad(C), phi component
      type (spec) :: sC, sC_

      call var_coll_meshmult(OC_dom%dr(1), codensity, cC)

      call var_coll2spec(codensity, sC, c2=cC,s2=sC_) 

      call tra_spec2phys(sC_, grad_C_r)

      call tra_grad(sC, grad_C_th, grad_C_phi)
      
   end subroutine cod_transform


!------------------------------------------------------------------------
!> Codensity Predictor 
!!
!! - N  := cN,                         save N at time t
!! - C  := Y cC + cN,   C := invX C,   get prediction C*
!! - cC := C                           copy prediction
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine cod_predictor(codensity_n, BC_Re, BC_Im, Codensity)
      use timestep, only : tim_multY, tim_invx
      use variables, only : var_coll_copy
      type(coll), intent(in) :: codensity_n !< codensity nonlinear term
      double precision, intent(in) :: BC_Re(:,0:) !< Real part of boundary conditions (harmonics)
      double precision, intent(in) :: BC_Im(:,0:) !< Imaginary part of boundary conditions (harmonics)
      type(coll), intent(inout) :: Codensity  !< codensity in spectral space
      call var_coll_copy (codensity_n, NC)
      call tim_multY     (.true.,YC,Codensity,codensity_n, C)
      call cod_setbc     (BC_Re, BC_Im, C)
      call tim_invX      (.true.,XC, C)
      call var_coll_copy (C, Codensity)
   end subroutine cod_predictor
   
   
!------------------------------------------------------------------------
!>  Codensity Corrector
!!
!!  - C  := c (cN - N),   	using N* get nlin correction
!!  - N  := cN		save last N
!!  - C  := invX C,   	get correction to C, correction has 0 bc
!!  - cC := cC + C		update correction
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine cod_corrector(codensity_n, codensity)
      use timestep, only : tim_nlincorr, tim_zerobc, tim_invx, tim_addcorr
      use variables, only : var_coll_copy
      type(coll), intent(in) :: codensity_n !< codensity nonlinear term
      type(coll), intent(inout) :: codensity !< codensity in spectral space
      call tim_nlincorr (codensity_n,NC, C)
      call var_coll_copy(codensity_n,NC)
      call tim_zerobc   (C)
      call tim_invX     (.true.,XC, C)
      call tim_addcorr  (C, codensity)

   end subroutine cod_corrector
   
   
!*************************************************************************
 end module codensity
!*************************************************************************
