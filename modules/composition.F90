!> Subroutines for composition, used for 2-component convection
module composition
   use variables, only : coll, var_H
   use meshs, only : lumesh, mesh, rdom
   use parameters, only : i_L1
   implicit none
   save

   type (coll), private :: C !< temporary (coll) for computations 
   type (coll), private :: NC !< temporary (coll) for computations
   
   type (lumesh), private :: XC(0:i_L1) !< cod lu mesh
   type (mesh),   private :: YC(0:i_L1) !< cod Y mesh

   type (coll), private :: cC !< temporary (coll) for computations
   
 contains

!------------------------------------------------------------------------
!>  initialise Composition
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine comp_precompute(OC_dom, composition, comp_Nc, source, BC_Re, BC_Im)
      use variables, only : var_coll_init
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(coll), intent(out) :: composition !< composition in spectral space
      type(coll), intent(out) :: comp_Nc !< composition nonlinear term
      double precision, intent(out) :: source(:) !< composition source
      double precision, intent(out) :: BC_Re(:,0:) !< Real part of boundary conditions
      double precision, intent(out) :: BC_Im(:,0:) !< imaginary part of boundary conditions
      call var_coll_init(OC_dom,var_H, composition)
      call var_coll_init(OC_dom,var_H, comp_NC)
      source    = 0d0
      BC_Re = 0d0
      BC_Im = 0d0
   end subroutine comp_precompute


!------------------------------------------------------------------------
!>  precomputation of composition timestepping matrices
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine comp_matrices(OC_dom)
      use parameters, only : d_Sc, d_Pm, b_boussinesq, smart_hd_comp, qh_comp, ldiff_comp
      use timestep, only : tim_lumesh_x, tim_mesh_Y
      use meshs, only : mes_lu_find
      type(rdom), intent(in) :: OC_dom !< outer core domain
      integer :: l
      double precision :: c1, c2, hyp_factor, c2_default
      c1 = d_Sc/d_Pm
      c2_default = 1d0

      do l = 0, i_L1
         ! smart hyperdiffusion for composition
         if (b_boussinesq .and. smart_hd_comp) then
            if (l>=ldiff_comp) then
               hyp_factor = qh_comp**(l-ldiff_comp)
               c2 = hyp_factor * c2_default
            else
               c2 = c2_default
            end if
         else
            c2 = c2_default
         end if

         call tim_lumesh_X      ( OC_dom, c1, c2, l, XC(l) )
         call comp_bc_C          ( OC_dom, XC(l), l)
         call mes_lu_find       ( OC_dom, XC(l) )
         call tim_mesh_Y        ( OC_dom, c1, c2, l, YC(l) )
      end do

   end subroutine comp_matrices


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!> Set boundary conditions on Composition
!>
!! Boundary conditions set on composition lumesh
!! depending on parameter i_comp_bc 
!!________________________________________________________________________
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   subroutine comp_bc_C(D,A,l)
      use parameters, only : i_comp_bc, i_Kl
      type (rdom),   intent(in)  :: D !< outer/inner core domain
      type (lumesh), intent(out) :: A !< lu mesh (for composition)
      integer,       intent(in)  :: l
      integer :: j

      if(i_comp_bc==1 .or. i_comp_bc==2) then
         A%M(2*i_KL+1,1)   = 1d0
      else
         do j = 1, 1+i_KL
            A%M(2*i_KL+1+1-j,j) = D%dr(1)%M(i_KL+1+1-j,j)
         end do
      end if
      
      if(i_comp_bc==1 .or. i_comp_bc==3) then
         A%M(2*i_KL+1,D%N) = 1d0
      else
         do j = D%N-i_KL, D%N
            A%M(2*i_KL+1+D%N-j,j) = D%dr(1)%M(i_KL+1+D%N-j,j)
         end do
      end if
      
   end subroutine comp_bc_C


!-------------------------------------------------------------------------
!>  set the RHS for the boundary condition = 0
!!________________________________________________________________________
!-------------------------------------------------------------------------
   subroutine comp_setbc(comp_bcRe, comp_bcIm,a)
      use parameters, only : i_inhombc_Comp, i_comp_bc, d_qi_comp, d_qo_comp
      use dyn_mpi, only : mpi_rnk
      double precision, intent(in) :: comp_bcRe(:,0:) !< Real part of boundary conditions (harmonics)
      double precision, intent(in) :: comp_bcIm(:,0:) !< imaginary part of boundary conditions (harmonics)
      type (coll), intent(inout) :: a !< collocated variable to apply BCs to
      integer :: nh, n

      n = a%D%N

      ! If i_inhombc_comp == 2 or 3 then the boundary values of the composition or its
      ! gradient are read in from either compBC.cdf.in or compBC.txt.in, which must be present in the run
      ! directory. Loaded in at io_load_compBC

      if ( (i_inhombc_comp==2) .or. (i_inhombc_comp==3) .or. (i_inhombc_comp==4)) then
         do nh = 0, a%H%pH1
            a%Re( 1, nh ) = comp_bcRe(1,nh)
            a%Re( n, nh ) = comp_bcRe(2,nh)
            a%Im( 1, nh ) = comp_bcIm(1,nh)
            a%Im( n, nh ) = comp_bcIm(2,nh)
         end do

      else if (i_inhombc_comp==1) then
         do nh = 0, a%H%pH1
            ! Initialise the right-hand side to zero
            a%Re( 1, nh ) = 0.d0
            a%Re( n, nh ) = 0.d0
            a%Im( 1, nh ) = 0.d0
            a%Im( n, nh ) = 0.d0
         end do
          ! If i_inhomBC == 1 then compostion or its gradient must be
          ! given a uniform value on a boundary.


          ! Fixed flux is controlled by parameters d_qi and d_qo
          ! given in parameters.F90. io_comp_bc =1 fixed composition
          ! on r_i, r_o. io_comp_bc=2, fixed composition 1 on r_i fixed flux (gradient) d_qo on r_o
          ! io_comp_bc=3, fixed composition = 0 on r_o, fixed flux d_qi on r_i,
          ! io_comp_bc=4 fixed gradient d_qi at r_i, d_qo at r_o.
          if (mpi_rnk==0) then
            ! Set bottom boundary
            if (i_comp_bc==1 .or. i_comp_bc==2) then ! inner boundary fixed composition
               a%Re(1, 0) = 1.0d0
            else
               if (i_comp_bc==3 .or. i_comp_bc==4) then ! inner boundary FF
                 a%Re(1, 0) = d_qi_comp
               endif
            endif
            ! Set top boundary
            ! Now setting top boundary to be d_qo if FFFF boundary conditions.
            if (i_comp_bc==2 .or. i_comp_bc==4 ) then ! outer boundary FF
               a%Re(n, 0) = d_qo_comp
            else
               a%Re(n, 0) = 0.d0 ! outer boundary Fixed composition
            endif
          endif
      end if

   end subroutine comp_setbc


!------------------------------------------------------------------------
!>  Get r, theta, phi components of gradient of composition
!!
!! Get radial gradient from finite difference, then transform to phys
!! Change composition coll to spec, then calculate th, phi gradients using
!! tra_grad().
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine comp_transform(OC_dom, composition, grad_comp_r, grad_comp_th, grad_comp_phi)
      use variables, only : phys, spec, coll, var_coll2spec, var_coll_meshmult
      use transform, only : tra_spec2phys, tra_grad
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(coll), intent(in) :: composition  !< composition in spectral space
      type(phys), intent(out) :: grad_comp_r !< grad(comp), radial component
      type(phys), intent(out) :: grad_comp_th !< grad(comp), theta component
      type(phys), intent(out) :: grad_comp_phi !< grad(comp), phi component
      type (spec) :: sC, sC_
      call var_coll_meshmult(OC_dom%dr(1),composition, cC)

      call var_coll2spec(composition,sC, c2=cC,s2=sC_)

      call tra_spec2phys(sC_, grad_comp_r)

      call tra_grad(sC, grad_comp_th,grad_comp_phi)
      
   end subroutine comp_transform


!------------------------------------------------------------------------
!> Composition predictor
!!
!!  N  := cN,                         save N at time t
!!  C  := Y cC + cN,   C := invX C,   get prediction C*
!!  cC := C                           copy prediction
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine comp_predictor(composition_n, BC_Re, BC_Im, composition)
      use timestep, only : tim_multY, tim_invx
      use variables, only : var_coll_copy
      type(coll), intent(in) :: composition_n !< composition nonlinear term
      double precision, intent(in) :: BC_Re(:,0:) !< Real part of boundary conditions (harmonics)
      double precision, intent(in) :: BC_Im(:,0:) !< Imaginary part of boundary conditions (harmonics)
      type(coll), intent(inout) :: composition !< composition in spectral space

      call var_coll_copy (composition_n, NC)
      call tim_multY     (.true.,YC,composition,composition_n, C)
      call comp_setbc     (BC_Re, BC_Im, C)
      call tim_invX      (.true.,XC, C)
      call var_coll_copy (C, composition)
      
   end subroutine comp_predictor
   
   
!------------------------------------------------------------------------
!> Composition Corrector
!!
!!  C  := c (cN - N),   	using N* get nlin correction
!!  N  := cN		save last N
!!  C  := invX C,   	get correction to C, correction has 0 bc
!!  cC := cC + C		update correction
!!________________________________________________________________________
!------------------------------------------------------------------------
   subroutine comp_corrector(composition_n, composition)
      use timestep, only : tim_nlincorr, tim_zerobc, tim_invx, tim_addcorr
      use variables, only : var_coll_copy
      type(coll), intent(in) :: composition_n !< composition nonlinear term
      type(coll), intent(inout) :: composition !< composition in spectral space
      call tim_nlincorr (composition_n,NC, C)
      call var_coll_copy(composition_n,NC)
      call tim_zerobc   (C)
      call tim_invX     (.true.,XC, C)
      call tim_addcorr  (C, composition)

   end subroutine comp_corrector
   
   
!*************************************************************************
 end module composition
!*************************************************************************
