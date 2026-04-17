#include "../parallel.h"
!*************************************************************************
!> Infrastructure to deal with the magnetic field and its boundary conditions
!!*************************************************************************
 module magnetic
!*************************************************************************
   use variables, only : coll, spec
   use meshs, only : iclumesh, mesh, rdom
   use parameters, only : i_L1
   !use timestep
   implicit none
   save

   ! _p_m stands for 'private' 'magnetic'
   type (coll), private :: Tor_p_m !< Temporary coll for computation
   type (coll), private :: Pol_p_m !< Temporary coll for computation
   type (coll), private :: NTor !< Temporary coll for computation
   type (coll), private :: NPol !< Temporary coll for computation
   type (coll), private :: icTor !< Temporary coll for computation
   type (coll), private :: icPol !< Temporary coll for computation
   type (coll), private :: icNTor !< Temporary coll for computation
   type (coll), private :: icNPol !< Temporary coll for computation
   
   type (iclumesh), private :: XTor(0:i_L1) !< XTor iclumesh for computation
   type (iclumesh), private :: XPol(0:i_L1) !< XPol iclumesh for computation
   type (mesh),     private :: YTor(0:i_L1) !< YTor mesh for computation
   type (mesh),     private :: YPol(0:i_L1) !< YPol mesh for computation
   type (mesh),     private :: icYTor(0:i_L1) !< icYTor mesh for computation
   type (mesh),     private :: icYPol(0:i_L1) !< icYPol mesh for computation

   type (coll), private :: cq !< coll for computation (qst, q)
   type (coll), private :: cs !< coll for computation (qst, s)
   type (coll), private :: ct !< coll for computation (qst, t)
   type (spec), private :: sq !< spec for computation (qst, q)
   type (spec), private :: ss !< spec for computation (qst, s)
   type (spec), private :: st !< spec for computation (qst, t)
   
 contains

!------------------------------------------------------------------------
!>  initialise magnetic field variables
!------------------------------------------------------------------------
   subroutine mag_precompute(OC_dom, IC_dom, tor, pol,&
                     ictor, icpol, BC_Re, BC_Im,&
                     BC_PolRe, BC_PolIm)
      use variables, only : var_coll_init, var_H
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(rdom), intent(in) :: IC_dom !< inner core domain
      type(coll), intent(inout) :: tor !< toroidal part of field
      type(coll), intent(inout) :: pol !< poloidal part of field
      type(coll), intent(inout) :: ictor !< toroidal part of field (inner core)
      type(coll), intent(inout) :: icpol !< poloidal part of field (outer core)
      double precision, intent(inout) :: BC_Re(:,0:) !< real part of boundary conditions
      double precision, intent(inout) :: BC_Im(:,0:) !< imaginary part of boundary conditions
      double precision, intent(inout) :: BC_PolRe(:,0:) !< real part of poloidal boundary conditions
      double precision, intent(inout) :: BC_PolIm(:,0:) !< imaginary part of poloidal boundary conditions
      call var_coll_init(OC_dom,var_H, tor)
      call var_coll_init(OC_dom,var_H, pol)
      call var_coll_init(IC_dom,var_H, ictor)
      call var_coll_init(IC_dom,var_H, icpol)
      BC_Re = 0d0
      BC_Im = 0d0
      BC_PolRe = 0d0
      BC_PolIm = 0d0
      
   end subroutine mag_precompute
   
   
!------------------------------------------------------------------------
!>  precomputation of magnetic field timestepping matrices
!------------------------------------------------------------------------
   subroutine mag_matrices(OC_dom, IC_dom)
      use parameters, only : b_boussinesq, smart_hd_mag, qh_mag, ldiff_mag, b_diff_profile
      use timestep, only : tim_iclumesh_Pol_X, tim_mesh_mag_Pol_Y, tim_mesh_Y, tim_iclumesh_Tor_x,&
      & tim_mesh_mag_tor_Y
      use meshs, only : mes_iclu_find
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(rdom), intent(in) :: IC_dom !< inner core domain
      integer :: l
      double precision :: c1, c2, c2_default, hyp_factor

      c1 = 1d0
      c2_default = 1d0

      do l = 1, i_L1
         ! smart HD, toroidal
         if (b_boussinesq .and. smart_hd_mag) then
            if (l>=ldiff_mag) then
               hyp_factor = qh_mag**(l-ldiff_mag)
               c2 = hyp_factor * c2_default
            else
               c2 = c2_default
            end if
         else
            c2 = c2_default
         end if
         call tim_iclumesh_Tor_X    (OC_dom, IC_dom, c1, c2, l, XTor(l) )
         call mag_bc_Tor        (OC_dom, IC_dom, XTor(l), l)
         call mes_iclu_find     (OC_dom, IC_dom, XTor(l) )
         if (.not.b_boussinesq.or.b_diff_profile) then
           call tim_mesh_mag_Tor_Y        (OC_dom, OC_dom, c1, c2, l,   YTor(l) )
           call tim_mesh_mag_Tor_Y        (OC_dom, IC_dom, c1, c2, l, icYTor(l) )
         else
           call tim_mesh_Y        ( OC_dom, c1, c2, l,   YTor(l) )
           call tim_mesh_Y        ( IC_dom, c1, c2, l, icYTor(l) )
         end if
      end do

      do l = 1, i_L1
         ! smart HD, poloidal
         if (b_boussinesq .and. smart_hd_mag) then
            if (l>=ldiff_mag) then
               hyp_factor = qh_mag**(l-ldiff_mag)
               c2 = hyp_factor * c2_default
            else
               c2 = c2_default
            end if
         else
            c2 = c2_default
         end if
         call tim_iclumesh_Pol_X    (OC_dom, IC_dom, c1, c2, l, XPol(l) )
         call mag_bc_Pol        (OC_dom, IC_dom, XPol(l), l)
         call mes_iclu_find     (OC_dom, IC_dom, XPol(l) )
         if (.not.b_boussinesq.or.b_diff_profile) then
           call tim_mesh_mag_Pol_Y        ( OC_dom, c1, c2, l,   YPol(l) )
           call tim_mesh_mag_Pol_Y        ( IC_dom, c2, c2, l, icYPol(l) )
         else
           call tim_mesh_Y        ( OC_dom, c1, c2, l,   YPol(l) )
           call tim_mesh_Y        ( IC_dom, c1, c2, l, icYPol(l) )
         end if 
      end do

   end subroutine mag_matrices

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!> Magentic field toroidal boundary conditions (BCs)
   subroutine mag_bc_Tor(OC_dom, IC_dom, A,l)
      use parameters, only : i_Kl, b_diff_profile, b_boussinesq, i_Nic
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(rdom), intent(in) :: IC_dom !< inner core domain
      type (iclumesh), intent(out) :: A !< mesh to apply BCs to
      integer,         intent(in)  :: l !< harmonic degree to apply BCs to
      integer :: j,j_, N_, Ni, No
      
      Ni  = IC_dom%N
      No  = OC_dom%N
      N_  = Ni+No-1
                                ! Insulating interior, BTor = 0
      A%M(2*i_KL+1,1)   = 1d0
                                ! insulating exterior, BTor = 0
      A%M(2*i_KL+1,N_) = 1d0
         			! Jump in dr BTor on icb for cic
      if(Ni==1) return
      if (.not.b_diff_profile.and.b_boussinesq) then 
        do j = 2, 1+i_KL
           j_ = j+Ni-1
           A%M(2*i_KL+1+Ni-j_,j_) =   OC_dom%dr(1)%M(i_KL+1+ 1-j,j)
        end do
        do j = Ni-i_KL, Ni-1
           A%M(2*i_KL+1+Ni-j,j) = - IC_dom%dr(1)%M(i_KL+1+Ni-j,j)
        end do
        A%M(2*i_KL+1,Ni) = OC_dom%dr(1)%M(i_KL+1, 1)  &
                       - IC_dom%dr(1)%M(i_KL+1,Ni)
!
!   if b_diff_profile is .true. or the anelastic reference state is set,
!   and there is an inner core, dT/dr is only continuous across the ICB if the
!   diffusivities are the same on both sides of the ICB and there is no-slip at
!   the ICB.  If either of these is not true, then the boundary condition
!   eta_oc 1/r d/dr rT_oc - eta_ic 1/r  d/dr rT_ic =   (r/l(l+1)) rhat dot curl E is
!   used.
!   Here eta_oc is diffusivity on the outer core side of the ICB, eta_ic the
!   diffusivity on the inner core side: not necessarily the same. T itself is always 
!   continuous across ICB. curl E is non-zero if there is a velocity jump across
!   the ICB. E_theta= u_theta B_r andi E_phi = (u_phi- omega r sin theta) B_r.
!   These are calculated in nonlinear.F90, non-magnetic, stored in mag_bcRe and
!   magbcIm for use in mag_set_bcTor.
!
      else
         do j = 2, 1+i_KL
           j_ = j+Ni-1
           A%M(2*i_KL+1+Ni-j_,j_) =   OC_dom%eta(1,0)*OC_dom%dr(1)%M(i_KL+1+ 1-j,j)
         end do                      
         A%M(2*i_KL+1,Ni) =  OC_dom%eta(1,0)*OC_dom%dr(1)%M(i_KL+1,1)+ &
                               OC_dom%r(1,-1)*OC_dom%eta(1,0)           
         do j = Ni-i_KL, Ni-1
           A%M(2*i_KL+1+Ni-j,j) = - IC_dom%eta(i_Nic,0)* IC_dom%dr(1)%M(i_KL+1+Ni-j,j)
         end do
         A%M(2*i_KL+1,Ni) = A%M(2*i_KL+1,Ni) -IC_dom%eta(i_Nic,0)*IC_dom%dr(1)%M(i_KL+1,Ni)  & 
                                    - IC_dom%r(i_Nic,-1)*IC_dom%eta(i_Nic,0)
      endif

   end subroutine mag_bc_Tor
   
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!> Magnetic field poloidal boundary conditions
   subroutine mag_bc_Pol(OC_dom, IC_dom, A,l)
      use parameters, only : b_inhom_magPol, b_fixed_at_outer, b_fixed_at_inner,&
                             & i_Kl
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(rdom), intent(in) :: IC_dom !< inner core domain
      type (iclumesh), intent(out) :: A !< mesh to apply BCs to
      integer,         intent(in)  :: l !< harmonic degree to apply BCs to
      integer :: j,j_, N_, Ni, No
      
      Ni  = IC_dom%N
      No  = OC_dom%N
      N_  = Ni+No-1
      if (b_inhom_magPol) then
        if (b_fixed_at_outer.and.b_fixed_at_inner) then
          A%M(2*i_KL+1,N_) = 1d0
          A%M(2*i_KL+1,1) = 1d0
        end if 
        if (b_fixed_at_outer.and.(.not.b_fixed_at_inner)) then
          A%M(2*i_KL+1,N_) = 1d0
          if(Ni==1) then
            do j = 1, 1+i_KL
              A%M(2*i_KL+1+1-j,j) = OC_dom%dr(1)%M(i_KL+1+1-j,j)
            end do
            A%M(2*i_KL+1,1) = A%M(2*i_KL+1,1) - OC_dom%r(1,-1)*l
          else
            do j = 1, 1+i_KL
              A%M(2*i_KL+1+1-j,j) = IC_dom%dr(1)%M(i_KL+1+1-j,j)
            end do
            A%M(2*i_KL+1,1) = A%M(2*i_KL+1,1) - IC_dom%r(1,-1)*l
            do j = 2, 1+i_KL
              j_ = j+Ni-1
              A%M(2*i_KL+1+Ni-j_,j_) =   OC_dom%dr(1)%M(i_KL+1+ 1-j,j)
            end do
            do j = Ni-i_KL, Ni-1
              A%M(2*i_KL+1+Ni-j,j) = - IC_dom%dr(1)%M(i_KL+1+Ni-j,j)
            end do
               A%M(2*i_KL+1,Ni) = OC_dom%dr(1)%M(i_KL+1, 1)  &
                       - IC_dom%dr(1)%M(i_KL+1,Ni)
          
          end if

!         if (mpi_rnk==0) print*,' GOT TO 3'
        end if 
        if ((.not.b_fixed_at_outer).and.b_fixed_at_inner) then  
          do j = No-i_KL, No
            j_ = j+Ni-1 
            A%M(2*i_KL+1+N_-j_,j_) = OC_dom%dr(1)%M(i_KL+1+No-j,j)
          end do
          A%M(2*i_KL+1,N_) = A%M(2*i_KL+1,N_) + OC_dom%r(No,-1)*(l+1)
          A%M(2*i_KL+1,1) = 1d0
        endif
      else
                          !  Homogeneous poloidal field boundary conditions here
                                ! Insulating interior,
                                ! (dr - l/r) BPol = 0
      if(Ni==1) then
         do j = 1, 1+i_KL
            A%M(2*i_KL+1+1-j,j) = OC_dom%dr(1)%M(i_KL+1+1-j,j)
         end do
         A%M(2*i_KL+1,1) = A%M(2*i_KL+1,1) - OC_dom%r(1,-1)*l
      else
         do j = 1, 1+i_KL
            A%M(2*i_KL+1+1-j,j) = IC_dom%dr(1)%M(i_KL+1+1-j,j)
         end do
         A%M(2*i_KL+1,1) = A%M(2*i_KL+1,1) - IC_dom%r(1,-1)*l
      end if
                                ! insulating exterior,
                                ! (dr + (l+1)/r) BPol = 0
      do j = No-i_KL, No
         j_ = j+Ni-1 
         A%M(2*i_KL+1+N_-j_,j_) = OC_dom%dr(1)%M(i_KL+1+No-j,j)
      end do
      A%M(2*i_KL+1,N_) = A%M(2*i_KL+1,N_) + OC_dom%r(No,-1)*(l+1)
         			! dr BPol cts on icb for cic
      if(Ni==1) return
      do j = 2, 1+i_KL
         j_ = j+Ni-1
         A%M(2*i_KL+1+Ni-j_,j_) =   OC_dom%dr(1)%M(i_KL+1+ 1-j,j)
      end do
      do j = Ni-i_KL, Ni-1
         A%M(2*i_KL+1+Ni-j,j) = - IC_dom%dr(1)%M(i_KL+1+Ni-j,j)
      end do
      A%M(2*i_KL+1,Ni) = OC_dom%dr(1)%M(i_KL+1, 1)  &
                       - IC_dom%dr(1)%M(i_KL+1,Ni)
      end if 

   end subroutine mag_bc_Pol


!-------------------------------------------------------------------------
!>  set the RHS for the toroidal magnetic boundary condition
!-------------------------------------------------------------------------
   subroutine mag_setbc_Tor(OC_dom, BC_Re, BC_Im, a)
      use timestep, only : tim_it
      type(rdom), intent(in) :: OC_dom !< outer core domain
      double precision, intent(inout) :: BC_Re(:,0:) !< Real part of boundary conditions
      double precision, intent(inout) :: BC_Im(:,0:) !< imaginary part of boundary conditions
      type (coll), intent(out) :: a !< field to apply BCs to
      integer :: N_
      N_ = OC_dom%N
      a%Re(N_,:) = 0d0
      a%Im(N_,:) = 0d0
      a%Re( 1,:) = BC_Re(1,:)
      a%Im( 1,:) = BC_Im(1,:)
      if(tim_it/=0) then
         a%Re(1,:) = a%Re(1,:) - BC_Re(2,:)
         a%Im(1,:) = a%Im(1,:) - BC_Im(2,:)
      end if
      BC_Re(2,:) = BC_Re(1,:)
      BC_Im(2,:) = BC_Im(1,:)
   end subroutine mag_setbc_Tor

!-------------------------------------------------------------------------
!>  set the RHS for the poloidal magnetic boundary condition
!!  Used to impose a magnetic field at the boundaries: magnetoconvection
!-------------------------------------------------------------------------
   subroutine mag_setbc_Pol(OC_dom, mag_bc_PolRe, mag_bc_PolIm, a)
      type(rdom), intent(in) :: OC_dom !< outer core domain
      double precision, intent(in) :: mag_bc_PolRe(:,0:) !< Real part of boundary conditions
      double precision, intent(in) :: mag_bc_PolIm(:,0:) !< imaginary part of boundary conditions
      type (coll), intent(out) :: a !< field to apply BCs to
      integer :: N_
      N_ = OC_dom%N
! a set to zero here if b_inhom_magPol is false
      a%Re( 1,:) = mag_bc_PolRe(1,:)
      a%Im( 1,:) = mag_bc_PolIm(1,:)
      a%Re(N_,:) = mag_bc_PolRe(2,:)
      a%Im(N_,:) = mag_bc_PolIm(2,:)

   end subroutine mag_setbc_Pol


!------------------------------------------------------------------------
!>  get physical B and curl B from toroidal, poloidal
!!
!! calculate B(r), B(th), B(phi),
!! curl_B(r), curl_B(th), and curl_B(phi) in physical space from
!! polidal and toroidal components of B in spectral (coll) space.
!!___
!------------------------------------------------------------------------
   subroutine mag_transform(tor, pol, rad, theta, phi,&
                         curl_rad, curl_theta, curl_phi)
      use transform, only : tra_qst2rtp
      use variables, only : phys, var_coll2spec, var_coll_copy, var_coll_TorPol2qst_mag,&
                            & var_coll_qstcurl
      type(coll), intent(in) :: tor !< toriodal part of field
      type(coll), intent(in) :: pol !< poloidal part of field
      type(phys), intent(inout) :: rad !< radial component of field
      type(phys), intent(inout) :: theta !< theta component of field
      type(phys), intent(inout) :: phi !< phi component of field
      type(phys), intent(inout) :: curl_rad !< radial component of curled field
      type(phys), intent(inout) :: curl_theta !< theta component of curled field
      type(phys), intent(inout) :: curl_phi !< phi component of curled field
      !_loop_lm_vars

      call var_coll_copy(tor, Tor_p_m)
      call var_coll_copy(pol, Pol_p_m)
      call var_coll_TorPol2qst_mag(Tor_p_m,Pol_p_m, cq,cs,ct)
      call var_coll2spec(cq,sq, c2=cs,s2=ss, c3=ct,s3=st) 
      call tra_qst2rtp(sq,ss,st, rad,theta,phi)  
      call var_coll_qstcurl(cq,cs,ct, cq,cs,ct)
      call var_coll2spec(cq,sq, c2=cs,s2=ss, c3=ct,s3=st)
      call tra_qst2rtp(sq,ss,st, curl_rad,curl_theta,curl_phi)

   end subroutine mag_transform


!------------------------------------------------------------------------
!> Predictor for magnetic field
!!
!!  - N  := mN,                         save N at time t
!!  - B  := Y mB + mN,   B := invX B,   get prediction B*
!!  - mB := B                           copy prediction
!------------------------------------------------------------------------
   subroutine mag_predictor(OC_dom, IC_dom, mag_NTor, mag_NPol, mag_icNTor, mag_icNPol,&
                   mag_bcRe, mag_bcIm, mag_bc_PolRe, mag_bc_PolIm,&
                   mag_Tor, mag_Pol, mag_icTor, mag_icPol)
      use variables, only : var_coll_copy
      use timestep, only : tim_multY, tim_zerobc, tim_icinvX
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(rdom), intent(in) :: IC_dom !< inner core domain
      type(coll), intent(in) :: mag_NTor !< toroidal field nonlinear term
      type(coll), intent(in) :: mag_NPol !< poloidal field nonlinear term
      type(coll), intent(in) :: mag_icNTor !< toroidal field nonlinear term (inner core)
      type(coll), intent(in) :: mag_icNPol !< poloidal field nonlinear term (inner core)
      double precision, intent(inout) :: mag_bcRe(:,0:) !< real part toroidal BCs 
      double precision, intent(inout) :: mag_bcIm(:,0:) !< imaginary part toroidal BCs
      double precision, intent(in) :: mag_bc_PolRe(:,0:) !< real part poloidal BCs
      double precision, intent(in) :: mag_bc_PolIm(:,0:) !< imaginary part poloidal BCs
      type(coll), intent(inout) :: mag_Tor !< toroidal field
      type(coll), intent(inout) :: mag_Pol !< poloidal field
      type(coll), intent(inout) :: mag_icTor !< toroidal field (inner core)
      type(coll), intent(inout) :: mag_icPol !< poloidal field (inner core)

      call var_coll_copy (mag_NTor, NTor)
      call tim_multY     (.false.,YTor,mag_Tor,mag_NTor, Tor_p_m)
      call var_coll_copy (mag_icNTor, icNTor)
      call tim_multY     (.false.,icYTor,mag_icTor,mag_icNTor, icTor)
      call tim_zerobc    (icTor)
      call mag_setbc_Tor (OC_dom, mag_bcRe, mag_bcIm,Tor_p_m)
      call tim_icinvX    (OC_dom,IC_dom,.false.,XTor, Tor_p_m,icTor)
      call var_coll_copy (Tor_p_m, mag_Tor)
      call var_coll_copy (icTor, mag_icTor)

      call var_coll_copy (mag_NPol, NPol)
      call tim_multY     (.false.,YPol,mag_Pol,mag_NPol, Pol_p_m)
      call var_coll_copy (mag_icNPol, icNPol)
      call tim_multY     (.false.,icYPol,mag_icPol,mag_icNPol, icPol)
      call tim_zerobc    (icPol)
!     call tim_zerobc    (Pol_p_m)
      call mag_setbc_Pol (OC_dom, mag_bc_PolRe, mag_bc_PolIm, Pol_p_m)
      call tim_icinvX    (OC_dom,IC_dom, .false.,XPol, Pol_p_m,icPol)
      call var_coll_copy (Pol_p_m, mag_Pol)
      call var_coll_copy (icPol, mag_icPol)

   end subroutine mag_predictor
   
   
!------------------------------------------------------------------------
!> Corrector for magnetic field
!!
!!  - B  := c (mN - N),   	using N* get nlin correction
!!  - N  := mN,
!!  - B  := invX B,   	get correction to B
!!  - mB := mB + B		update correction
!------------------------------------------------------------------------
   subroutine mag_corrector(OC_dom, IC_dom, mag_ntor, mag_npol, mag_icNTor, mag_icNPol,&
                       mag_bcRe, mag_bcIm,&
                       mag_tor, mag_pol, mag_icTor, mag_icPol)
      use timestep, only : tim_nlincorr, tim_zerobc, tim_icinvX, tim_addcorr
      use variables, only : var_coll_copy
      type(rdom), intent(in) :: OC_dom !< outer core domain
      type(rdom), intent(in) :: IC_dom !< inner core domain
      type(coll), intent(in) :: mag_NTor !< toroidal field nonlinear term 
      type(coll), intent(in) :: mag_NPol !< poloidal field nonlinear term
      type(coll), intent(in) :: mag_icNTor !< toroidal field nonlinear term (inner core)
      type(coll), intent(in) :: mag_icNPol !< poloidal field nonlinear term (inner core)
      double precision, intent(inout) :: mag_bcRe(:,0:) !< real part toroidal BCs 
      double precision, intent(inout) :: mag_bcIm(:,0:) !< imaginary part toroidal BCs
      type(coll), intent(inout) :: mag_Tor !< toroidal field
      type(coll), intent(inout) :: mag_Pol !< poloidal field
      type(coll), intent(inout) :: mag_icTor !< toroidal field (inner core)
      type(coll), intent(inout) :: mag_icPol !< poloidal field (inner core)


      call tim_nlincorr (mag_NTor,NTor, Tor_p_m)
      call var_coll_copy(mag_NTor,NTor)
      call tim_nlincorr (mag_icNTor,icNTor, icTor)
      call var_coll_copy(mag_icNTor,icNTor)
      call tim_zerobc   (icTor)
      call mag_setbc_Tor(OC_dom, mag_bcRe, mag_bcIm, Tor_p_m)
      call tim_icinvX   (OC_dom,IC_dom,.false.,XTor, Tor_p_m, icTor)
      call tim_addcorr  (Tor_p_m, mag_Tor)
      call tim_addcorr  (icTor, mag_icTor)

      call tim_nlincorr (mag_NPol,NPol, Pol_p_m)
      call var_coll_copy(mag_NPol,NPol)
      call tim_nlincorr (mag_icNPol,icNPol, icPol)
      call var_coll_copy(mag_icNPol,icNPol)
      call tim_zerobc   (icPol)
      call tim_zerobc   (Pol_p_m)
      call tim_icinvX   (OC_dom,IC_dom,.false.,XPol, Pol_p_m, icPol)
      call tim_addcorr  (Pol_p_m, mag_Pol)
      call tim_addcorr  (icPol, mag_icPol)

   end subroutine mag_corrector
   
   
!*************************************************************************
 end module magnetic
!*************************************************************************
