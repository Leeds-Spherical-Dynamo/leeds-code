#include "../parallel.h"
!*************************************************************************
!> Time-stepping routines and variables.
!!   X f* = Y f + N
!!*************************************************************************
 module timestep
!*************************************************************************
   
   implicit none
   save

   double precision :: tim_t !< Instantaneous time
   double precision :: tim_dt !< Instantaneous time step
   double precision :: tim_dterr
   integer          :: tim_it
   integer          :: tim_step
   double precision :: tim_corr_dt
   double precision :: tim_cfl_dt
   logical          :: tim_new_dt
   integer          :: tim_pcn

 contains

!------------------------------------------------------------------------
!>  set the initial time and timestep, 
!!  possibly overwritten by loading a previously saved state
!------------------------------------------------------------------------
   subroutine tim_precompute()

      use parameters, only : d_time, d_timestep

      tim_t       = d_time
      tim_dt      = d_timestep
      tim_dterr   = 0d0
      tim_it      = 0
      tim_step    = 0
      tim_corr_dt = 0d0
      tim_cfl_dt  = 0d0
      tim_new_dt  = .false.
      tim_pcn     = -1
   end subroutine tim_precompute


!------------------------------------------------------------------------
!>  Form the velocity timestepping matrix  (c1/dt - c2 impl Laplace), for each l
!------------------------------------------------------------------------
   subroutine tim_lumesh_vel_X(D,c1,c2,l,A)

      use meshs, only : rdom, lumesh
      use parameters, only : d_implicit, i_Kl, d_Hyp

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (lumesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_, ll1
      integer :: j

      c1_ =  c1 / tim_dt
      c2_ =  c2 * d_implicit
      if(d_implicit==0d0) c2_=1d0
      ll1 = dble(l*(l+1))


!
!   Toroidal/Poloidal velocity equations are
!   dT /dt - d^2 T / dr^2 - (2 xi + 2/r)dT /dr - (2 xi / r + xi^2 + dxi /dr)T
!   + L^2 T /r^2  + d_Hyp*L^4 T /r^4 = R.H.S
!

      A%M(i_KL+1:, :) = - c2_ * D%CompLap%M

!
!  Add in horizontal parts and remaining compressible terms
!
      A%M(2*i_KL+1,:) = A%M(2*i_KL+1,:) + c1_ + c2_*(ll1*D%r(:,-2) &
        + d_Hyp*(ll1*D%r(:,-2))**2 &
        - 2.0d0*D%rho(:,1)*D%r(:,-1) - (D%rho(:,1))**2 - D%rho(:,2))


      do j = 1, 1+i_KL
         A%M(2*i_KL+1+1-j,j) = 0d0
      end do
      do j = D%N-i_KL, D%N
         A%M(2*i_KL+1+D%N-j,j) = 0d0
      end do

   end subroutine tim_lumesh_vel_X

!------------------------------------------------------------------------
!>  Form the Greens function timestepping matrix  
!!  (c1/dt - c2 impl Laplace), for each l
!------------------------------------------------------------------------
   subroutine tim_lumesh_Gre_X(D,c1,c2,l,A)

      use meshs, only : rdom, lumesh
      use parameters, only : d_implicit, i_Kl

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (lumesh),    intent(out) :: A !< greens timestepping matrix
      double precision :: c1_, c2_, ll1
      integer :: j

      c1_ =  c1 / tim_dt
      c2_ =  c2 * d_implicit
      if(d_implicit==0d0) c2_=1d0
      ll1 = dble(l*(l+1))

!
!  Add in horizontal parts and remaining compressible terms
!  Green's function equation is
!  -  d^2 g / dr^2 - ( xi + 2/r)dg /dr - ( xi / r + dxi /dr)g
!   + L^2 g /r^2 = R.H.S
!

      A%M(i_KL+1:, :) = - c2_ * D%GreLap%M

      A%M(2*i_KL+1,:) = A%M(2*i_KL+1,:) + c1_ + c2_*(ll1*D%r(:,-2) &
        - D%rho(:,1)*D%r(:,-1) - D%rho(:,2))

      do j = 1, 1+i_KL
         A%M(2*i_KL+1+1-j,j) = 0d0
      end do
      do j = D%N-i_KL, D%N
         A%M(2*i_KL+1+D%N-j,j) = 0d0
      end do

   end subroutine tim_lumesh_Gre_X

!------------------------------------------------------------------------
!>  Form the entropy timestepping matrix  
!!  (c1/dt - c2 impl Laplace), for each l
!------------------------------------------------------------------------
   subroutine tim_lumesh_X(D,c1,c2,l,A)
      use meshs, only : rdom, lumesh
      use parameters, only : d_implicit, i_Kl, d_Hyp

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (lumesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_, ll1
      integer :: j

      c1_ =  c1 / tim_dt
      c2_ =  c2 * d_implicit
      if(d_implicit==0d0) c2_=1d0
      ll1 = dble(l*(l+1))

!
!  Note: Radial derivative part of entropy diffusion same as
!  for Greens function so EntLap=GreLap
!
!  Entropy  equation is
!  ## Piotr ## Pr -> 1/q
!  q^(-1) dS/dt -  d^2 S / dr^2 - ( xi + 2/r)dS /dr + L^2 S /r^2
!     - d_Hyp L^4 S / r^4 = R.H.S
!

      A%M(i_KL+1:, :) = - c2_ * D%GreLap%M
      A%M(2*i_KL+1,:) = A%M(2*i_KL+1,:) + c1_ + c2_*(ll1*D%r(:,-2) &
                      + d_Hyp*(ll1*D%r(:,-2))**2 )

      do j = 1, 1+i_KL
         A%M(2*i_KL+1+1-j,j) = 0d0
      end do
      do j = D%N-i_KL, D%N
         A%M(2*i_KL+1+D%N-j,j) = 0d0
      end do

   end subroutine tim_lumesh_X

!------------------------------------------------------------------------
!>  Form the flayer velocity  timestepping matrix  
!!  (c1/dt - c2 impl Laplace), for each l
!------------------------------------------------------------------------
   subroutine tim_lumesh_flayer_X(D,c1,c2,l,A)
      use meshs, only : rdom, lumesh
      use parameters, only : d_implicit, i_Kl, d_Hyp

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (lumesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_, ll1
      integer :: j

      c1_ =  c1 / tim_dt
      c2_ =  c2 * d_implicit
      if(d_implicit==0d0) c2_=1d0
      ll1 = dble(l*(l+1))


      A%M(i_KL+1:, :) = - c2_ * D%radLap_flayer%M
      A%M(2*i_KL+1,:) = A%M(2*i_KL+1,:) + c1_ + c2_*(ll1*D%r(:,-2) &
                      + d_Hyp*(ll1*D%r(:,-2))**2 ) * D%f_factor(:)

      do j = 1, 1+i_KL
         A%M(2*i_KL+1+1-j,j) = 0d0
      end do
      do j = D%N-i_KL, D%N
         A%M(2*i_KL+1+D%N-j,j) = 0d0
      end do

   end subroutine tim_lumesh_flayer_X



!------------------------------------------------------------------------
!>  Form the timestepping matrix  (c1/dt - c2 impl Laplace), for each l
!------------------------------------------------------------------------
   subroutine tim_iclumesh_Tor_X(mes_oc, mes_ic, c1,c2,l,A)

      use meshs, only : rdom, iclumesh
      use parameters, only : d_implicit, i_Kl, b_boussinesq, d_Hyp
      
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (iclumesh),  intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_, ll1
      integer :: j, N_, Ni, No
      
      Ni  = mes_ic%N
      No  = mes_oc%N
      N_  = Ni+No-1
      c1_ =  c1 / tim_dt
      c2_ =  c2 * d_implicit
      if(d_implicit==0d0) c2_=1d0
      ll1 = dble(l*(l+1))

      A%M = 0d0
      A%M(i_KL+1:, :Ni) = - c2_ * mes_ic%radLap%M(:,:Ni)
      A%M(2*i_KL+1,:Ni) = A%M(2*i_KL+1,:Ni)  &
                        + c1_ + c2_*ll1*mes_ic%r(:Ni,-2)

      if(Ni==1) A%M=0d0
 

      A%M(i_KL+1:, Ni:N_) = A%M(i_KL+1:, Ni:N_)  &
                          - c2_*mes_oc%MagTorLap%M(:,:No)
!
!   MagTorLap is  eta d^2  / dr^2 + (d eta/dr + 2 eta/r)d /dr ;
!   rest is added below
!
      if (.not.b_boussinesq) then
      A%M(2*i_KL+1,Ni:N_) = A%M(2*i_KL+1,Ni:N_)  &
            + c1_ + c2_*(ll1*mes_oc%eta(:No,0)*mes_oc%r(:No,-2) &
          - mes_oc%eta(:No,1)*mes_oc%r(:No,-1)+ d_Hyp*(ll1*mes_oc%r(:,-2))**2)
      else
      A%M(2*i_KL+1,Ni:N_) = A%M(2*i_KL+1,Ni:N_)  &
                          + c1_ + c2_*ll1*mes_oc%r(:No,-2)
      end if
      
      do j = 1, 1+i_KL
         A%M(2*i_KL+1+1-j,j) = 0d0
      end do
      do j = N_-i_KL, N_
         A%M(2*i_KL+1+N_-j,j) = 0d0
      end do
      if(Ni==1) return
      do j = Ni-i_KL, Ni+i_KL
         A%M(2*i_KL+1+Ni-j,j) = 0d0
      end do

   end subroutine tim_iclumesh_Tor_X


!------------------------------------------------------------------------
!>  Form the timestepping matrix  (c1/dt - c2 impl Laplace), for each l
!------------------------------------------------------------------------
   subroutine tim_iclumesh_Pol_X(mes_oc, mes_ic, c1,c2,l,A)
      use meshs, only : rdom, iclumesh
      use parameters, only : d_implicit, i_Kl, b_boussinesq, d_Hyp
      
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (iclumesh),  intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_, ll1
      integer :: j, N_, Ni, No
      
      Ni  = mes_ic%N
      No  = mes_oc%N
      N_  = Ni+No-1
      c1_ =  c1 / tim_dt
      c2_ =  c2 * d_implicit
      if(d_implicit==0d0) c2_=1d0
      ll1 = dble(l*(l+1))

      A%M = 0d0
      A%M(i_KL+1:, :Ni) = - c2_ * mes_ic%radLap%M(:,:Ni)
      A%M(2*i_KL+1,:Ni) = A%M(2*i_KL+1,:Ni)  &
                        + c1_ + c2_*ll1*mes_ic%r(:Ni,-2)


!
!   Poloidal magnetic equations are
!   dT /dt - eta d^2 T / dr^2 - (2 eta/r)dT /dr 
!   + eta L^2 T /r^2  + d_Hyp*L^4 T /r^4 = R.H.S
!

      if(Ni==1) A%M=0d0
 

      A%M(i_KL+1:, Ni:N_) = A%M(i_KL+1:, Ni:N_)  &
                          - c2_*mes_oc%MagPolLap%M(:,:No)

!
!   MagTorLap is  eta d^2  / dr^2 + (d eta/dr + 2 eta/r)d /dr ;
!   rest is added below
!

      if (.not.b_boussinesq) then
      A%M(2*i_KL+1,Ni:N_) = A%M(2*i_KL+1,Ni:N_)  &
           + c1_ + c2_*(ll1*mes_oc%eta(:No,0)*mes_oc%r(:No,-2) &
                    + d_Hyp*(ll1*mes_oc%r(:,-2))**2)
      else
      A%M(2*i_KL+1,Ni:N_) = A%M(2*i_KL+1,Ni:N_)  &
                       + c1_ + c2_*ll1*mes_oc%r(:No,-2)
      end if

      do j = 1, 1+i_KL
         A%M(2*i_KL+1+1-j,j) = 0d0
      end do
      do j = N_-i_KL, N_
         A%M(2*i_KL+1+N_-j,j) = 0d0
      end do
      if(Ni==1) return
      do j = Ni-i_KL, Ni+i_KL
         A%M(2*i_KL+1+Ni-j,j) = 0d0
      end do

   end subroutine tim_iclumesh_Pol_X


!------------------------------------------------------------------------
!>  Form the matrix  c1/dt + c2 (1-impl) Laplace, for each l
!------------------------------------------------------------------------
   subroutine tim_mesh_vel_Y(D,c1,c2,l,A)

      use meshs, only : rdom, mesh
      use parameters, only : d_implicit, i_Kl, d_Hyp

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (mesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_

      c1_ =  c1 / tim_dt
      c2_ =  c2 * (1d0-d_implicit)


!
!  Add in horizontal parts and remaining compressible terms
!
!
!   Toroidal/Poloidal velocity equations are
!   dT /dt - d^2 T / dr^2 - (2 xi + 2/r)dT /dr - (2 xi / r + xi^2 + dxi /dr)T
!   + L^2 T /r^2  - d_Hyp*L^4 T /r^4 = R.H.S
!

      A%M = c2_ * D%CompLap%M

      A%M(i_KL+1,:) = A%M(i_KL+1,:) + c1_ - c2_*(l*(l+1)*D%r(:,-2) &
        + d_Hyp*(l*(l+1)*D%r(:,-2))**2 &
! #CAJ Mar 2012 #  xi(:0) replaced by rho(:,1)
        - 2.0d0*D%rho(:,1)*D%r(:,-1) - (D%rho(:,1))**2 - D%rho(:,2))

   end subroutine tim_mesh_vel_Y


!------------------------------------------------------------------------
!>  Form the matrix  c1/dt + c2 (1-impl) Laplace, for each l
!------------------------------------------------------------------------
   subroutine tim_mesh_Gre_Y(D,c1,c2,l,A)
      use meshs, only : rdom, mesh
      use parameters, only : d_implicit, i_Kl

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (mesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_

      c1_ =  c1 / tim_dt
      c2_ =  c2 * (1d0-d_implicit)

!
!  Add in horizontal parts and remaining compressible terms
!  Green's function equation is
!  -  d^2 g / dr^2 - ( xi + 2/r)dg /dr - ( xi / r + dxi /dr)g
!   + L^2 g /r^2 = R.H.S
!
      A%M = c2_ * D%GreLap%M

      A%M(i_KL+1,:) = A%M(i_KL+1,:) + c1_ - c2_*(l*(l+1)*D%r(:,-2) &
! #CAJ Mar 2012 #  xi(:0) replaced by rho(:,1), xi(:1) replaced by rho(:,2)
        - D%rho(:,1)*D%r(:,-1) - D%rho(:,2))


   end subroutine tim_mesh_Gre_Y


!------------------------------------------------------------------------
!>  Form the matrix  c1/dt + c2 (1-impl) Laplace, for each l
!------------------------------------------------------------------------
   subroutine tim_mesh_Y(D,c1,c2,l,A)
      use meshs, only : rdom, mesh
      use parameters, only : d_implicit, i_Kl, d_Hyp

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (mesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_

      c1_ =  c1 / tim_dt
      c2_ =  c2 * (1d0-d_implicit)

!
!
!  Note: Radial derivative part of entropy diffusion same as
!  for Greens function so EntLap=GreLap
!
!  Entropy  equation is
!  ## Piotr ## Pr ->1/q
!  q^(-1) dS/dt -  d^2 S / dr^2 - ( xi + 2/r)dS /dr + L^2 S /r^2
!     - d_Hyp L^4 S / r^4 = R.H.S
!
      A%M = c2_ * D%GreLap%M
      A%M(i_KL+1,:) = A%M(i_KL+1,:) + c1_ - c2_*(l*(l+1)*D%r(:,-2) &
                            + d_Hyp*(l*(l+1)*D%r(:,-2))**2 )

   end subroutine tim_mesh_Y

!------------------------------------------------------------------------
!>  Form the matrix  c1/dt + c2 (1-impl) Laplace, for each l
!------------------------------------------------------------------------
   subroutine tim_mesh_flayer_Y(D,c1,c2,l,A)
      use meshs, only : rdom, mesh
      use parameters, only : d_implicit, i_Kl, d_Hyp

      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (mesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_

      c1_ =  c1 / tim_dt
      c2_ =  c2 * (1d0-d_implicit)

      A%M = c2_ * D%radLap_flayer%M
      A%M(i_KL+1,:) = A%M(i_KL+1,:) + c1_ - c2_*( l*(l+1)*D%r(:,-2) &
                            + d_Hyp*(l*(l+1)*D%r(:,-2))**2 ) &
                            * D%f_factor(:) ! f_layer enhanced diffsuion

   end subroutine tim_mesh_flayer_Y



!------------------------------------------------------------------------
!>  Form the matrix  c1/dt + c2 (1-impl) Laplace, for each l
!------------------------------------------------------------------------
   subroutine tim_mesh_mag_Tor_Y(mes_oc, D,c1,c2,l,A)

      use meshs, only : rdom, mesh
      use parameters, only : d_implicit, i_Kl, d_Hyp

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (mesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_

      c1_ =  c1 / tim_dt
      c2_ =  c2 * (1d0-d_implicit)

!
!  ## Piotr ## (03/08/2009_
!  Note: Radial derivative part of magnetic diffusion is now NOT the same as
!  for Greens function so MagLap != GreLap
!
!
!   Toroidal magnetic equations are
!   dT /dt - eta d^2 T / dr^2 - (d eta/dr + 2 eta/r)dT /dr - (d eta/dr )T/r
!   + L^2 T /r^2  + d_Hyp*L^4 T /r^4 = R.H.S
!
      A%M = c2_ * D%MagTorLap%M
      A%M(i_KL+1,:) = A%M(i_KL+1,:) + c1_ - c2_*(l*(l+1)*D%eta(:,0)*D%r(:,-2) &
                      - mes_oc%eta(:,1)*mes_oc%r(:,-1) + d_Hyp*(l*(l+1)*D%r(:,-2))**2 )

   end subroutine tim_mesh_mag_Tor_Y


!------------------------------------------------------------------------
!>  Form the matrix  c1/dt + c2 (1-impl) Laplace, for each l
!------------------------------------------------------------------------
   subroutine tim_mesh_mag_Pol_Y(D,c1,c2,l,A)

      use meshs, only : rdom, mesh
      use parameters, only : d_implicit, i_Kl, d_Hyp
      
      type (rdom),      intent(in)  :: D !< inner/outer core domain
      double precision, intent(in)  :: c1 !< left hand side factor
      double precision, intent(in)  :: c2 !< right hand side factor
      integer,          intent(in)  :: l !< spherical harmonic degree
      type (mesh),    intent(out) :: A !< timestepping matrix
      double precision :: c1_, c2_

      c1_ =  c1 / tim_dt
      c2_ =  c2 * (1d0-d_implicit)


!
!  ## Piotr ## (03/08/2009_
!  Note: Radial derivative part of magnetic diffusion is now NOT the same as
!  for Greens function so MagLap != GreLap
!
!
!   Poloidal magnetic equations are
!   dT /dt - eta d^2 T / dr^2 - (2 eta/r)dT /dr 
!   + eta L^2 T /r^2  + d_Hyp*L^4 T /r^4 = R.H.S
!
      A%M = c2_ * D%MagPolLap%M
      A%M(i_KL+1,:) = A%M(i_KL+1,:) + c1_ - c2_*(l*(l+1)*D%eta(:,0)*D%r(:,-2) &
                       + d_Hyp*(l*(l+1)*D%r(:,-2))**2 )

   end subroutine tim_mesh_mag_Pol_Y



!------------------------------------------------------------------------
!>  solve system Ax=y for x, replaces b,  includes l=m=0 mode if(zero)
!!  uses lapack routine dgbtrs()
!!
!!  Arguments:
!!  - 'N'   = form of equations (no transpose)
!!  - b%D%N = Order of A
!!  - i_KL  = Number of sub-diagonals of A
!!  - i_KL  = Number of super-diagonals of A
!!  - 1     = Number of columns in the RHS
!!  - A%M   = Details of LU factorisation of A
!!  - 3*i_KL+1 = leading dimension of A%M
!!  - A%piv = Number of pivot indices
!!  - b%Re  = The RHS
!!  - b%D%N = Dimension of RHS
!------------------------------------------------------------------------
   subroutine tim_invX(zero,A, b)

      use variables, only : coll
      use meshs, only : lumesh
      use parameters, only : i_L1, i_Kl

      logical,       intent(in)    :: zero !< include l=m=0 mode?
      type (coll),   intent(inout) :: b !< y(in), x(out) for Ax=y
      type (lumesh), intent(in)    :: A(0:i_L1) !< Matrix A for Ax=y
      integer :: info
      _loop_lm_vars

      _loop_lm_begin(b)
         if( l==0 .and. .not.zero) cycle
         call dgbtrs('N', b%D%N, i_KL, i_KL, 1, A(l)%M, 3*i_KL+1,  &
                     A(l)%ipiv, b%Re(1,nh), b%D%N, info )
         if(info/=0) stop 'tim_invX.1'
         call dgbtrs('N', b%D%N, i_KL, i_KL, 1, A(l)%M, 3*i_KL+1,  &
                     A(l)%ipiv, b%Im(1,nh), b%D%N, info )
         if(info/=0) stop 'tim_invX.2'
      _loop_lm_end
   
   end subroutine tim_invX


!------------------------------------------------------------------------
!>  solve system Ax=y for x, replaces b,  includes l=m=0 mode if(zero)
!------------------------------------------------------------------------
   subroutine tim_icinvX(mes_oc,mes_ic,zero,A, oc,ic)

      use meshs, only : rdom, iclumesh
      use variables, only : coll
      use parameters, only : i_L1, i_KL

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      logical,         intent(in)    :: zero !< include l=m=0 mode?
      type (coll),     intent(inout) :: oc !< outer core field
      type (coll),     intent(inout) :: ic !< inner core field
      type (iclumesh), intent(in)    :: A(0:i_L1) !< Mesh for Ax=y
      double precision :: b(mes_ic%N+mes_oc%N-1)
      integer :: info, N_,Ni,No
      _loop_lm_vars
      
      Ni = mes_ic%N
      No = mes_oc%N
      N_ = Ni+No-1

      _loop_lm_begin(oc)
         if( l==0 .and. .not.zero) cycle
            
         b(1:Ni-1) = ic%Re(1:Ni-1,nh)
         b(Ni:N_)  = oc%Re(1:No,  nh)
         call dgbtrs('N', N_, i_KL, i_KL, 1, A(l)%M, 3*i_KL+1,  &
                     A(l)%ipiv, b, N_, info )
         if(info/=0) stop 'tim_invX.1'
         ic%Re(1:Ni,nh) = b(1:Ni)
         oc%Re(1:No,nh) = b(Ni:N_)

         b(1:Ni-1) = ic%Im(1:Ni-1,nh)
         b(Ni:N_)  = oc%Im(1:No,  nh)
         call dgbtrs('N', N_, i_KL, i_KL, 1, A(l)%M, 3*i_KL+1,  &
                     A(l)%ipiv, b, N_, info )
         if(info/=0) stop 'tim_invX.2'
         ic%Im(1:Ni,nh) = b(1:Ni)
         oc%Im(1:No,nh) = b(Ni:N_)
      _loop_lm_end

   end subroutine tim_icinvX


!------------------------------------------------------------------------
!>  multiply d = Ab + c,  includes l=m=0 mode if(zero)
!------------------------------------------------------------------------
   subroutine tim_multY(zero,A,b,c, d)

      use variables, only : coll, var_coll_copy
      use meshs, only : mesh
      use parameters, only : i_L1, i_Kl

      logical,     intent(in)  :: zero !< include l=m=0 mode?
      type (coll), intent(in)  :: b !< b, in: d=Ab+c
      type (coll), intent(in)  :: c !< c, in: d=Ab+c
      type (mesh), intent(in)  :: A(0:i_L1) !< A, in: d=Ab+c
      type (coll), intent(out) :: d !< d, in: d=Ab+c
      integer :: j,nl,nr,n
      _loop_lm_vars

      call var_coll_copy(c, d)

      _loop_lm_begin(b)
         if( l==0 .and. .not.zero) cycle
         do j = 1, d%D%N
            nl = max(1,j-i_KL)
            nr = min(j+i_KL,d%D%N)
            do n = nl, nr
               d%Re(n,nh)  &
                  = d%Re(n,nh) + A(l)%M(i_KL+1+n-j,j) * b%Re(j,nh)
               d%Im(n,nh)  &
                  = d%Im(n,nh) + A(l)%M(i_KL+1+n-j,j) * b%Im(j,nh)
            end do
         end do
      _loop_lm_end
 
   end subroutine tim_multY


!-------------------------------------------------------------------------
!>  set the RHS for the boundary condition = 0
!-------------------------------------------------------------------------
   subroutine tim_zerobc(a)
      use variables, only : coll
      type (coll), intent(inout) :: a
      integer :: N_, n
      N_  = a%D%N
      do n = 0, a%H%pH1
         a%Re(  1, n ) = 0d0
         a%Re( N_, n ) = 0d0
         a%Im(  1, n ) = 0d0
         a%Im( N_, n ) = 0d0
      end do
   end subroutine tim_zerobc
   
   

!-------------------------------------------------------------------------
!>  get the nonlinear correction:  nc := c ( n - np )
!!
!! where 'c' is d_implicit from parameters.F90
!-------------------------------------------------------------------------
   subroutine tim_nlincorr(n,np, nc)
      use variables, only : coll
      use parameters, only : d_implicit
      type (coll), intent(in)  :: n !< n, in: nc=c(n-np)
      type (coll), intent(in)  :: np !< np, in: nc=c(n-np)
      type (coll), intent(out) :: nc !< nc, in: nc=c(n-np)
      integer :: r, nhm
      nc%D => n%D
      nc%H => n%H
      r   = n%D%N
      nhm = n%H%pH1
      nc%Re(1:r,0:nhm) = d_implicit * (n%Re(1:r,0:nhm) - np%Re(1:r,0:nhm))
      nc%Im(1:r,0:nhm) = d_implicit * (n%Im(1:r,0:nhm) - np%Im(1:r,0:nhm))
   end subroutine tim_nlincorr
   

!-------------------------------------------------------------------------
!>  add the correction:  a := a + ac
!-------------------------------------------------------------------------
   subroutine tim_addcorr(ac, a)
      use variables, only : coll
      use dyn_mpi
      type (coll), intent(in)    :: ac
      type (coll), intent(inout) :: a
      double precision :: dterr, d
      integer :: nh, n
      dterr = 0d0
      do nh = 0, a%H%pH1
         do n = 1, a%D%N
            dterr = dterr + ac%Re(n,nh)*ac%Re(n,nh)
            dterr = dterr + ac%Im(n,nh)*ac%Im(n,nh)
            a%Re(n,nh) = a%Re(n,nh) + ac%Re(n,nh)
            a%Im(n,nh) = a%Im(n,nh) + ac%Im(n,nh)
         end do
      end do
#ifdef _MPI
      call mpi_allreduce(dterr, d, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      dterr = d
#endif
      tim_dterr = tim_dterr + dsqrt(dterr)
   end subroutine tim_addcorr


!-------------------------------------------------------------------------
!>  check for convergence via the 2-norm of the correction
!!  -- see tim_addcorr()
!-------------------------------------------------------------------------
   subroutine tim_check_cgce()
      use parameters, only : i_time_step_control, i_maxtstep, d_timestep, i_save_rate2, d_dterr, d_courant
      use dyn_mpi, only : mpi_rnk
      double precision, save :: lasterr
      
      if(i_time_step_control==1) then
          if(tim_it==1) then
    !
    !  Note that tim_corr_dt is now defined simply as the correction
    !  change between the predicted and corrected values
    !
             tim_corr_dt = tim_dterr
             lasterr = 1d99
          end if
          if(tim_dt<1d-12 .and. tim_step>30) then
             if(mpi_rnk==0) print*, 'tim_check_cgce: dt --> 0 !!!?'
             tim_step = i_maxtstep-1
             tim_it = 0
          else if(tim_it>10) then
             if(mpi_rnk==0) print*, 'tim_check_cgce: too many its!!!'
             tim_step = i_maxtstep-1
             tim_it = 0
          else if(tim_dterr>lasterr) then
             if(mpi_rnk==0) print*, 'tim_check_cgce: increasing error!!!'
             if(tim_dterr>1d0) tim_step = i_maxtstep-1
    !         if(tim_dterr<2d0*d_dterr) tim_corr_dt = tim_dt/(1d1*d_courant)
             tim_it = 0
          else if(tim_dterr> 1d0) then
             lasterr = tim_dterr
             tim_it = tim_it + 1
          else
             if(mpi_rnk==0 .and. modulo(tim_step,i_save_rate2)==0) then
                if(d_timestep> 0d0) print*,' step=',tim_step,' its=',tim_it
                if(d_timestep<=0d0) print*,' step=',tim_step,' dt=',real(tim_dt), &
                & 'tim_corr_dt',tim_corr_dt
    ! flushing

    ! flushing
                if(d_timestep> 0d0) write(66,*) ' step=',tim_step,' its=',tim_it
                if(d_timestep<=0d0) write(66,*) ' step=',tim_step,' dt=',real(tim_dt), &
                'tim_corr_dt',tim_corr_dt
              CALL flush(66)
             end if
             tim_it = 0
          end if
          tim_dterr = 0d0
      
      else if(i_time_step_control==2) then
          if(tim_it==1) then
             tim_corr_dt = tim_dt * dsqrt( d_dterr/tim_dterr )
             lasterr = 1d99
          end if         

          if(tim_dt<1d-12 .and. tim_step>30) then
             if(mpi_rnk==0) print*, 'tim_check_cgce: dt --> 0 !!!?'
             tim_step = i_maxtstep-1
             tim_it = 0      
          else if(tim_it==tim_pcn) then
             tim_it = 0
             if(mpi_rnk==0 .and. modulo(tim_step,i_save_rate2)==0) then
                if(d_timestep> 0d0) print*,' step=',tim_step
                if(d_timestep<=0d0) print*,' step=',tim_step,' dt=',real(tim_dt)
             end if
          else if(tim_it>10) then
             if(mpi_rnk==0) print*, 'tim_check_cgce: too many its!!!'
             tim_step = i_maxtstep-1
             tim_it = 0
          else if(tim_dterr>lasterr) then
             if(mpi_rnk==0) print*, 'tim_check_cgce: increasing error!!!'
             if(tim_dterr>2d0*d_dterr) tim_step = i_maxtstep-1
             if(tim_dterr<2d0*d_dterr) tim_corr_dt = tim_dt/(1d1*d_courant)
             tim_it = 0
          else if(tim_dterr>d_dterr) then
             lasterr = tim_dterr
             tim_it = tim_it + 1
          else          
             if(mpi_rnk==0 .and. modulo(tim_step,i_save_rate2)==0) then
                if(d_timestep> 0d0) print*,' step=',tim_step,' its=',tim_it
                if(d_timestep<=0d0) print*,' step=',tim_step,' dt=',real(tim_dt)
             end if
             tim_it = 0
          end if
          tim_dterr = 0d0
      end if

      
   end subroutine tim_check_cgce

!-------------------------------------------------------------------------
!> calculate dt
!!
!! stored in tim_dt
!-------------------------------------------------------------------------
   subroutine tim_new_tstep()

      use parameters, only : i_time_step_control, d_corr_lower, d_corr_upper, d_courant

      double precision :: dt
      
      if(i_time_step_control==1) then
          dt = tim_dt   
    !
    !  New timestep controller goes in here
    !
          if ( tim_corr_dt .gt.d_corr_upper) dt = 0.5d0*dt
          if ( tim_corr_dt .lt.d_corr_lower) dt = 1.273974d0*dt
    !      if(tim_corr_dt>0d0)  dt = min(tim_corr_dt, 1d-03)
    !    dt = min(dt,1d-03)
    !      if(tim_cfl_dt >0d0)  dt = min(tim_cfl_dt, dt)
          if(tim_step == 0  )  dt = min(dt/1d03, 1d-08)
    !      dt = dt * d_courant

          tim_new_dt = (dt<tim_dt*0.97d0 .or. dt>tim_dt*1.03d0)
          if(tim_new_dt)  tim_dt = dt

      else if(i_time_step_control==2) then
          dt = tim_dt   
          if(tim_corr_dt>0d0)  dt = min(tim_corr_dt, dt*2d0/d_courant)
          if(tim_cfl_dt >0d0)  dt = min(tim_cfl_dt, dt)
          if(tim_step == 0  )  dt = min(dt/1d3, 1d-10)
          dt = dt * d_courant
          
          tim_new_dt = (dt<tim_dt*0.97d0 .or. dt>tim_dt*1.03d0)
          if(tim_new_dt)  tim_dt = dt
      end if

   end subroutine tim_new_tstep


!*************************************************************************
 end module timestep
!*************************************************************************

