 
!*************************************************************************
! parameters
!***************************************************************************
!> Parameters
!!
!!  - N		n in [1,N] radial r_n
!!  - L		l in [0,L) theta
!!  - M		m in [0,M) phi
!!  - Mp	index m=0,1,2,.. corresponds to m=0,Mp,2Mp,... so Mp imposes Mp-fold symmetry in the phi direction
!***************************************************************************
#include "../parallel.h"
 module parameters

!***************************************************************************
   implicit none
   save

   integer,     parameter :: i_N           = 64 !< Number of outer core radial points
   integer,     parameter :: i_Nic         = 1   !< Number of inner core radial points (for insulating set `=1`)
   integer,     parameter :: i_L           = 64 !< Maximum degree \f$l \in [0,L) \f$
   integer,     parameter :: i_M           = 64 !< Maximum order index, \f$m \in [0,M)\f$
   integer,     parameter :: i_Mp          = 1   !< Azimuthal periodicity i.e. \f$ m=0,M_p,2M_p,...,(M-1)M_p\f$

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!  Input parameters based on midpoint values. See anelastic notes for precise
!  definitions.  d_Hyp is the hyperviscosity parameter, zero being no
!  hyperviscosity.

   !> Modified Rayleigh Number
   double precision, parameter :: d_Ra          = 6.5d04
   !> Prandtl Number
   double precision, parameter :: d_Pr          = 1.0d0
   !> magnetic prandtl number
   double precision, parameter :: d_Pm          = 1.0d0
   !> Ekman Number
   double precision, parameter :: d_E           = 1.0d-03
   !> Ratio of inner/ outer radius
   double precision, parameter :: d_rratio      = 0.35d0
   !> Hyperviscosity parameter
   double precision, parameter :: d_Hyp         = 0.0d-02
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!> Toggle for boussinesq/Anelastic
!!
!!   - For boussinesq: b_boussinesq = .true. (slightly faster)
!!   - For anelastic:  b_boussinesq = .false.
!!
!!  If Anelastic: The file planet_model.txt must be in the run directory.
!!  Matlab scripts are available to generate
!!  planet_model.txt for Jupiter-like models or for polytropic models.
   logical,          parameter :: b_boussinesq  = .true.
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!> Smart Hyperdiffusivity
!!
!!  Smart Hyperdiffusivity, similar to Schaeffer 2017. Consists of a boolean switch
!!  to turn smart hyperdiffusivity on or off, a cut off degree value ldiff: for spherical
!!  harmonic degree > ldiff, hyperdiffusion kicks in. Then the amount of hyperdiffusivity
!!  qhyper, which should be in region 1.001 <= qhyper <= 1.5, higher numbers mean higher
!!  degrees are cut off more quickly
!!  Different to d_Hyp hyperdiffusivity, and only implemented for boussinesq simulations.
   logical,          parameter :: smart_hyperdiffusion = .false.
   !> degree at which smart hyperdiffusion becomes effective
   integer,          parameter :: ldiff         = 30
   !> base rate of hyperdiffusion
   double precision, parameter :: qhyper        = 1.0325
!
!  In the b_boussinesq=.true. case there is an option to read in electrical diffusivity profiles
!  b_diff_profile = .true. for outer core diffusivity and its radial derivative to be read in
!  from the file diffusivity.txt.in. This  must contain i_N entry pairs for the outer
!  core, starting at the inner boundary and going outward, followed by i_Nic
!  diffusivity entry pairs  starting at the centre and finishing at the inner core
!  boundary. If b_diff_profile = .false. no diffusivity.txt.in is needed.
!
   !> Toggle to read in electrical diffusivity profiles from diffusivity.txt.in
   logical,          parameter :: b_diff_profile = .false.
   !> b_rot_step normally false is no inner core, i.e. i_Nic=1
   logical,          parameter :: b_rot_tstep   = .false.    !.true.
   !> Timestep the velocity field
   logical,          parameter :: b_vel_tstep   = .true.
   !> Timestep the codensity
   logical,          parameter :: b_cod_tstep   = .true.
   !> Timestep magnetic field. .false. sets magnetic field to zero throughout run.
   logical,          parameter :: b_mag_tstep   = .false.
   !> Initialise with noise
   logical,          parameter :: b_init_with_noise = .false.
!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

   !> Velocity boundary conditions
   !!
   !! - `=1` is both boundaries no-slip.
   !! - `=2` is no-slip on inner boundary, stress-free on outer boundary
   !! - `=3` is no-slip on the outer boundary stress-free on the inner boundary
   !! - `=4` is stress-free on both boundaries
   !! If both boundaries are stress-free, angular momentum is set to zero.
   integer,          parameter :: i_vel_bc      = 1

   !> Homoegeneous codensity (entropy) boundary conditions
   !!
   !! Only takes effect if b_inhomobc_C = .false.
   !! - `=1` is codensity (entropy) set to 1 on the inner boundary, set to zero on the outer boundary.
   !! - `=2` codensity fixed at 1 on inner boundary, d C/dr is fixed to d_qo at top boundary.
   !! - `=3` codensity fixed to 0 at top boundary and dC/dr fixed to d_qi at the inner bottom boundary.
   !! - `=4` dC/dr fixed to 0 for all l except l=0, and C is set to zero for l=0,  at the top. dC/dr set to  d_qi at bottom.
   integer,          parameter :: i_cod_bc      = 1

   !> Inhomogeneous Boundary conditions for codensity (entropy)
   !!
   !! If .true. then the file codBC.txt.in must be present in the run directory, and it must
   !!  have the format of three integers, the first being 1 for the inner core
   !!  boundary, 2 for the outer core boundary, the next two integers being l and m, the l
   !!  and m harmonic numbers, followed by two reals, the real and imaginary parts
   !!  of the spectral harmonic of cod_C or d C / dr depending on the value of
   !!  i_cod_bc. Currently diffusivity must be continuous across ICB.
   logical,          parameter :: b_inhombc_C = .false.

   !> Codensity (entropy) source term
   double precision, parameter :: d_source = 0.0d0
   !> Codensity (entropy) flux, outer core
   double precision, parameter :: d_qo  = -1.0d0
   !> Codensity (entropy) flux, inner core
   double precision, parameter :: d_qi  = 0.d0

   !> toggle for magnetoconvetion
   !!
   !! If .true. poloidal magnetic field set to mag_bc_PolRe and mag_bc_PolIm on the boundaries r_inner annd r_outer.
   !! These are read in from mag_PolBC.txt.in which must be in the run directory. \n
   !! Also, at least one (or both) of  b_fixed_at_outer, b_fixed_at_inner has to be .true.
   logical,          parameter :: b_inhom_magPol = .false.

   !> Magnetic field fixed at outer core
   !!
   !! If b_inhom_magPol = .true., either b_fixed_at_outer or b_fixed_at_inner
   !! must be true.
   logical,          parameter :: b_fixed_at_outer = .true.
   !> Magnetic field fixed at inner core
   !!
   !! If b_inhom_magPol = .true., either b_fixed_at_outer or b_fixed_at_inner
   !! must be true.
   logical,          parameter :: b_fixed_at_inner = .false.

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   !> Save frequency for snapshot output files
   integer,          parameter :: i_save_rate1  = 1000
   !> Save frequency for temporal output files
   integer,          parameter :: i_save_rate2  = 20
   !> Maximum number of timesteps (no limit if set <0)
   integer,          parameter :: i_maxtstep    = 1500000
   !> Maximum number of cpu hours
   double precision, parameter :: d_cpuhours    = 1d0
   !> Start time (taken from `state.cdf.in` if set < 0)
   double precision, parameter :: d_time        = 0.d0
   !> Fixed timestep (dynamically controlled if set < 0)
   double precision, parameter :: d_timestep    = -5d-05
   double precision, parameter :: d_corr_upper  = 0.020d0   !0.02d0
   double precision, parameter :: d_corr_lower  = 0.00200d0 !0.002d0
   !> Implicitness \f$c\f$
   double precision, parameter :: d_implicit    = 0.5d0
   double precision, parameter :: d_init_pert   = 1d-08 !1d-4 !*******

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   !> Number of adjacent points (finite difference stencil)
   integer,          parameter :: i_KL  = 3
   integer,          parameter :: i_L1  = i_L-1  !< Max degree -1,  (L-1)
   integer,          parameter :: i_M1  = i_M-1  !< Max order - 1,  (m-1)
   integer,          parameter :: i_Ma  = i_M*3/2 !< (3/2) x max order, (3xi_M/2)
   !> Total number of Harmonics, LM - MpM(M-1)/2 -1,  i_L*i_M - i_Mp*(i_M-1)*i_M/2
   integer,          parameter :: i_H   = i_L*i_M - i_Mp*i_M1*i_M/2
   integer,          parameter :: i_Th  = i_L*3/2 !< (3/2) x max degree, (3xi_L/2)
   integer,          parameter :: i_Ph  = i_Ma*2 !< 3 x max order, (3xi_M/2)x2
   integer,          parameter :: i_H1  = i_H/_Nth - 1 !< i_H/_Nth - 1
   integer,          parameter :: i_pH1 = (_Nr+i_H1)/_Nr-1 !< (_Nr + (i_H/_Nth -1) )/_Nr -1
   integer,          parameter :: i_pN  = (_Nr+i_N-1)/_Nr !< (_Nr+i_N-1)/_Nr
   integer,          parameter :: i_pTh = (_Nth+i_Th-1)/_Nth !< (_Nth+i_Th-1)/_Nth
   double precision, parameter :: d_PI  = 3.1415926535897931d0

 contains

!---------------------------------------------------------------------------
!> Check Parameters
!!
!! Checks that @link i_L @endlink, @link i_M @endlink, and @link _Nth @endlink
!! are viable
!---------------------------------------------------------------------------
   subroutine par_precompute()
      if(i_L<=i_M1*i_Mp) stop 'set L > (M-1)*Mp'
      if(i_M>1 .and. modulo(i_M,2)==1) stop 'if(M>1) set M even'
      if(modulo((i_M+1)/2,_Nth)/=0) stop 'require _Nth divides (M+1)/2'
      if(i_N<_Nr) stop 'set i_N>=_Nr'
   end subroutine par_precompute


!***************************************************************************
 end module parameters
!***************************************************************************


