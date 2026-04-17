 
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
   integer,     parameter :: i_M           = 16 !< Maximum order index, \f$m \in [0,M)\f$
   integer,     parameter :: i_Mp          = 4   !< Azimuthal periodicity i.e. \f$ m=0,M_p,2M_p,...,(M-1)M_p\f$
   double precision, parameter :: d_alpha  = -1.0d0 !< alpha_map, parameter to reduce clustering in radial grid 

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!  Input parameters based on midpoint values. See anelastic notes for precise
!  definitions.  d_Hyp is the hyperviscosity parameter, zero being no
!  hyperviscosity.

   !> Rayleigh Number
   double precision, parameter :: d_Ra          = 6.5d04
   !> Compositional Rayleigh Number
   double precision, parameter :: d_Ra_comp     = 1e0
   !> Prandtl Number
   double precision, parameter :: d_Pr          = 1.0d0
   !> Schmidt Number (Compositional Prandtl Number)
   double precision, parameter :: d_Sc          = 1.0d0
   !> magnetic prandtl number
   double precision, parameter :: d_Pm          = 5.0d0
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

   logical, parameter :: smart_hd_vel = .false.
   logical, parameter :: smart_hd_cod = .false.
   logical, parameter :: smart_hd_mag = .false.
   logical, parameter :: smart_hd_comp = .false.

   integer, parameter :: ldiff = 30
   double precision, parameter :: qh = 1.0325

   ! The following are for different amounts
   ! of hyperdiffusion for different quantities.
   ! All fields will use ldiff and qh as default,
   ! change the follwing params to positive values
   ! to specify them separately.
   integer :: ldiff_vel = -1
   double precision :: qh_vel = -1
   integer :: ldiff_cod = -1
   double precision :: qh_cod = -1
   integer :: ldiff_mag = -1
   double precision :: qh_mag = -1
   integer :: ldiff_comp = -1
   double precision :: qh_comp = -1

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
   logical,          parameter :: b_rot_tstep   = .false.
   !> Force no rotation of inner core, still calculate torques
   logical,          parameter :: b_rot_force_0 = .false.
   !> Timestep the velocity field
   logical,          parameter :: b_vel_tstep   = .true.
   !> Timestep the codensity
   logical,          parameter :: b_cod_tstep   = .true.
   !> Timestep the Composition
   logical,          parameter :: b_comp_tstep  = .false.
   !> Timestep magnetic field. .false. sets magnetic field to zero throughout run.
   logical,          parameter :: b_mag_tstep   = .true.
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
   !! - `=1` is codensity (entropy) set to 1 on the inner boundary, set to zero on the outer boundary.
   !! - `=2` codensity fixed at 1 on inner boundary, d C/dr is fixed to d_qo at top boundary.
   !! - `=3` codensity fixed to 0 at top boundary and dC/dr fixed to d_qi at the inner bottom boundary.
   !! - `=4` dC/dr fixed to 0 for all l except l=0, and C is set to zero for l=0,  at the top. dC/dr set to  d_qi at bottom.
   integer,          parameter :: i_cod_bc      = 1
   !> Inhomogeneous Boundary conditions (codensity)
   !! 
   !! 3 options here:
   !! 1: Homogeneous boundary conditions
   !! 2: inhomogeneous boundary conditions, loaded from codBC.cdf.in
   !! 3: inhomogeneous boundary conditions, loaded from codBC.txt.in
   integer,          parameter :: i_inhombc_C = 1
   !> Codensity (entropy) flux, outer core, if i_inhombc_C=1
   double precision, parameter :: d_qo  = -1.0d0
   !> Codensity (entropy) flux, inner core, if i_inhombc_C=1
   double precision, parameter :: d_qi  = 0.d0
   !> Codensity source term
   !!
   !! 1: uniform, load from d_source
   !! 2: load from codBC.cdf.in
   !! 3: load from codS.txt.in
   integer,          parameter :: i_source_load = 1
   !> Codensity (entropy) source term
   !! This is only used if i_source_load = 1 
   double precision, parameter :: d_source = 0.0d0


   
   !> Homoegeneous compositonal boundary conditions
   !!
   !! - `=1` is composition set to 1 on the inner boundary, set to zero on the outer boundary.
   !! - `=2` composition fixed at 1 on inner boundary, d Comp/dr is fixed to d_qo_comp at top boundary.
   !! - `=3` composition fixed to 0 at top boundary and dComp/dr fixed to d_qi_comp at the inner bottom boundary.
   !! - `=4` dcomp/dr fixed to 0 for all l except l=0, and Composition is set to zero for l=0,  at the top. dComp/dr set to  d_qi_comp at bottom.
   integer,          parameter :: i_comp_bc      = 1
   !> Inhomogeneous Boundary conditions (compostion) 
   !! 
   !! 4 options here:
   !! 1: Homogeneous boundary conditions
   !! 2: inhomogeneous boundary conditions, loaded from compBC.cdf.in
   !! 3: inhomogeneous boundary conditions, loaded from compBC.txt.in
   !! 4: flayer
   integer,          parameter :: i_inhombc_comp = 1
   !> Codensity (entropy) flux, outer core
   double precision, parameter :: d_qo_comp  = -1.0d0
   !> Codensity (entropy) flux, inner core
   double precision, parameter :: d_qi_comp  = 0.d0
   !> Codensity source term
   !!
   !! 1: uniform, load from d_source_comp
   !! 2: load from codBC.cdf.in
   !! 3: load from codS.txt.in
   !! 4: flayer
   integer,          parameter :: i_source_comp_load = 1
   !> Codensity (entropy) source term
   double precision, parameter :: d_source_comp = 0.0d0

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
! Flayer parameters

   double precision, parameter :: d_enhance_diff = 1.0 !< enhanced diffusion parameter for f-layer
   double precision, parameter :: d_flayer_radius = -1.0d0 !< ri<x<ro,  r-position of flayer. Set negative to turn off.
   double precision, parameter :: d_transition_width = 0.05 !< smoothing width for hyperdiffusion transition  
   double precision, parameter :: d_fl_q = -2d0 !< flayer icb flux
   double precision, parameter :: d_fl_kl = 0.1d0 !< enhanced diffusivity (lower)
   double precision, parameter :: d_fl_ku = 1d0 !< enhanced diffusivity (upper)
   double precision, parameter :: d_fl_Cps = 2.51952632905013d0 !< composition gradient at rs  

!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
   !> write outputs after x timesteps (1) or y diffusion times (2)
   !!
   !! Be careful using i_save_option = 2, you could easily run
   !! an entire simulation and produce no outputs !!!
   integer,          parameter :: i_save_option = 1
   !> time between snapshots, in simulation time
   double precision, parameter :: d_save_state = 1e-2
   !> time between spectra, profiles, in sim time
   double precision, parameter :: d_save_spec = 1e-2
   !> time between time series writes, in sim time
   double precision, parameter :: d_save_series = 1e-4
   !> Save frequency for snapshot output files
   integer,          parameter :: i_save_rate1  = 2000
   !> Save frequency for temporal output files
   integer,          parameter :: i_save_rate2  = 20
   !> Save frequency for spectra output files
   integer,          parameter :: i_save_rate3  = 1000
   !> Maximum number of timesteps (no limit if set <0)
   integer,          parameter :: i_maxtstep    = 1500000
   !> Maximum number of cpu hours
   double precision, parameter :: d_cpuhours    = 1d0
   !> Start time (taken from `state.cdf.in` if set < 0)
   double precision, parameter :: d_time        = -1.0d0
   !> timestep control, 1=unified version, 2=legacy (Willis) version
   integer,          parameter :: i_time_step_control = 1
   double precision, parameter :: d_courant     = 0.8 !< courant timestep safety factor
   double precision, parameter :: d_dterr       = 0.025d0  !< timestep convergence criteria
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
   integer,          parameter :: i_L1  = i_L-1  !< Max degree, (i_L-1)
   integer,          parameter :: i_M1  = i_M-1  !< Max order, (i_m-1)
   integer,          parameter :: i_Ma  = i_M*3/2 !< (3/2) x max order, (3*i_M/2)
   !> Total number of Harmonics, LM - MpM(M-1)/2 -1,  i_L*i_M - i_Mp*(i_M-1)*i_M/2
   integer,          parameter :: i_H   = i_L*i_M - i_Mp*i_M1*i_M/2
   integer,          parameter :: i_Th  = i_L*3/2 !< global num points in theta, (3*i_L/2)
   integer,          parameter :: i_Th2 = (i_Th+1)/2 !< (global num points in theta)/2
   integer,          parameter :: i_Ph  = i_Ma*2 !< global num points in phi, (3*i_M)
   integer,          parameter :: i_H1  = i_H/_Nth - 1 !< local num harmonics in type(spec)
   integer,          parameter :: i_pH1 = (_Nr+i_H1)/_Nr-1 !< local num harmonics in type(coll)
   integer,          parameter :: i_pN  = (_Nr+i_N-1)/_Nr !< local num radial points in spec, phys
   integer,          parameter :: i_pTh = (_Nth+i_Th-1)/_Nth !< local num theta points in type(phys)
   integer,          parameter :: dp    = selected_real_kind(15)
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
      ! set smart HD params
      if(ldiff_vel<0) ldiff_vel = ldiff
      if(ldiff_mag<0) ldiff_mag = ldiff
      if(ldiff_cod<0) ldiff_cod = ldiff
      if(ldiff_comp<0) ldiff_comp = ldiff
      if(qh_vel<0) qh_vel=qh
      if(qh_mag<0) qh_mag=qh
      if(qh_cod<0) qh_cod=qh
      if(qh_comp<0) qh_comp=qh
   end subroutine par_precompute


!***************************************************************************
 end module parameters
!***************************************************************************


