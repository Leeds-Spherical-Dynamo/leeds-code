#include "../parallel.h"
#include "../version.h"
!**************************************************************************
!> IN/OUT 
!!**************************************************************************
 module io
!**************************************************************************
   use variables, only : coll

   implicit none
   save

   character(12), parameter :: code_version = _version !< code version defined in VERSION.h
   character(200)        :: io_statefile !< filename, statefile
   character(200)        :: io_codBCfile !< filename, codensity BCs, netcdf
   character(200)        :: io_refstatefile !< filename, reference, state
   character(200)        :: io_codBCtxt !< filename, codnsity BCs, txt
   character(200)        :: io_codStxt !< filename, codensity source, txt
   character(200)        :: io_diffusfile !< filename, diffusivity profile, txt
   character(200)        :: io_magB0file !< filename, mag conditions, cdf
   character(200)        :: io_magPolBCfile !< filename, poloidal mag BCs, txt
   character(200)        :: io_compBCfile !< filename, composition BCs, cdf
   character(200)        :: io_compBCtxt !< filename, composition BCs, txt
   character(200)        :: io_compStxt !< filename, composition source, txt
   integer,     private  :: io_save1 !< statefile serial number
   integer,     private  :: io_save2 !< time series serial number
   integer,     private  :: io_save3 !< spectrum/profile file serial number
   double precision, private :: io_savetime1 !< time of next statefile write
   double precision, private :: io_savetime2 !< time of next time series write
   double precision, private :: io_savetime3 !< time of next spectra/profile write
   integer,     private  :: io_vE !< unit identifier, vel_energy.dat
   integer,     private  :: io_mE !< unit identifier, mag_energy.dat
   integer,     private  :: io_mIC !< unit identifier, mag_ic.dat
   integer,     private  :: io_mMB !< unit identifier, mag_cmb.dat
   integer,     private  :: io_mpar !< unit identifier, mag_partition.dat
   integer,     private  :: io_vpar !< unit identifier, vel_partition.dat
   integer,     private  :: io_dynOut !< unit identifier, dyn_outputs.dat
   integer,     private  :: io_cNu !< unit identifier, cod_nusselt.dat
   integer,     private  :: io_rTq !< unit identifier, rot_torque.dat
   integer,     private  :: io_tdt !< unit identifier, tim_step.dat
   integer,     private  :: io_pts !< unit identifier, pts_bench.dat
   integer,     private  :: io_compNu !< unit identifier, comp_nusselt.dat
   !type (coll), private  :: cq,cs,ct
   logical,     private  :: do_random_init
   type (coll), private  :: cq,cs,ct
   double precision, private :: io_KinEn !< total kinetic energy
   double precision, private :: io_MagEn !< total magnetic energy
   double precision, private :: io_Nus !< nusselt number
   double precision, private :: io_ViscDiss !< viscous dissipation
   double precision, private :: io_OhmDiss !< ohmic dissipation
   double precision, private :: io_fdip !< fraction of ohmic dissipation
   double precision, private :: io_thdip !< dipole tilt angle
   double precision, private :: io_Nus_comp !< sherwood number (compositional nusselt number)
   double precision, private :: io_vol !< non-dimensional spherical shell volume

 contains
 
 
!--------------------------------------------------------------------------
!>  initialiser fn
!!
!! Puts file names and unit identifiers into variables and sets
!! all counters to 0.
!--------------------------------------------------------------------------
   subroutine io_precompute(mes_ic)
      use meshs, only : rdom
      use parameters, only : d_rratio, d_Pi, dp, b_cod_tstep, b_boussinesq,&
                             b_comp_tstep, b_vel_tstep, b_vel_tstep, b_mag_tstep,&
                             b_rot_tstep, d_timestep
      
      type(rdom), intent(in) :: mes_ic !< inner core domain

      double precision :: ri, ro
      ri = d_rratio/(1-d_rratio)
      ro = 1/(1-d_rratio)
      io_vol = (4.0_dp/3.0_dp * d_PI) * (ro**3 - ri**3)
      io_statefile = 'state.cdf.in'
      io_codBCfile = 'codBC.cdf.in'
      io_codBCtxt = 'codBC.txt.in'
      io_codStxt = 'codS.txt.in'
      io_compBCfile = 'compBC.cdf.in'
      io_compBCtxt = 'compBC.txt.in'
      io_compStxt = 'compS.txt.in'
      io_magPolBCfile = 'mag_PolBC.txt.in'
      io_refstatefile = 'planet_model.txt.in'
      io_diffusfile = 'diffusivity.txt.in' 
      io_magB0file = 'magB0.cdf.in'
      io_save1 = 0
      io_save2 = 0
      io_save3 = 0
      io_vE  = 0
      io_vpar = 0
      io_mE  = 0
      io_mpar = 0
      io_mIC = 0
      io_mMB = 0
      io_cNu = 0
      io_compNu = 0
      io_rTq = 0
      io_tdt = 0
      io_pts = 0
      if(b_cod_tstep .and.  &
         b_boussinesq) io_cNu  = 50
      if (b_comp_tstep.and.  &
         b_boussinesq) io_compNu = 55
      if(b_vel_tstep)                   io_vE   = 30
      if(b_vel_tstep)                   io_vpar  = 35
      if(b_mag_tstep)                   io_mE   = 40
      if(b_mag_tstep .and. mes_ic%N>1 ) io_mIC  = 41
      if(b_mag_tstep)                   io_mMB  = 42
      if(b_mag_tstep)                   io_mpar = 45
      if(b_rot_tstep)                   io_rTq  = 60
      if(d_timestep<0d0)                io_tdt  = 70
      io_dynOut = 75
   end subroutine io_precompute
 
 

!--------------------------------------------------------------------------
!>  Open files written to every ... steps runtime
!--------------------------------------------------------------------------
   subroutine io_openfiles()
      use dyn_mpi, only : mpi_rnk
      if(mpi_rnk/=0)  return
      if(io_mpar /=0)  open(io_mpar,  status='unknown', file='mag_partition.dat')
      if(io_vpar /=0)  open(io_vpar,  status='unknown', file='vel_partition.dat') 
      if(io_vE /=0)   open(io_vE,  status='unknown', file='vel_energy.dat')
      if(io_mE /=0)   open(io_mE,  status='unknown', file='mag_energy.dat')
      if(io_mIC/=0)   open(io_mIC, status='unknown', file='mag_ic.dat')
      if(io_mMB/=0)   open(io_mMB, status='unknown', file='mag_cmb.dat')
      if(io_cNu/=0)   open(io_cNu, status='unknown', file='cod_nusselt.dat')
      if(io_compNu/=0)   open(io_compNu, status='unknown', file='comp_nusselt.dat')
      if(io_rTq/=0)   open(io_rTq, status='unknown', file='rot_torque.dat')
      if(io_tdt/=0)   open(io_tdt, status='unknown', file='tim_step.dat'  )
      if(io_pts/=0)   open(io_pts, status='unknown', file='pts_bench.dat'  )
      open(io_dynOut, status='unknown', file='dynamo_outputs.dat')
   end subroutine io_openfiles


!--------------------------------------------------------------------------
!>  Open files written to every ... steps runtime
!--------------------------------------------------------------------------
   subroutine io_force_flush_files()
      use dyn_mpi, only : mpi_rnk
      if(mpi_rnk/=0) return
      open(io_dynOUT, status='unknown', file='dynamo_outputs.dat')
      close(io_dynOut)
      if(io_vE /=0) then
          close(io_vE)
          open(io_vE, position='append', status='unknown', file='vel_energy.dat')
          close(io_vpar)
          open(io_vpar, position='append', status='unknown', file='vel_partition.dat')
      endif
      if(io_mE /=0) then
          close(io_mE)
          open(io_mE, position='append', status='unknown', file='mag_energy.dat')
          close(io_mpar)
          open(io_mpar, position='append', status='unknown', file='mag_partition.dat')
      endif
      if(io_mIC/=0) then
          close(io_mIC)
          open(io_mIC, position='append', status='unknown', file='mag_ic.dat')
      endif
      if(io_mMB/=0) then
          close(io_mMB)
          open(io_mMB, position='append', status='unknown', file='mag_cmb.dat')
      endif
      if(io_cNu/=0)  open(io_cNu, status='unknown', file='cod_nusselt.dat')
      if(io_rTq/=0) then
          close(io_rTq)
          open(io_rTq, position='append', status='unknown', file='rot_torque.dat')
      endif
      if(io_tdt/=0) then
          close(io_tdt)
          open(io_tdt, position='append', status='unknown', file='tim_step.dat'  )
      endif
      if(io_pts/=0) then
          close(io_pts)
          open(io_pts, position='append', status='unknown', file='pts_bench.dat'  )
      endif
   end subroutine io_force_flush_files



!--------------------------------------------------------------------------
!>  Close files written to during runtime
!--------------------------------------------------------------------------
   subroutine io_closefiles()
      use dyn_mpi, only : mpi_rnk
      if(mpi_rnk/=0) return
      if(io_vE /=0) close(io_vE)
      if(io_mpar/=0) close(io_mpar)
      if(io_vpar/=0) close(io_vpar)
      if(io_mE /=0) close(io_mE)
      if(io_mIC/=0) close(io_mIC)
      if(io_mMB/=0) close(io_mMB)
      if(io_cNu/=0) close(io_cNu)
      if(io_compNu/=0) close(io_compNu)
      if(io_rTq/=0) close(io_rTq)
      if(io_tdt/=0) close(io_tdt)
      if(io_pts/=0) close(io_pts)
      close(io_dynOut)
   end subroutine io_closefiles


!--------------------------------------------------------------------------
!>  Write to files
!--------------------------------------------------------------------------
   subroutine io_write2files(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol, mag_icTor, mag_icPol,&
                          rot_omega, rot_velTorq, rot_magTorq, cod_C, comp_C)

      use meshs, only : rdom
      use parameters, only : d_Pm, d_E, b_boussinesq, i_save_option, i_save_rate1, i_save_rate2,&
                             i_save_rate3, d_save_state, d_save_state, d_save_series, d_save_spec
      use timestep, only : tim_step, tim_t, tim_new_dt
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      type(coll), intent(in) :: vel_tor !< toroidal component of velocity
      type(coll), intent(in) :: vel_pol !< poloidal component of velocity
      type(coll), intent(in) :: mag_tor !< toroidal component of magnetic field
      type(coll), intent(in) :: mag_pol !< poloidal component of magnetic field
      type(coll), intent(in) :: mag_icTor !< toroidal component of mag field (inner core)
      type(coll), intent(in) :: mag_icPol !< toroidal component of mag field (inner core)
      type(coll), intent(in) :: cod_C !< codensity
      type(coll), intent(in) :: comp_C !< composition
      double precision, intent(in) :: rot_omega !< inner core rotation rate
      double precision, intent(in) :: rot_velTorq !< viscous torque on inner core
      double precision, intent(in) :: rot_magTorq !< magnetic torque on inner core

      double precision :: KE !< pre-factor for kinetic energy
      double precision :: VD !< pre-factor for viscous dissipation
      double precision :: ME !< pre-factor for magnetic energy
      double precision :: OD !< pre-factor for ohmic dissipation
      logical :: write_state, write_spec, write_series

      KE = 0.5d0
      VD = d_Pm
!  ## Piotr ## March 18, 2008
!  d_Ro -> d_E/d_Pm
      ME = 0.5d0*d_Pm/d_E
      OD = 1d0*d_Pm/d_E
      if(b_boussinesq) then
         ME = d_Pm/d_E
         OD = 2d0*d_Pm/d_E
      end if

      write_state=.false.
      write_spec=.false.
      write_series=.false.

      if (i_save_option==1) then 
         ! save based on numer of time steps
         if(modulo(tim_step,i_save_rate1)==0) write_state=.true.
         if(modulo(tim_step,i_save_rate2)==0) write_series=.true.
         if(modulo(tim_step,i_save_rate3)==0) write_spec=.true.
      else if (i_save_option == 2) then 
         ! save based on in-simulation time
         if (tim_step == 0) then 
            io_savetime1 = tim_t 
            io_savetime2 = tim_t 
            io_savetime3 = tim_t
            write_state=.true.
            write_spec=.true.
            write_series=.true.
         else
            if (tim_t > (io_savetime1+d_save_state )) then 
               write_state=.true.
               io_savetime1=tim_t
               ! truncate time
               io_savetime1=dble(int(io_savetime1/d_save_state))*d_save_state
            end if
            if (tim_t>(io_savetime2 + d_save_series)) then 
               write_series=.true.
               io_savetime2=tim_t
               io_savetime2=dble(int(io_savetime2/d_save_series))*d_save_series
            endif
            if (tim_t > (io_savetime3 + d_save_spec)) then 
               write_spec=.true.
               io_savetime3=tim_t
               io_savetime3=dble(int(io_savetime3/d_save_spec))*d_save_spec
            end if
         end if
      end if

      ! Save state
      if(write_state) then
         call io_save_state(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol, mag_icTor, mag_icPol,cod_C, comp_c, rot_omega)
         io_save1 = io_save1+1
      endif

      ! Save spectra and profiles
      if(write_spec) then
         !call io_save_state()
         if (b_boussinesq) then
           call io_save_spectrum('vel',vel_Tor,vel_Pol,KE)
         else 
           call io_save_spectrum_anelastic('vel',vel_Tor,vel_Pol,KE)
         endif
         call io_save_spectrum ('mag',mag_Tor,mag_Pol,ME)
         call io_save_spectrum1('cod',cod_C,1d0)
         call io_save_spectrum1('comp',comp_C,1d0)
         call io_save_profile(mes_oc,'temp',cod_C)
         call io_save_profile(mes_oc,'comp',comp_C)
         call io_rad_profile('vel', vel_Tor, Vel_Pol, KE, mes_oc)
         call io_rad_profile('mag', mag_Tor, mag_Pol, ME, mes_oc)
         io_save3 = io_save3+1
      endif

      ! Save time series
      if(write_series) then
         if(io_vE /=0) call io_write_vel_energy(mes_oc,io_vE, io_vpar,vel_Tor,vel_Pol,cod_C,KE,VD,tim_step, cod_C)
         if(io_mE /=0) call io_write_mag_energy(mes_oc,io_mE, io_mpar, mag_Tor,  mag_Pol,  ME,OD, mag_tor, mag_pol)
         if(io_mIC/=0) call io_write_mag_energy(mes_oc,io_mIC, io_mpar, mag_icTor,mag_icPol,ME,OD, mag_tor, mag_pol)
         if(io_mMB/=0) call io_write_cmb   (mes_oc,io_mMB,mag_Tor,  mag_Pol,  ME)
         if(io_cNu/=0) call io_write_nussel(mes_oc,io_cNu,cod_C,1)
         if(io_compNu/=0) call io_write_nussel(mes_oc,io_compNu,comp_C,2)
         if(io_rTq/=0) call io_write_torque(io_rTq,rot_omega, rot_velTorq, rot_magTorq)
         ! Write non-dimensional standard dynamo outputs eg Re, Rm.
         call io_write_dyn_outputs(io_dynOut)
!        if(io_pts/=0) call io_write_bench (io_pts)
! ## Piotr ## (16/07/2009)
! Flushing the content of vel_energy.dat and mag_energy.dat by closing and reopening these files.
         if( (mpi_rnk .eq. 0) .and. (modulo(tim_step,i_save_rate2*100)==0) ) then
           if(io_vE /=0)  then
             close(io_vE)
             open(io_vE,  position='append', status='unknown', file='vel_energy.dat')
           end if
           if(io_mE /=0)  then
             close(io_mE)
             open(io_mE,  position='append', status='unknown', file='mag_energy.dat')
           end if
         end if


        io_save2 = io_save2+1
      end if

      if(tim_new_dt) call io_write_timestep(io_tdt)
   
   end subroutine io_write2files



!--------------------------------------------------------------------------
!>  Load state - start from previous solution
!--------------------------------------------------------------------------
   subroutine io_load_state(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol, mag_icTor,& 
                    mag_icPol, cod_C, comp_C, rot_omega)

      use netcdf
      use meshs, only : rdom
      use parameters, only : i_KL, d_time, b_mag_tstep, d_timestep
      use dyn_mpi
      use timestep, only : tim_t, tim_dt

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      type(coll), intent(inout) :: vel_tor !< toroidal component of velocity
      type(coll), intent(inout) :: vel_pol !< poloidal component of velocity
      type(coll), intent(inout) :: mag_tor !< toroidal component of magnetic field
      type(coll), intent(inout) :: mag_pol !< poloidal component of magnetic field
      type(coll), intent(inout) :: mag_icTor !< toroidal component of mag field (inner core)
      type(coll), intent(inout) :: mag_icPol !< toroidal component of mag field (inner core)
      type(coll), intent(inout) :: cod_C !< codensity
      type(coll), intent(inout) :: comp_C !< composition
      double precision, intent(out) :: rot_omega !< inner core rotation rate

      integer :: e, e1, e2, f, i, rd, icrd, N_
      integer :: n, No,Ni, No_,Ni_
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: Ao(lda,lda,mes_oc%N)
      double precision :: Ai(lda,lda,mes_ic%N)
      double precision, allocatable :: roc(:), ric(:)
      logical :: into, inti
!for angular momentum removal
      double precision :: w,f1(vel_Tor%D%N), f2(vel_Tor%D%N),Lz,Lz2,omz
      double precision :: Lx,Lx2,omx,Ly,Ly2,omy
      _loop_lm_vars


      e=nf90_open(io_statefile, nf90_nowrite, f)      
      if(e/=nf90_noerr) then
         if(mpi_rnk==0) open(99, file='PRECOMPUTING')
         if(mpi_rnk==0) close(99, status='delete')
         stop 'state file not found!'
      end if
!  ## Piotr ##
!  initialize random numbers generator
!  (to perturb modes not present in the initail state file
!   with a random noise of small amplitude)

      call random_seed
      e=nf90_get_att(f, nf90_global,  't', d)
!      print*,d
      if(d_time<0d0) tim_t = d

      No = mes_oc%N
      Ni = mes_ic%N
      e=nf90_inq_dimid(f,    'r', rd)
      e=nf90_inq_dimid(f,  'icr', icrd)
      e=nf90_inquire_dimension(f,   rd, len=No_)
      e=nf90_inquire_dimension(f, icrd, len=Ni_)

      allocate(roc(No_))
      allocate(ric(Ni_))
      e=nf90_inq_varid(f, 'r', i)
      e=nf90_get_var(f, i, roc)
      e=nf90_inq_varid(f, 'icr', i)
      e=nf90_get_var(f, i, ric)
      into = (No_/=No)
      if(.not.into) into = (maxval(dabs(mes_oc%r(1:No,1)-roc))>1d-8)
      inti = (Ni_/=Ni)
      if(.not.inti) inti = (maxval(dabs(mes_ic%r(1:Ni,1)-ric))>1d-8)
      if(into) then
         if(mpi_rnk==0) print*,' N  :', No_, ' --> ',No 
         call io_interp_wts(No_,roc,No,mes_oc%r(1,1), Ao)
      end if
       if(inti .and. Ni/=1 .and. Ni_/=1) then
         if(mpi_rnk==0) print*,' Nic:', Ni_, ' --> ',Ni 
         call io_interp_wts(Ni_,ric,Ni,mes_ic%r(1,1), Ai)
      end if

      e=nf90_inq_varid(f, 'dt', i)
      e=nf90_get_var(f, i, d)
!      print*,d
      if(d_timestep<0d0) tim_dt = d
!      print*,'TIMESTEPPING:',tim_dt,tim_t
      e=nf90_inq_varid(f, 'icOmega', i)
      e=nf90_get_var(f, i, rot_omega)
      do_random_init = .true.
      call io_load_coll(f,'uT',into,No_,roc,Ao, vel_Tor)
      call io_load_coll(f,'uP',into,No_,roc,Ao, vel_Pol)
      call io_load_coll(f, 'C',into,No_,roc,Ao, cod_C)
      call io_load_coll(f,'BT',into,No_,roc,Ao, mag_Tor)
      call io_load_coll(f,'BP',into,No_,roc,Ao, mag_Pol)
      e=nf90_inq_varid(f, 'icBT', i)
      e1=nf90_inq_varid(f, 'icBP', i)
      if((inti .and. Ni_==1).or. e/=nf90_NoErr .or. e1/=nf90_NoErr) then
         call io_ic_extend(mes_oc,mes_ic,'icBT',mag_Tor,mag_icTor)
         call io_ic_extend(mes_oc,mes_ic,'icBP',mag_Pol,mag_icPol)
      else if(Ni/=1) then
         call io_load_coll(f,'icBT',inti,Ni_,ric,Ai, mag_icTor)
         call io_load_coll(f,'icBP',inti,Ni_,ric,Ai, mag_icPol)      
      end if
      ! check if state file has compostional information, if so, load it.
      ! if not, load codensity field from statefile into composition, as a first
      ! guess at initial conditions.
      e2=nf90_inq_varid(f, 'Comp', i)
      if (e2 == nf90_NoErr) then
          call io_load_coll(f,'Comp',into,No_,roc,Ao, comp_C)
          if(mpi_rnk==0)  print*, "...loaded composition"
      else
          call io_load_coll(f, 'C', into,No_,roc,Ao, comp_C)
          if(mpi_rnk==0)  print*, "...loaded composition from codensity"
      end if

      deallocate(roc)
      deallocate(ric)
!WD 10/9/2013, for nonmagnetic runs (b_mag_tstep==false) setting all magnetic energies to zero. They anyway are not updated
    if (.not. b_mag_tstep) then
     _loop_lm_begin(mag_Tor)
     N_=vel_Tor%D%N
     do n=1,N_
        mag_Tor%Re(n,nh)   = 0d0
        mag_Tor%Im(n,nh)   = 0d0 
        mag_Pol%Re(n,nh)   = 0d0
        mag_Pol%Im(n,nh)   = 0d0
        mag_icTor%Re(n,nh) = 0d0
        mag_icTor%Im(n,nh) = 0d0
        mag_icPol%Re(n,nh) = 0d0
        mag_icPol%Im(n,nh) = 0d0
     end do
     _loop_lm_end
    end if

      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
   end subroutine io_load_state



!--------------------------------------------------------------------------
!>  Load all .txt files
!!
!! Depending on parameters chosen, will look for diffusivity.txt.in, 
!! planet_model.txt.in, and confirm that they are correct.
!--------------------------------------------------------------------------
   subroutine io_load_radial_txtfiles(mes_oc, mes_ic)

      use meshs, only : rdom
      use parameters, only : b_boussinesq, b_diff_profile, i_Nic, i_N
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(inout) :: mes_oc !< outer core domain
      type(rdom), intent(inout) :: mes_ic !< inner core domain

      integer :: n, sta_input_file
      double precision :: readf1,readf2
 
!  Reads in the reference state file for anelastic cases where b_boussinesq =
!  .false.
!
!  If b_boussinesq = .true. and b_diff_profile = .true. read in diffusivity and its r-derivative
!  from diffusivity.txt. If b_boussinesq = .true. and b_diff_profile = .false. set diffusivity
!  to unity and its derivative to zero. 
      if (b_boussinesq) then
        if (b_diff_profile) then
          n=0
          open(unit=12, status='old', file=io_diffusfile,ACTION='read',iostat=sta_input_file)
          if (sta_input_file /= 0)    stop 'No diffusivity.txt.in file found in  &
                      the run directory'
          do while(.true.)
            read(12,*,end=9) readf1, readf2
            n=n+1
          end do
    9     close(12)
          if (mpi_rnk==0) print*,'number of rows read from diffusivity.txt.in', n
          if (i_Nic==0) then
            if (n/=i_N) then
                  stop 'incorrect diffusivity.txt file found in the run directory'
            end if
          else
            if (n/=(i_N+i_Nic)) then
                 stop 'incorrect diffusivity.txt file found in the run directory'
            end if
          end if 
          open(unit=12,file=io_diffusfile,iostat=sta_input_file,status='old')
          if (sta_input_file/=0) then
               stop 'No diffusivity.txt file found in the run directory'
          end if 
          do n=1,i_N
            read(12,*) mes_oc%eta(n,0), mes_oc%eta(n,1)
          end do
          if (i_Nic/=1) then
            do n=1,i_Nic
              read(12,*) mes_ic%eta(n,0), mes_ic%eta(n,1)
            end do
          end if 
          close(12)
          if (mpi_rnk==0) print*, ' loaded diffusivity file'
        else
          do n=1,i_N
            mes_oc%eta(n,0) = 1d0
            mes_oc%eta(n,1) = 0d0
          end do
          if (i_Nic/=1) then
            do n=1,i_Nic
              mes_ic%eta(n,0) = 1d0
              mes_ic%eta(n,1) = 0d0
            end do
          end if
        end if
      end if 
!  Read in the reference state data on every processor from planet_model.txt. This file must have i_N lines
!  if no inner core, i_N + i_Nic lines if there is an inner core.
!  In the Boussinesq case, no planet_model.txt file is needed.
!  If the run is anelastic,  the run will stop if no planet_model.txt is detected
!  in the run directory.
!  The first i_N lines of planet_model.txt must contain the seven quantities
!  rho, xi = (1/rho) d rho / dr, d xi / dr, Temperature, d Temp /  dr, eta and d
!  eta / dr at every radial mesh point in the outer core.
!  The next i_Nic lines (in the case of an inner core i_Nic > 1) must contain
!  eta and d eta / dr in the inner core, at every radial mesh point in the inner
!  core.
!  Matlab scripts are available for generating the reference state.
      if (.not.b_boussinesq) then

        n=0
        open(unit=11, status='old', file=io_refstatefile,ACTION='read',iostat=sta_input_file)
        if (sta_input_file /= 0)  stop 'No planet_model.txt.in file found in the run directory'
        do while(.true.)
          read(11,*,end=8) readf1 
          n=n+1
        end do
    8   close(11)
        if (mpi_rnk==0) print*,'number of rows read from planet_mag.txt.in', n
        if (i_Nic==1) then
          if (n/=i_N) then
             stop 'incorrect number of rows in planet_model.txt.in in the run directory'
          end if
        else
          if (n/=(i_N + i_Nic)) then
             stop 'incorrect number of rows in planet_model.txt.in in the run directory'
          end if
        end if 
        open(unit=11,file=io_refstatefile,iostat=sta_input_file,status='old')
        if (sta_input_file==0) then
           if(mpi_rnk==0) print*,'Loaded anelastic reference state file ...'
        else
!          if (mpi_rnk==0) print*,'Non-Bousssinesq runs must have a &
!               valid planet_model.txt in the run directory'
          stop 'Non-Bousssinesq runs must have a valid planet_model.txt.in in the run directory'
        end if
        do n=1,i_N
           read(11,*) mes_oc%rho(n,0), mes_oc%rho(n,1), mes_oc%rho(n,2), &
           mes_oc%Temp0(n,0), mes_oc%Temp0(n,1), mes_oc%eta(n,0), mes_oc%eta(n,1)
        end do
        if (i_Nic/=1) then
           do n=1,i_Nic
             read(11,*) mes_ic%eta(n,0), mes_ic%eta(n,1)
           end do
        end if
        close(11)
      end if




   end subroutine io_load_radial_txtfiles


!--------------------------------------------------------------------------
!>  Load codensity boundary condition and source info
!--------------------------------------------------------------------------
   subroutine io_load_codBC(mes_oc, cod_C, cod_S, cod_bcRe, cod_bcIm)

      use netcdf
      use meshs, only : rdom
      use parameters, only : i_Kl
      use dyn_mpi
      use variables, only : var_H, var_coll_init

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(in) :: cod_C !< outer core domain
      double precision, intent(out) :: cod_S(:) !< codensity source term, S(r)
      double precision, intent(out) :: cod_bcRe(:,0:) !< codensity BCs, real part
      double precision, intent(out) :: cod_bcIm(:,0:) !< codensity BCs, imaginary part

      integer :: e, f, i, rd, i9,j9
      integer :: n, No, No_
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: Ao(lda,lda,mes_oc%N)
      double precision, allocatable :: roc(:), fn(:)
      logical :: into, intxt, intxt1
      
      ! first open the file, with file handle f
      ! e is the status returned by the open attempt.
      e=nf90_open(io_codBCfile, nf90_nowrite, f)  

      ! if e is not noerr, there was an error opening
      ! codBC.cdf.in, so stop program
      if(e/=nf90_noerr) then
         stop 'codBC file not found,but inhomogeneous boundaries specified!'
         return

      ! if e == noerr, then codBC.cdf.in was opened ok,
      ! so print confirmation and continue.
      else
         if(mpi_rnk==0) print *, 'found codBC.cdf.in, loading...'
      end if
      ! codBC.cdf.in now loaded, with id 'f'
      
      ! No now equals number of points in outer core mesh
      No = mes_oc%N

      ! for file f (codBC.cdf.in), find the id (rd) of the dimension
      ! 'r' (radius)
      e=nf90_inq_dimid(f,    'r', rd)
      ! for file f (codBC.cdf.in), find length of dimemnsion rd
      ! ('r'=radius), and put in No_
      e=nf90_inquire_dimension(f,   rd, len=No_)

      ! allocate fn with length No_
      allocate(fn (No_))
      ! allocate roc with length No_
      allocate(roc(No_))

      ! put id of varialbe 'r' from file 'f' into 'i'
      e=nf90_inq_varid(f, 'r', i)
      ! put values of 'r' from 'f' into roc
      e=nf90_get_var(f, i, roc)
      ! does radial grid size specified in parameters match radial
      ! grid size in codBC.cdf.in ?
      into = (No_/=No)

      ! if No_ and No are different, interpolate the weights of the 
      ! finite difference discretisation
      if(.not.into) into = (maxval(dabs(mes_oc%r(1:No,1)-roc))>1d-8)
      if(into) then
         if(mpi_rnk==0) print*,' N  :', No_, ' --> ',No 
         call io_interp_wts(No_,roc,No,mes_oc%r(1,1), Ao)
      end if

      ! initialise fn(No_) to 0
      fn = 0d0

      ! load 'S' from codBC.cdf.in and put it in
      ! cod_S, interpolating if necessary
      e=nf90_inq_varid(f, 'S', i)
      e=nf90_get_var(f, i, fn)
      if(into) then
         call io_interp(No_,roc,fn,Ao,No,mes_oc%r(1,1), cod_S)
      else
         cod_S = fn
      end if

      ! Load 'Cbc' from 'codBC.cdf.in', and put it in cq. 
      do_random_init = .false.
      ! first initialise cq as mes_oc object.
      call var_coll_init(mes_oc, var_H, cq)

      ! now load 'Cbc' from file id f into cq%Re and cq%Im
      call io_load_coll(f, 'Cbc', .false., 2, roc, Ao, cq)
      cod_bcRe = cq%Re(:2,:)
      cod_bcIm = cq%Im(:2,:)

      deallocate(fn )
      deallocate(roc)

      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
      call io_save_codBCtxt(cod_C, cod_bcRe, cod_bcIm)
      call io_save_codStxt(mes_oc, cod_S)

   end subroutine io_load_codBC


!--------------------------------------------------------------------------
!>  Load codBC.txt.in (if it exists) 
!--------------------------------------------------------------------------
   subroutine io_load_codBCtxt(cod_C, cod_bcRe, cod_bcIm)
      use dyn_mpi
      type(coll), intent(in) :: cod_C !< codenesity
      double precision, intent(out) :: cod_bcRe(:,0:) !< codensity BCs, real part
      double precision, intent(out) :: cod_bcIm(:,0:) !< codensity BCs, imaginary part

      double precision :: dRe, dIm1
      integer :: l__h,m__h,riro
      logical :: intxt
      _loop_lm_vars

      inquire(file=io_codBCtxt, exist=intxt)
      if ( .not. intxt ) then
         stop ' codBC.txt.in not found '
         return
      end if
      cod_bcRe=0d0
      cod_bcIm=0d0
      open(99, status='old', file=io_codBCtxt)
      do while(.true.)
        read(99,*,end=9) riro, l__h, m__h, dRe, dIm1
        if (riro<1.or.riro>2) stop 'codBCtxt: riro error'
        _loop_lm_begin(cod_C)
          if (l__h == l.and. m__h==m_) then 
            cod_bcRe(riro,nh) = dRe
            cod_bcIm(riro,nh) = dIm1
          end if
        _loop_lm_end
      end do
    9 close(99)

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_codBCtxt


!--------------------------------------------------------------------------
!>  Save codBC.txt.in (if it exists) 
!-------------------------------------------------------------------------
   subroutine io_save_codBCtxt(cod_C, cod_bcRe, cod_bcIm)
      type(coll), intent(in) :: cod_C !< codensity
      double precision, intent(in) :: cod_bcRe(:,0:) !< codensity BCs, real part
      double precision, intent(in) :: cod_bcIm(:,0:) !< codnsity BCs, imaginary part

      double precision :: dRe, dIm
      _loop_lm_vars

      if(_Np/=1) return  !implemented for _Np=1 only
      open(99,status='unknown',file='codBC.txt.in')
      _loop_lm_begin(cod_C)
         if(cod_bcRe(1,nh)/=0d0 .or. cod_bcIm(1,nh)/=0d0)  &
            write(99,'(I2,2I6,2e28.20)')  &
               1, l, m_, cod_bcRe(1,nh), cod_bcIm(1,nh)
         if(cod_bcRe(2,nh)/=0d0 .or. cod_bcIm(2,nh)/=0d0)  &
            write(99,'(I2,2I6,2e28.20)')  &
               2, l, m_, cod_bcRe(2,nh), cod_bcIm(2,nh)
      _loop_lm_end
      close(99)
           
   end subroutine io_save_codBCtxt

!--------------------------------------------------------------------------
!>  Load internal heating from codS.txt.in
!--------------------------------------------------------------------------
   subroutine io_load_codStxt(mes_oc, cod_S)
      use meshs, only : rdom
      use parameters, only : i_KL
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      double precision, intent(out) :: cod_S(:) !< codensity source term, S(r)

      integer :: n, No, No_
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: Ao(lda,lda,mes_oc%N)
      double precision, allocatable :: roc(:), fn(:)
      logical :: into, intxt

      inquire(file=io_codStxt, exist=intxt)
      if ( .not. intxt) then
          stop ' codS.txt.in not found '
          return
      end if

      open(99,status='old',file='codS.txt.in')
      read(99,*) n
      if(n==0) return 		! no internal heating
      if(n==1) then		! uniform heating
         read(99,*) d
         cod_S = d
         return
      end if
      				! S=S(r)
      No_= n
      No = mes_oc%N
      allocate(fn (No_))
      allocate(roc(No_))
      do n = 1, No_
         read(99,*) roc(n), fn(n)
      end do

      into = (No_/=No)
      if(.not.into) into = (maxval(dabs(mes_oc%r(1:No,1)-roc))>1d-8)
      if(into) then
         if(mpi_rnk==0) print*,' N  :', No_, ' --> ',No 
         call io_interp_wts(No_,roc,No,mes_oc%r(1,1), Ao)
         call io_interp(No_,roc,fn,Ao,No,mes_oc%r(1,1), cod_S)
      else
         cod_S = fn
      end if
      
      if(.not.into) return
      if(mpi_rnk==0) call system('cp codS.txt.in codS.txt.old')
      call io_save_codStxt(mes_oc, cod_S)
   
   end subroutine io_load_codStxt 


!--------------------------------------------------------------------------
!>  Save internal heating to codS.txt.in
!-------------------------------------------------------------------------
   subroutine io_save_codStxt(mes_oc, cod_S)
      use meshs, only : rdom
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      double precision, intent(in) :: cod_S(:) !< codensity source term, S(r)

      double precision :: Smn, Smx
      integer :: n
      
      if(mpi_rnk/=0) return
      open(99,status='unknown',file='codS.txt.in')
      Smn = minval(cod_S)
      Smx = maxval(cod_S)
      if(Smn==0d0 .and. Smx==0d0) then
         write(99,*) 0
      else if(Smn==Smx) then
         write(99,*) 1
         write(99,'(e28.20)') cod_S(1)
      else
         write(99,*) mes_oc%N
         do n = 1, mes_oc%N
            write(99,'(2e28.20)') mes_oc%r(n,1), cod_S(n)
         end do
      end if
      close(99)
   
   end subroutine io_save_codStxt


!--------------------------------------------------------------------------
!>  Load composition boundary condition and source info
!--------------------------------------------------------------------------
   subroutine io_load_compBC(mes_oc, comp_C, comp_S, comp_bcRe, comp_bcIm)
      use netcdf
      use meshs, only : rdom
      use parameters, only : i_KL
      use dyn_mpi
      use variables, only : var_H, var_coll_init

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(in) :: comp_C !< composition
      double precision, intent(out) :: comp_S(:) !< composition source term, S(r)
      double precision, intent(out) :: comp_bcRe(:,0:) !< compostion BCs, real part
      double precision, intent(out) :: comp_bcIm(:,0:) !< composition BCs, imaginary part

      integer :: e, f, i, rd, i9,j9
      integer :: n, No, No_
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: Ao(lda,lda,mes_oc%N)
      double precision, allocatable :: roc(:), fn(:)
      logical :: into, intxt, intxt1

      ! first open the file, with file handle f
      ! e is the status returned by the open attempt.
      e=nf90_open(io_compBCfile, nf90_nowrite, f)

      ! if e is not noerr, there was an error opening
      ! compBC.cdf.in, so stop program
      if(e/=nf90_noerr) then
         stop 'compBC file not found,but inhomogeneous boundaries specified!'
         return

      ! if e == noerr, then compBC.cdf.in was opened ok,
      ! so print confirmation and continue.
      else
         if(mpi_rnk==0) print *, 'found compBC.cdf.in, loading...'
      end if
      ! compBC.cdf.in now loaded, with id 'f'

      ! No now equals number of points in outer core mesh
      No = mes_oc%N

      ! for file f (codBC.cdf.in), find the id (rd) of the dimension
      ! 'r' (radius)
      e=nf90_inq_dimid(f,    'r', rd)
      ! for file f (compBC.cdf.in), find length of dimemnsion rd
      ! ('r'=radius), and put in No_
      e=nf90_inquire_dimension(f,   rd, len=No_)

      ! allocate fn with length No_
      allocate(fn (No_))
      ! allocate roc with length No_
      allocate(roc(No_))

      ! put id of varialbe 'r' from file 'f' into 'i'
      e=nf90_inq_varid(f, 'r', i)
      ! put values of 'r' from 'f' into roc
      e=nf90_get_var(f, i, roc)
      ! does radial grid size specified in parameters match radial
      ! grid size in compBC.cdf.in ?
      into = (No_/=No)

      ! if No_ and No are different, interpolate the weights of the
      ! finite difference discretisation
      if(.not.into) into = (maxval(dabs(mes_oc%r(1:No,1)-roc))>1d-8)
      if(into) then
         if(mpi_rnk==0) print*,' N  :', No_, ' --> ',No
         call io_interp_wts(No_,roc,No,mes_oc%r(1,1), Ao)
      end if

      ! initialise fn(No_) to 0
      fn = 0d0

      ! load 'S' from compBC.cdf.in and put it in
      ! comp_S, interpolating if necessary
      e=nf90_inq_varid(f, 'S', i)
      e=nf90_get_var(f, i, fn)
      if(into) then
         call io_interp(No_,roc,fn,Ao,No,mes_oc%r(1,1), comp_S)
      else
         comp_S = fn
      end if

      ! Load 'Cbc' from 'compBC.cdf.in', and put it in cq.
      do_random_init = .false.
      ! first initialise cq as mes_oc object.
      call var_coll_init(mes_oc, var_H, cq)

      ! now load 'Compbc' from file id f into cq%Re and cq%Im
      call io_load_coll(f, 'Compbc', .false., 2, roc, Ao, cq)
      comp_bcRe = cq%Re(:2,:)
      comp_bcIm = cq%Im(:2,:)

      deallocate(fn )
      deallocate(roc)

      e=nf90_close(f)
#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
      call io_save_compBCtxt(comp_C, comp_bcRe, comp_bcIm)
      call io_save_compStxt(mes_oc, comp_S)

   end subroutine io_load_compBC


!--------------------------------------------------------------------------
!>  Load compBC.txt.in (if it exists)
!--------------------------------------------------------------------------
   subroutine io_load_compBCtxt(comp_C, comp_bcRe, comp_bcIm)
      use dyn_mpi

      type(coll), intent(in) :: comp_C  !< composition
      double precision, intent(out) :: comp_bcRe(:,0:) !< composition BCs, real part
      double precision, intent(out) :: comp_bcIm(:,0:) !< composition BCs, imaginary part

      double precision :: dRe, dIm1
      integer :: l__h,m__h,riro
      logical :: intxt
      _loop_lm_vars

      inquire(file=io_compBCtxt, exist=intxt)
      if ( .not. intxt ) then
         stop ' compBC.txt.in not found '
         return
      end if
      comp_bcRe=0d0
      comp_bcIm=0d0
      open(99, status='old', file=io_compBCtxt)
      do while(.true.)
        read(99,*,end=9) riro, l__h, m__h, dRe, dIm1
        if (riro<1.or.riro>2) stop 'compBCtxt: riro error'
        _loop_lm_begin(comp_C)
          if (l__h == l.and. m__h==m_) then
            comp_bcRe(riro,nh) = dRe
            comp_bcIm(riro,nh) = dIm1
          end if
        _loop_lm_end
      end do
    9 close(99)

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif

   end subroutine io_load_compBCtxt


!--------------------------------------------------------------------------
!>  Save compBC.txt.in (if it exists)
!-------------------------------------------------------------------------
   subroutine io_save_compBCtxt(comp_C, comp_bcRe, comp_bcIm)
      type(coll), intent(in) :: comp_C !< composition
      double precision, intent(in) :: comp_bcRe(:,0:) !< composition BCs, real part
      double precision, intent(in) :: comp_bcIm(:,0:) !< composition BCs, imaginary part

      double precision :: dRe, dIm
      _loop_lm_vars

      if(_Np/=1) return  !implemented for _Np=1 only
      open(99,status='unknown',file='compBC.txt.in')
      _loop_lm_begin(comp_C)
         if(comp_bcRe(1,nh)/=0d0 .or. comp_bcIm(1,nh)/=0d0)  &
            write(99,'(I2,2I6,2e28.20)')  &
               1, l, m_, comp_bcRe(1,nh), comp_bcIm(1,nh)
         if(comp_bcRe(2,nh)/=0d0 .or. comp_bcIm(2,nh)/=0d0)  &
            write(99,'(I2,2I6,2e28.20)')  &
               2, l, m_, comp_bcRe(2,nh), comp_bcIm(2,nh)
      _loop_lm_end
      close(99)

   end subroutine io_save_compBCtxt

!-------------------------------------------------------------------------
!>   Flayer compostion boundary conditions and source
!-------------------------------------------------------------------------
   subroutine io_comp_flayer(mes_oc, comp_C, comp_S, comp_bcRe, comp_bcim)
      use meshs, only : rdom
      use parameters, only : d_flayer_radius, d_transition_width, d_fl_q,&
                             d_fl_kl, d_fl_ku, d_fl_cps, i_N
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(in) :: comp_C !< composition
      double precision, intent(out) :: comp_S(:) !< composition source term, S(r)
      double precision, intent(out) :: comp_bcRe(:,0:) !< compostion BCs, real part
      double precision, intent(out) :: comp_bcim(:,0:) !< compostion BCs, imaginary part

      double precision :: Sl, Su, width, q, kl, ku, Cps, rs, ro, ri, r
      integer :: N, l__h, m__h
      _loop_lm_vars
      N = 0
      ri = mes_oc%r(1,1)
      ro = mes_oc%r(i_N,1)
      rs = d_flayer_radius
      width = d_transition_width
      q = d_fl_q
      kl = d_fl_kl
      ku = d_fl_ku
      Cps = d_fl_Cps
      Sl = -3.0d0*(q*ri**2 - kl*Cps*rs**2)/(rs**3 - ri**3) / kl
      Su = -3.0d0*(          ku*Cps*rs**2)/(ro**3 - rs**3) / ku
!      do n = 1, i_N
!         r=mes_oc%r(n,1)
!         comp_S(n) = ((Sl - Su)*(1.0d0 - tanh((r-rs)/width)) / 2.0d0) + Su 
!      enddo 
!      do n = 1, i_N
!         r = mes_oc%r(n,1)
!         if (r < rs) then
!            comp_S(n) = Sl
!         else
!            comp_S(n) = Su
!       endif
!      enddo

      if (width > 0.0d0) then
      ! Smooth transition with tanh
        do n = 1, i_N
           r = mes_oc%r(n,1)
           comp_S(n) = ((Sl - Su)*(1.0d0 - tanh((r - rs)/width)) / 2.0d0) + Su
        end do
      else
      ! Sharp Heaviside step
        do n = 1, i_N
           r = mes_oc%r(n,1)
           if (r <= rs) then
              comp_S(n) = Sl
           else
              comp_S(n) = Su
           end if
        end do
      end if
      if (mpi_rnk==0) print *, "created flayer comp_S"
      _loop_lm_begin(comp_C)
        if (l == 0 .and. m_==0) then
          comp_bcRe(1,nh) = -q/kl 
          comp_bcRe(2,nh) = 0.0d0
        end if
      _loop_lm_end
      if (mpi_rnk==0) print *, "set compBC to flayer values: ", comp_bcRe(1,0), comp_bcRe(2,0)
      call io_save_compStxt(mes_oc, comp_S)
      call io_save_compBCtxt(comp_c, comp_bcre, comp_bcIm)
   end subroutine io_comp_flayer
      
!--------------------------------------------------------------------------
!>  Load internal composition source from compS.txt.in
!--------------------------------------------------------------------------
   subroutine io_load_compStxt(mes_oc, comp_S)
      use meshs, only : rdom
      use parameters, only : i_KL
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      double precision, intent(out) :: comp_S(:) !< composition source term, S(r)

      integer :: n, No, No_
      double precision :: d
      integer, parameter :: lda = 2*i_KL+1
      double precision :: Ao(lda,lda,mes_oc%N)
      double precision, allocatable :: roc(:), fn(:)
      logical :: into, intxt

      inquire(file=io_compStxt, exist=intxt)
      if ( .not. intxt ) then
         stop ' compS.txt.in not found '
         return
      end if
      open(99,status='old',file='compS.txt.in')
      read(99,*) n
      if(n==0) return           ! no internal heating
      if(n==1) then             ! uniform heating
         read(99,*) d
         comp_S = d
         return
      end if
      ! S=S(r)
      No_= n
      No = mes_oc%N
      allocate(fn (No_))
      allocate(roc(No_))
      do n = 1, No_
         read(99,*) roc(n), fn(n)
      end do

      into = (No_/=No)
      if(.not.into) into = (maxval(dabs(mes_oc%r(1:No,1)-roc))>1d-8)
      if(into) then
         if(mpi_rnk==0) print*,' N  :', No_, ' --> ',No
         call io_interp_wts(No_,roc,No,mes_oc%r(1,1), Ao)
         call io_interp(No_,roc,fn,Ao,No,mes_oc%r(1,1), comp_S)
      else
         comp_S = fn
      end if

      if(.not.into) return
      if(mpi_rnk==0) call system('cp compS.txt.in compS.txt.old')
      call io_save_compStxt(mes_oc, comp_S)

   end subroutine io_load_compStxt


!--------------------------------------------------------------------------
!>  Save internal heating to compS.txt.in
!-------------------------------------------------------------------------
   subroutine io_save_compStxt(mes_oc, comp_S)
      use meshs, only : rdom
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      double precision, intent(in) :: comp_S(:) !< composition source term, S(r)

      double precision :: Smn, Smx
      integer :: n

      if(mpi_rnk/=0) return
      open(99,status='unknown',file='compS.txt.in')
      Smn = minval(comp_S)
      Smx = maxval(comp_S)
      if(Smn==0d0 .and. Smx==0d0) then
         write(99,*) 0
      else if(Smn==Smx) then
         write(99,*) 1
         write(99,'(e28.20)') comp_S(1)
      else
         write(99,*) mes_oc%N
         do n = 1, mes_oc%N
            write(99,'(2e28.20)') mes_oc%r(n,1), comp_S(n)
         end do
      end if
      close(99)

   end subroutine io_save_compStxt


!--------------------------------------------------------------------------
!>  Load poloidal magnetic field boundary condition 
!--------------------------------------------------------------------------
   subroutine io_load_magPolBC(mag_pol, mag_bc_PolRe, mag_bc_PolIm)
      use dyn_mpi
      type(coll), intent(in) :: mag_pol !< poloidal component of magnetic field
      double precision, intent(out) :: mag_bc_PolRe(:,:) !< poloidal magnetic BCs, real part
      double precision, intent(out) :: mag_bc_PolIm(:,:) !< poloidal magnetic BCs, imaginary part

      double precision :: dRe, dIm1
      integer :: l__h,m__h,riro
      _loop_lm_vars

      mag_bc_PolRe=0d0
      mag_bc_PolIm=0d0
      open(99, status='old', file=io_magPolBCfile)
      do while(.true.)
        read(99,*,end=9) riro, l__h, m__h, dRe, dIm1
        if (riro<1.or.riro>2) stop 'codBCtxt: riro error'
        _loop_lm_begin(mag_Pol)
          if (l__h == l.and. m__h==m_) then 
            mag_bc_PolRe(riro,nh) = dRe
            mag_bc_PolIm(riro,nh) = dIm1
          end if
        _loop_lm_end
      end do
    9 close(99)

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
#endif
  

   end subroutine io_load_magPolBC


!--------------------------------------------------------------------------
!>  index of l,m mode in state file
!--------------------------------------------------------------------------
   subroutine io_state_nh(l,m_,iL,iM,iMp, nh)
      integer, intent(in)  :: l !< spherical harmonic degree
      integer, intent(in)  :: m_ !< actual spherical harmonic order
      integer, intent(in)  :: iL !< max spherical harmonic degree
      integer, intent(in)  :: iM !< max spherical harmonic order
      integer, intent(in)  :: iMp !< azimuthal periodicity
      integer, intent(out) :: nh !< harmonic index
      integer :: m
      if(modulo(m_,iMp)/=0 .or. l>=iL .or. m_>=iM*iMp) then
         nh = -1
         return
      end if
      m  = m_/iMp
      nh = iL*m - iMp*(m-1)*m/2 + l-m_!+1 - 1
   end subroutine io_state_nh


!--------------------------------------------------------------------------
!>  Load coll variable
!--------------------------------------------------------------------------
  subroutine io_load_coll(f,nm,interp,N__,r,W, a)
      use netcdf
      use parameters, only : i_KL, i_pH1, i_L, i_M, i_Mp, d_init_pert
      use dyn_mpi, only : mpi_rnk

      type (coll),      intent(inout) :: a !< scalar field in spectral space
      integer,          intent(in) :: f !< netcdf file id
      character(*),     intent(in) :: nm !< netcdf variable name
      logical,          intent(in) :: interp !< interpolate to new grid size?
      integer,          intent(in) :: N__ !< radial grid size in netcdf file
      double precision, intent(in) :: r(N__) !< radial grid in netcdf file
      double precision, intent(in) :: W(2*i_KL+1,2*i_KL+1,*) !< interpolation weights

      integer :: i !< netcdf variable id

      double precision :: dRe(N__,i_pH1+1), dIm(N__,i_pH1+1)
      integer :: nh_,nh__,n
      integer :: L__, M__, Mp__      
      integer :: e,j
      _loop_lm_vars
         
      e=nf90_inq_varid(f,nm, i)
      e=nf90_get_att(f,i, 'L',  L__)
      e=nf90_get_att(f,i, 'M',  M__)
      e=nf90_get_att(f,i,'Mp', Mp__)
      if(mpi_rnk==0) then
         if(L__/=i_L) print*, nm, ' L :', L__,' --> ',i_L 
         if(M__/=i_M) print*, nm, ' M :', M__,' --> ',i_M 
         if(Mp__/=i_Mp .and. M__/=1)  &
                      print*, nm, ' Mp:',Mp__,' --> ',i_Mp     
      end if

      a%Re = 0d0
      a%Im = 0d0
      n = 1

      _loop_lm_begin(a)      
         call io_state_nh(l,m_,L__,M__,Mp__, nh__)
         if(nh__==-1) then
            ! mode is not in state file
            if(do_random_init .and. d_init_pert>0d0) then
               ! add noise to this mode
               call random_number(a%Re(:,nh))
               call random_number(a%Im(:,nh))
               a%Re(:,nh)=2d0*d_init_pert*(a%Re(:,nh)-0.5d0)
               a%Im(:,nh)=2d0*d_init_pert*(a%Im(:,nh)-0.5d0)
            end if
            n = 1
         else
            ! mode is in state file
            if(nh>=a%H%pH1) nh_ = -1
            if(nh <a%H%pH1) call io_state_nh(  &
               a%H%l(a%H%pH0+nh+1),a%H%m(a%H%pH0+nh+1),L__,M__,Mp__, nh_)
            if(nh__+1==nh_) then
               ! mode is contiguous with next mode in state file
               n = n + 1
            else
               ! last n modes were contiguous, load them!
               if(.not.interp) then
                  j = nh__+1 - (n-1)
                  e=nf90_get_var(f,i, a%Re(1:N__,nh-n+1:nh),start=(/1,j,1/))
                  e=nf90_get_var(f,i, a%Im(1:N__,nh-n+1:nh),start=(/1,j,2/))
               else
                  j = nh__+1 - (n-1)
                  e=nf90_get_var(f,i, dRe(1:N__,1:n),start=(/1,j,1/))
                  e=nf90_get_var(f,i, dIm(1:N__,1:n),start=(/1,j,2/))
                  do j = 1, n
                     call io_interp(N__,r,dRe(1,j),W,  &
                        a%D%N,a%D%r(1,1), a%Re(1,nh-n+j))
                     call io_interp(N__,r,dIm(1,j),W,  &
                        a%D%N,a%D%r(1,1), a%Im(1,nh-n+j))
                  end do
               end if
               n = 1
            end if
         end if
      _loop_lm_end      
   end subroutine io_load_coll

!--------------------------------------------------------------------------
!>  interpolate; map both onto xrange [0,1] and return weights
!--------------------------------------------------------------------------
   subroutine io_interp_wts(ni,xi,no,xo, A)
      use parameters, only : i_KL
      use meshs, only : mes_mat_invert
      integer, parameter :: lda = 2*i_KL+1

      integer,          intent(in)  :: ni !< number of grid points of incoming data
      integer,          intent(in)  :: no !< number of grid points to interpolate to
      double precision, intent(in)  :: xi(ni) !< grid of incoming data
      double precision, intent(in)  :: xo(no) !< grid to interpolate to
      double precision, intent(out) :: A(lda,lda,no) !< weights for interpolation

      double precision :: xi1(ni), xo1(no)
      integer :: n,nn,i,j,l,r
         			! map both onto [0,1] for convenience
      xi1 = (xi-xi(1))/(xi(ni)-xi(1))
      xo1 = (xo-xo(1))/(xo(no)-xo(1))
      
      do n = 1, no
         j = 1
         do while(xi1(j) < xo1(n)-1d-8)
           j = j+1
         end do
         l = max(1,j-i_KL)
         r = min(j+i_KL,ni)
         nn = r-l+1
         do i = 1, nn
            A(i,1,n) = 1d0
         end do
         do j = 2, nn
            do i = 1, nn
               A(i,j,n) = A(i,j-1,n) * (xi1(l+i-1)-xo1(n)) / dble(j-1) 
            end do
         end do
         call mes_mat_invert(nn,A(1,1,n),lda)
      end do

   end subroutine io_interp_wts


!--------------------------------------------------------------------------
!>  interpolate, given weights from io_interp_wts()
!--------------------------------------------------------------------------
   subroutine io_interp(ni,xi,fi,A,no,xo, fo)
      use parameters, only : i_KL
      integer, parameter :: lda = 2*i_KL+1

      integer,          intent(in)  :: ni !< number of grid points, incoming data
      integer,          intent(in)  :: no !< number of grid points, output
      double precision, intent(in)  :: xi(ni) !< grid of incoming data
      double precision, intent(in)  :: fi(ni) !< input data
      double precision, intent(in)  :: xo(no) !< grid of output
      double precision, intent(in)  :: A(lda,lda,no) !< weights for interpolation
      double precision, intent(out) :: fo(no) !< interpolated data on output grid

      double precision :: xi1(ni), xo1(no)
      integer :: n,nn,i,j,l,r
         			! map both onto [0,1] for convenience
      xi1 = (xi-xi(1))/(xi(ni)-xi(1))
      xo1 = (xo-xo(1))/(xo(no)-xo(1))
      
      do n = 1, no
         j = 1
         do while(xi1(j) < xo1(n)-1d-8)
           j = j+1
         end do
         l = max(1,j-i_KL)
         r = min(j+i_KL,ni)
         nn = r-l+1
         fo(n) = dot_product(A(1,1:nn,n),fi(l:r))
      end do

   end subroutine io_interp


!-------------------------------------------------------------------------
!> extend into inner core:  f = a r^3 + b r^2  
!! f=f'=0 at r=0, f,f' cts at icb
!-------------------------------------------------------------------------
   subroutine io_ic_extend(mes_oc,mes_ic,nm,oc,ic)
      use meshs, only : rdom
      use dyn_mpi, only : mpi_rnk
      use parameters, only : i_KL
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      character(*), intent(in)  :: nm !< name of field being extended
      type (coll),  intent(in)  :: oc !< outer core field
      type (coll),  intent(out) :: ic !< inner core field
      double precision :: ri,ri2,ri3, a,b,f,f1
      integer :: nhm,nh,j

      if(mpi_rnk==0)  print*, nm, ' extending B onto inner core'
      ri  = mes_oc%r(1,1)
      ri2 = mes_oc%r(1,2)
      ri3 = mes_oc%r(1,3)
      nhm = oc%H%pH1
      
      do nh = 0, nhm
         f = oc%Re(1,nh)
         f1 = 0d0
         do j = 1, 1+i_KL
            f1 = f1 + mes_oc%dr(1)%M(i_KL+1+1-j,j)*oc%Re(j,nh)
         end do
         a = ( f1*ri - 2d0*f ) / ri3
         b = ( f - a*ri3) / ri2
         ic%Re(:,nh) = a*mes_ic%r(:,3) + b*mes_ic%r(:,2)
      end do
      do nh = 0, nhm
         f = oc%Im(1,nh)
         f1 = 0d0
         do j = 1, 1+i_KL
            f1 = f1 + mes_oc%dr(1)%M(i_KL+1+1-j,j)*oc%Im(j,nh)
         end do
         a = ( f1*ri - 2d0*f ) / ri3
         b = ( f - a*ri3) / ri2
         ic%Im(:,nh) = a*mes_ic%r(:,3) + b*mes_ic%r(:,2)
      end do

   end subroutine io_ic_extend

!--------------------------------------------------------------------------
!>  Open a netcdf file with parallel access
!--------------------------------------------------------------------------
#ifdef _PNCDF
   subroutine io_create_par(fname, f)
      use netcdf
      use dyn_mpi
      character(*), intent(in)  :: fname !< filename
      integer,      intent(out) :: f !< netcdf file id
      integer :: e,info,cmode
      info =mpi_info_null
      cmode=ior(nf90_mpiio,nf90_clobber)
      cmode=ior(cmode,nf90_netcdf4)
      e=nf90_create_par(fname, cmode, mpi_comm_world, info, f)
      if(e/=nf90_noerr)  print*, nf90_strerror(e)
      if(e/=nf90_noerr)  stop
   end subroutine io_create_par
#endif


!--------------------------------------------------------------------------
!>  Save state
!--------------------------------------------------------------------------
   subroutine io_save_state(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol, mag_icTor, mag_icPol,&
                       cod_C, comp_C,rot_omega)
      use netcdf
      use meshs, only : rdom
      use parameters, only : i_L, i_M, i_Mp, d_E, d_Pm, d_Ra, d_Ra_comp, d_Pr, d_rratio, i_N, i_Nic, d_Sc
      use dyn_mpi, only : mpi_rnk
      use timestep, only : tim_t, tim_dt

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      type(coll), intent(in) :: vel_tor !< toroidal component of velocity
      type(coll), intent(in) :: vel_pol !< poloidal component of velocity
      type(coll), intent(in) :: mag_tor !< toroidal component of magnetic field
      type(coll), intent(in) :: mag_pol !< poloidal component of magnetic field
      type(coll), intent(in) :: mag_icTor !< toroidal component of mag field, (inner core)
      type(coll), intent(in) :: mag_icPol !< poloidal component of mag field, (inner core)
      type(coll), intent(in) :: cod_C !< codensity
      type(coll), intent(in) :: comp_C !< composition
      double precision, intent(in) :: rot_omega !< inner core rotation rate

      character(4) :: cnum
      character*35 :: fname
      integer :: e, f
      integer :: H
      integer :: rd,icrd, uTHd,uPHd,BTHd,BPHd, ReImd
      integer :: H_dim
      integer :: dt, r,icr, omega, uT,uP,C,Comp,BT,BP,icBT,icBP

      H = i_L*i_M - i_Mp*(i_M-1)*i_M/2 !- 1
            
      write(cnum,'(I4)') io_save1
      if(io_save1<1000) cnum(1:1) = '0'
      if(io_save1< 100) cnum(2:2) = '0'
      if(io_save1<  10) cnum(3:3) = '0'
      fname = 'state'//cnum//'.cdf.dat'

      if(mpi_rnk==0) then
         print*, ' saving state'//cnum//'  t=', tim_t

         write(66,*)  ' saving state'//cnum//'  t=', tim_t
         CALL flush(66)

#ifdef _PNCDF
      end if
      call io_create_par(fname, f)
#else
         e=nf90_create('state'//cnum//'.cdf.dat', nf90_clobber, f)
#endif
         e=nf90_put_att(f, nf90_global,  't', tim_t)
         e=nf90_put_att(f, nf90_global, 'Ro', d_E/d_Pm)
         e=nf90_put_att(f, nf90_global,  'E', d_E)
         e=nf90_put_att(f, nf90_global, 'Ra', d_Ra)
         e=nf90_put_att(f, nf90_global, 'Ra_comp', d_Ra_comp)
         e=nf90_put_att(f, nf90_global,  'Pr', d_Pr)
         e=nf90_put_att(f, nf90_global, 'Sc', d_Sc)
         e=nf90_put_att(f, nf90_global,  'Pm', d_Pm)
         e=nf90_put_att(f, nf90_global,  'q',  d_Pm/d_Pr)
         e=nf90_put_att(f, nf90_global, 'riro', d_rratio)
         e=nf90_put_att(f, nf90_global, 'LSD_version', code_version)

!          e=nf90_def_dim(f,   'r', mes_oc%N, rd)
!          e=nf90_def_dim(f, 'icr', mes_ic%N, icrd)
!          e=nf90_def_dim(f, 'uTH', H, uTHd)
!          e=nf90_def_dim(f, 'uPH', H, uPHd)
!          e=nf90_def_dim(f, 'BTH', H, BTHd)
!          e=nf90_def_dim(f, 'BPH', H, BPHd)
!          e=nf90_def_dim(f, 'ReIm',  2, ReImd)

         e=nf90_def_dim(f,   'r', mes_oc%N, rd)
         e=nf90_def_dim(f, 'icr', mes_ic%N, icrd)
         e=nf90_def_dim(f, 'H', H, H_dim)
         e=nf90_def_dim(f, 'ReIm',  2, ReImd)

         e=nf90_def_var(f,   'r', nf90_double, (/rd/), r)
         e=nf90_def_var(f, 'icr', nf90_double, (/icrd/), icr)
         e=nf90_def_var(f,       'dt', nf90_double, dt)
         e=nf90_def_var(f,  'icOmega', nf90_double, omega)

!          call io_define_coll(f,  'uT',(/rd,uTHd,ReImd/),vel_Tor, uT)
!          call io_define_coll(f,  'uP',(/rd,uPHd,ReImd/),vel_Pol, uP)
!          call io_define_coll(f,   'C',(/rd,uPHd,ReImd/),cod_C,   C)
!          call io_define_coll(f,  'BT',(/rd,BTHd,ReImd/),mag_Tor, BT)
!          call io_define_coll(f,  'BP',(/rd,BPHd,ReImd/),mag_Pol, BP)
!          call io_define_coll(f,'icBT',(/rd,BTHd,ReImd/),mag_icTor, icBT)
!          call io_define_coll(f,'icBP',(/rd,BPHd,ReImd/),mag_icPol, icBP)

         call io_define_coll(f,  'uT',(/rd,H_dim,ReImd/),vel_Tor, uT)
         call io_define_coll(f,  'uP',(/rd,H_dim,ReImd/),vel_Pol, uP)
         call io_define_coll(f,   'C',(/rd,H_dim,ReImd/),cod_C,   C)
         call io_define_coll(f,   'Comp',(/rd,H_dim,ReImd/),comp_C,   Comp)
         call io_define_coll(f,  'BT',(/rd,H_dim,ReImd/),mag_Tor, BT)
         call io_define_coll(f,  'BP',(/rd,H_dim,ReImd/),mag_Pol, BP)
         call io_define_coll(f,'icBT',(/rd,H_dim,ReImd/),mag_icTor, icBT)
         call io_define_coll(f,'icBP',(/rd,H_dim,ReImd/),mag_icPol, icBP)

         e=nf90_enddef(f)
#ifdef _PNCDF
      if(mpi_rnk==0) then
#endif
         e=nf90_put_var(f,     r, mes_oc%r(1:i_N,1))
         e=nf90_put_var(f,   icr, mes_ic%r(1:i_Nic,1))
         e=nf90_put_var(f,    dt, tim_dt)
         e=nf90_put_var(f, omega, rot_omega)
      end if

      call io_save_coll(f, uT,i_N,vel_Tor)
      call io_save_coll(f, uP,i_N,vel_Pol)
      call io_save_coll(f,  C,i_N,cod_C)
      call io_save_coll(f,  Comp,i_N,comp_C)
      call io_save_coll(f, BT,i_N,mag_Tor)
      call io_save_coll(f, BP,i_N,mag_Pol)
      call io_save_coll(f, icBT,i_Nic,mag_icTor)
      call io_save_coll(f, icBP,i_Nic,mag_icPol)

#ifndef _PNCDF
      if(mpi_rnk==0)  &
#endif
         e=nf90_close(f)

   end subroutine io_save_state
 


!--------------------------------------------------------------------------
!>  Define variable in netcdf file, ready for saving
!--------------------------------------------------------------------------
   subroutine io_define_coll(f,nm,dims,a, id)
      use netcdf
      use parameters, only : i_L, i_M, i_Mp
      integer,      intent(in) :: f !< netcdf file id
      integer,      intent(in) :: dims(3) !< array containing dimensions of variable
      character(*), intent(in) :: nm !< name of variable
      type (coll),  intent(in) :: a !< variable
      integer, intent(out) :: id !< netcdf variable id
      integer :: e
      e=nf90_def_var(f, nm, nf90_double, dims, id)
      e=nf90_put_att(f, id,  'L', i_L)      
      e=nf90_put_att(f, id,  'M', i_M)
      e=nf90_put_att(f, id, 'Mp', i_Mp)
   end subroutine io_define_coll

!--------------------------------------------------------------------------
!>  Save coll variable in netcdf file
!--------------------------------------------------------------------------
   subroutine io_save_coll(f,id,N__,a)
      use netcdf
      use parameters, only : i_pH1, i_N, i_L, i_M, i_Mp
      use dyn_mpi
      integer,     intent(in) :: f !< netcdf file id
      integer,     intent(in) :: id !< netcdf variable id
      integer,     intent(in) :: N__ !< number of radial points
      type (coll), intent(in) :: a !< scalar field to save
      integer :: e, n,j,nh_,nh__
      integer :: r,rr, pH0,pH1
      integer :: nh, l(0:i_pH1), m_(0:i_pH1)

#ifdef _PNCDF
      logical :: nf90coll
      nf90coll = .true.
      call mpi_barrier(mpi_comm_world, mpi_er)
      e=nf90_var_par_access(f,id, nf90_collective)
      if(e/=nf90_noerr)  print*, nf90_strerror(e)
      if(e/=nf90_noerr)  stop
      pH0 = a%H%pH0
      pH1 = a%H%pH1
      ct%Re = a%Re
      ct%Im = a%Im
      l (0:pH1) = a%H%l(pH0:pH0+pH1)
      m_(0:pH1) = a%H%m(pH0:pH0+pH1)
#else
      if(mpi_rnk==0) then
         do r = 0, mpi_sze-1
            rr  = modulo(r,_Nr)
            pH0 = a%H%pH0_(rr)
            pH1 = a%H%pH1_(rr)
!           print*,' TEST save_coll 1',r,pH0,pH1
            if(r==0) then
               ct%Re = a%Re
               ct%Im = a%Im
               l  = a%H%l(0:i_pH1)
               m_ = a%H%m(0:i_pH1)
!              print*,' TEST save_coll 2',r,pH0,pH1,l(7),m_(7) 
# ifdef _MPI
            else
               mpi_tg = r
               call mpi_recv( ct%Re, i_N*(pH1+1), mpi_double_precision,  &
                  r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
               call mpi_recv( ct%Im, i_N*(pH1+1), mpi_double_precision,  &
                  r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
               call mpi_recv( l,  pH1+1, mpi_integer,  &
                  r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
               call mpi_recv( m_, pH1+1, mpi_integer,  &
                  r, mpi_tg, mpi_comm_world, mpi_st, mpi_er)
# endif _MPI
            end if
#endif _PNCDF
            n = 1
            do nh = 0, pH1
               call io_state_nh(l(nh),m_(nh),i_L,i_M,i_Mp, nh__)
!              if (nh==7) print*,' TEST save_coll 3',l(7),m_(7),nh__,mpi_rnk
               if(nh>=pH1) nh_ = -1
               if(nh <pH1)  &
                  call io_state_nh(l(nh+1),m_(nh+1),i_L,i_M,i_Mp, nh_)
!                 if (nh==7) print*,' TEST save_coll 3B ',l(8),m_(8),nh_,mpi_rnk  
               if(nh__+1==nh_) then  ! mode contiguous with next file
                  n = n + 1
               else                  ! last n modes ctgs, save them!
                  j = nh__+1 - (n-1)
!                 if (mpi_rnk==0) print*, 'TEST save_coll 4', j, nh__, n
                  e=nf90_put_var(f,id, ct%Re(1:N__,nh-n+1:nh),start=(/1,j,1/))
                  e=nf90_put_var(f,id, ct%Im(1:N__,nh-n+1:nh),start=(/1,j,2/))
                  n = 1
#ifdef _PNCDF
                  if(nf90coll) e=nf90_var_par_access(f,id, nf90_independent)
                  nf90coll = .false.
#endif _PNCDF
               end if
            end do
#ifndef _PNCDF
         end do
# ifdef _MPI
      else
         mpi_tg = mpi_rnk
         call mpi_send( a%Re, i_N*(a%H%pH1+1), mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%Im, i_N*(a%H%pH1+1), mpi_double_precision,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%H%l(a%H%pH0), a%H%pH1+1, mpi_integer,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
         call mpi_send( a%H%m(a%H%pH0), a%H%pH1+1, mpi_integer,  &
            0, mpi_tg, mpi_comm_world, mpi_er)
# endif _MPI
      end if
#endif _PNCDF
   end subroutine io_save_coll


   subroutine io_write_vel_energy(mes_oc,ifile,ipfile,Tor,Pol,Ent,Efac,Dfac,tim_step,cod_C)
      use meshs, only : rdom
      use parameters, only : i_KL, i_L1, i_M1, b_boussinesq, d_PI, d_Pr
      use dyn_mpi, only : mpi_rnk
      use timestep, only : tim_t
      use variables, only : var_coll_normlm, var_coll_norm_comp, var_coll_norm_surf,&
                            var_coll_norm_veldiss, var_coll_qstcurl, var_coll_mom_comp,&
                            var_coll_norm_zon, var_coll_norm, var_coll_TorPol2qst, var_coll_normS

      type(rdom), intent(in) :: mes_oc !< outer core domain
      integer, intent(in) :: ifile !< file unit identifier for _energy.dat
      integer, intent(in) :: ipfile !< file unit identifier for _partition.dat
      type (coll), intent(in) :: Tor !< toroidal component of field
      type (coll), intent(in) :: Pol !< poloidal component of field
      type (coll), intent(in) :: Ent !< entropy/codensity (anelastic/boussinesq)
      double precision, intent(in) :: Efac !< prefactor for vel energy
      double precision, intent(in) :: Dfac !< prefactor for viscous dissipation
      integer, intent(in) :: tim_step !< time step
      type(coll), intent(in) :: cod_C !< codensity

      integer :: N_,pH0,pH1,L1,M1,pnh,nh,l,m,m_
      integer :: N,j, i
      double precision :: Ee,Eo,E1,E0,EDis,flx,dsi,dso,E, Esym, ETm0, El1
      double precision :: Nutop, Nubot,Dz0,c1,SS,topflx,botflx,Lx,Ly,Lz
      double precision :: Et,Etas,Etes,Etesas,ESin,ESout,ESin1,ESout1
      double precision :: Ep,Epas,Epes,Epesas,Eaxtot, Esurf, Ecomp
      ! arrays for calculating lpol and weighted sum
      double precision :: EQh(0:2*i_KL), EQl(0:i_L1), EQm(0:i_M1)
      double precision :: ETh(0:2*i_KL), ETl(0:i_L1), ETm(0:i_M1)
      double precision :: ESh(0:2*i_KL), ESl(0:i_L1), ESm(0:i_M1)
      double precision :: lpol, l_Eofl
      logical :: calc_length = .True.


      if(mpi_rnk==0) then
        N = cod_C%D%N
        flx = 0d0
        do j = N-i_KL, N
          flx = flx  -  cod_C%Re(j,0) * cod_C%D%dr(1)%M(i_KL+1+N-j,j)
        end do
! #CAJ 15th October gets topflx and botflx either anelastic or Boussinesq
!
        if (.not.b_boussinesq) then
          topflx=4.0d0*d_PI*cod_C%D%r(N,2)*cod_C%D%rho(N,0)*cod_C%D%Temp0(N,0)*flx
        else
          topflx=4.0d0*d_PI*cod_C%D%r(N,2)*flx
        end if
        flx = 0d0
        do j = 1,1+i_KL
           flx = flx  -  cod_C%Re(j,0) * cod_C%D%dr(1)%M(i_KL+2-j,j)
        end do
        if (.not.b_boussinesq) then
          botflx=4.0d0*d_PI*cod_C%D%r(1,2)*cod_C%D%rho(1,0)*cod_C%D%Temp0(1,0)*flx
        else
          botflx=4.0d0*d_PI*cod_C%D%r(1,2)*flx
        end if

! #CAJ 28th March replace xi with rho(:,1)
      end if

      call var_coll_TorPol2qst(Tor,Pol, cq,cs,ct)

      ! Calculating l-dependent lengthscales
      l_Eofl = 0
      lpol = 0
      if (calc_length) then
         call var_coll_normlm(ct, ETh,ETl,ETm)
         call var_coll_normlm(cq, EQh,EQl,EQm)
         call var_coll_normlm(cs, ESh,ESl,ESm)
         do i = 1, i_L1
            l_Eofl = l_Eofl + i * Efac*(ETl(i)+EQl(i)+ESl(i))
         end do
         lpol=maxloc(EQl+ESl,1)
      endif

      if (.not.b_boussinesq) then
        call var_coll_norm_comp(cq, Eo,Ee,E1,E0)
      else
        call var_coll_norm(cq, Eo,Ee,E1,E0)
      end if
      call var_coll_norm_surf(cq,ESin,ESout)
      if (.not.b_boussinesq) then 
        call var_coll_norm_veldiss(cq,Ecomp)
      else 
        Ecomp=0d0
      end if
      Ep     = Eo + Ee
      Epes   = Ee
      Epesas = E1
      Epas   = E0
      Esym   = Ee
      El1    = E1
      if (.not.b_boussinesq) then
        call var_coll_norm_comp(cs, Eo,Ee,E1,E0)
      else
         call var_coll_norm(cs, Eo,Ee,E1,E0)
      end if
      call var_coll_norm_surf(cs,ESin1,ESout1)
      Ep     = Ep   + Eo + Ee
      Epes   = Epes + Ee
      Epesas = Epesas   + E1
      Epas   = Epas + E0
      ESin   = ESin + ESin1
      ESout  = ESout + ESout1
      Esym   = Esym + Ee
      El1    = El1 + E1
      if (.not.b_boussinesq) then
        call var_coll_norm_comp(ct, Eo,Ee,E1,E0)
      else
        call var_coll_norm(ct, Eo,Ee,E1,E0)
      end if
      call var_coll_norm_surf(ct,ESin1,ESout1)
      Et     = Eo + Ee
      Etes   = Eo
      Etesas = E0-E1
      Etas   = E0
      ESin   = ESin + ESin1
      ESout  = ESout + ESout1
      Esym   = Esym + Eo
      El1    = EL1 + E1
      call var_coll_normS(Ent, SS)
      SS = SS * d_Pr
      call var_coll_qstcurl(cq,cs,ct, cq,cs,ct)
      if (.not.b_boussinesq) then 
        call var_coll_norm_comp(cq, Eo,Ee,E1,E0)
      else
        call var_coll_norm(cq, Eo,Ee,E1,E0)
      end if
      EDis = Eo + Ee
      if (.not.b_boussinesq) then
        call var_coll_norm_comp(cs, Eo,Ee,E1,E0)
      else
        call var_coll_norm(cs, Eo,Ee,E1,E0)
      end if
      EDis = EDis + Eo + Ee
      if (.not.b_boussinesq) then
        call var_coll_norm_comp(ct, Eo,Ee,E1,E0)
      else
        call var_coll_norm(ct, Eo,Ee,E1,E0)
      end if 
      EDis = EDis + Eo + Ee
      ETm0   = E0
      ETm0   = Efac * ETm0
!     E    = Efac * (Ep+Et) * d_Pr * c1 / d_Ra
      E      = Efac * (Ep+Et)
      Et     = Efac * Et
      Etas   = Efac * Etas
      Etes   = Efac * Etes
      Etesas = Efac * Etesas
      Ep     = Efac * Ep
      Epas   = Efac * Epas
      Epes   = Efac * Epes
      Epesas = Efac * Epesas
      Eaxtot = Epas + Etas
      EDis   = Dfac * EDis
      Esurf  = Dfac * 2d0 * (ESin - ESout)
      EDis   = EDis  + Esurf + Ecomp
      Esym   = Efac * Esym
      El1    = Efac * El1
      ! Saving value to compute Rm, Re later
      io_KinEn  = E
      io_ViscDiss  = EDis

! ## Piotr ## (28/07/2009)
! Added computation of angular momentum
! Why? It isn't written out anywhere!!!!
      call var_coll_mom_comp(mes_oc,Tor, Lx,Ly,Lz)
      call var_coll_norm_zon(Tor,E0)
      if(mpi_rnk/=0) return
      if (tim_step==0)  &
      write(ifile,'(12A20)') '#       Time        ','   Kinetic Energy  ','   Total Toroidal  ','   E_eq_sym   ', 'E_l=1', &
                             '  E_tor_m=0 ', 'Edis      ', 'Axisymmetric KE'
      write(ifile,'(12e20.12)') tim_t, E, Et, Esym, El1, ETm0, Edis, Eaxtot


      if (tim_step==0)  &

      write(ipfile,'(11A20)') '#       Time       ',' Total toroidal    ',' Total poloidal    ','  axisym. tor.     ', &
                             ' axisym. pol.      ',' equat.sym.tor.    ', ' equat.sym.pol.   ', &
                             'equat+axisym.tor.  ','equat.axisym.pol.  ', ' sum(l x u(l).u(l)) ', ' lpol '

      write(ipfile,'(11e20.12)') tim_t, Et, Ep, Etas, Epas, &
                                Etes,Epes,Etesas,Epesas, l_Eofl, lpol

! flushing
      CALL flush(ipfile)
      CALL flush(ifile)
   end subroutine io_write_vel_energy

   
   
!--------------------------------------------------------------------------
!>  write bulk magnetic energy to energy file
!--------------------------------------------------------------------------
   subroutine io_write_mag_energy(mes_oc,ifile,ipfile,Tor,Pol,Efac,Dfac,mag_tor,mag_pol)
      use meshs, only : rdom
      use dyn_mpi, only : mpi_rnk
      use variables, only : var_coll_TorPol2qst_mag, var_coll_norm, var_coll_qstcurl, var_coll_norm_magdiss
      use timestep, only : tim_t, tim_step

      type(rdom), intent(in) :: mes_oc !< outer core domain
      integer, intent(in) :: ifile !< unit file identifier for _energy.dat
      integer, intent(in) :: ipfile !< unit file identifier for _partition.dat
      type (coll), intent(in) :: Tor !< toroidal component of field
      type (coll), intent(in) :: Pol !< poloidal component of field 
      type (coll), intent(in) :: mag_tor !< toroidal component of mag field (outer core)
      type (coll), intent(in) :: mag_pol !< poloidal component of mag field (outer core)
      double precision, intent(in) :: Efac !< prefactor for magnetic energy
      double precision, intent(in) :: Dfac !< prefactor for ohmic dissipation

      double precision :: Ee,Eo,E1,E0, E,ETor,ESym,El1,ETm0, EDis,Dip,Quad,Bth,Br,Emer
      double precision :: Ep,Et,Epes,Etes,Epesas,Etesas,Epas,Etas,Eaxtot
      integer :: N_

      if(mpi_rnk==0) then
       N_ = mes_oc%N
       Dip = mag_Pol%Re(N_,1)/mes_oc%r(N_,1)
       Quad = 2*mag_Pol%Re(N_,2)/mes_oc%r(N_,1)
      endif

      call var_coll_TorPol2qst_mag(Tor,Pol, cq,cs,ct)
      call var_coll_norm(cq, Eo,Ee,E1,E0)
! Modify to split into the poloidal, toroidal, axisym poloidal, axisym toroidal,
! Equatorially symmetric poloidal (l+m even), Equatorially symmetric toroidal
! (l+m eveni), Equ symm and axisymm pol, Equ symm and axisymm tor.
      Ep    = Eo + Ee
      Epes  = Ee
      Epesas= E1
      Epas  = E0
      Esym  = Ee
      El1   = E1
      call var_coll_norm(cs, Eo,Ee,E1,E0)
      Ep    = Ep + Eo + Ee
      Epes  = Epes + Ee
      Epesas= Epesas + E1
      Epas  = Epas + E0
      Esym  = Esym + Ee
      El1   = El1 + E1
      call var_coll_norm(ct, Eo,Ee,E1,E0)
      Et    = Eo + Ee
      Etes  = Eo
      Etesas= E0-E1
      Etas  = E0
      ETm0  = E0
      Esym  = Esym + Eo
      El1   = El1 + E1

      call var_coll_qstcurl(cq,cs,ct, cq,cs,ct)
      call var_coll_norm_magdiss(cq, Eo,Ee,E1,E0)
      EDis = Eo + Ee
      call var_coll_norm_magdiss(cs, Eo,Ee,E1,E0)
      EDis = EDis + Eo + Ee
      call var_coll_norm_magdiss(ct, Eo,Ee,E1,E0)
      EDis = EDis + Eo + Ee
      
      E      = Efac * (Ep + Et)
      Et     = Efac * Et
      Etas   = Efac * Etas
      Etes   = Efac * Etes
      Etesas = Efac * Etesas
      Ep     = Efac * Ep
      Epas   = Efac * Epas
      Epes   = Efac * Epes
      Epesas = Efac * Epesas
      Eaxtot = Etas + Epas
      EDis   = Dfac * EDis
      ETm0   = Efac * ETm0
      Esym   = Efac * Esym
      El1    = Efac * El1
      ! Saving values to compute Elsasser, fohm
      io_MagEn  = E
      io_OhmDiss  = EDis

      if(mpi_rnk/=0) return

      if (tim_step==0)  &
      write(ifile,'(10A22)') '# Time ','Magnetic Energy ', 'Total Toroidal ', 'E_eq_sym ', 'E_l=1 ', 'E_tor_m=0 ', 'Mag. Dissipation ', 'Axial Dipole g10 ','Axial Quadrupole g20 ','Axisymmetric ME '
      write(ifile,'(10e22.12)')  tim_t, E, Et, Esym, El1, ETm0, EDis, Dip, Quad, Eaxtot

      if (tim_step==0)  &
      write(ipfile,'(9A20)') '# Time ','Total toroidal ', 'Total poloidal ','axisym. tor. ', 'axisym. pol. ','equat.sym.tor. ', ' equat.sym.pol. ', 'equat+axisym.tor. ', 'equat.axisym.pol. '
      write(ipfile,'(9e20.12)') tim_t, Et, Ep, Etas, Epas, Etes,Epes,Etesas,Epesas

! flushing
      CALL flush(ifile)
   end subroutine io_write_mag_energy


!--------------------------------------------------------------------------
!> Write magnetic energy at CMB
!!
!! Top surface magnetic energies go are written to mag_cmb.dat
!! 7 entries: time, surface magnetic energy, Sum of l+m even contributions to
!! the surface magnetic energy, l=1 contributions to the surface magnetic energy,
!! theta location of the dipole term N magnetic pole, phi location of the 
!! dipole term N magnetic pole, m=0 contributions to the surface magnetic energy
!--------------------------------------------------------------------------
   subroutine io_write_cmb(mes_oc,ifile,Tor,Pol,Efac)
      use meshs, only : rdom
      use dyn_mpi
      use parameters, only : i_Mp, i_M1, d_PI
      use timestep, only : tim_t, tim_step
      use variables, only : var_coll_norm_orig, var_coll_TorPol2qst

      type(rdom), intent(in) :: mes_oc !< outer core domain
      integer, intent(in) :: ifile !< file unit identifier for mag_cmb.dat
      type (coll), intent(in) :: Tor !< toroidal part of magnetic field
      type (coll), intent(in) :: Pol !< poloidal part of magnetic field
      double precision, intent(in) :: Efac !< prefactor for magnetic energy

      double precision :: Eo,Ee,E1,E0, EPol,ESym,El1, g10,g11,h11,th,ph
      double precision :: Epas, E12, El12
      integer :: N_, L_, r,rt,nh_
      _loop_lm_vars
      
      call var_coll_TorPol2qst(Tor,Pol, cq,cs,ct)
      N_ = mes_oc%N
      ! Set all values except cmb to zero
      cq%Re(:N_-1,:) = 0d0
      cq%Im(:N_-1,:) = 0d0
      cs%Re(:N_-1,:) = 0d0
      cs%Im(:N_-1,:) = 0d0
      call var_coll_norm_orig(cq, Eo,Ee,E1,E0,E12)
      EPol = Eo + Ee
      ESym = Ee
      El1  = E1
      Epas   = E0
      El12 = E12
      
      call var_coll_norm_orig(cs, Eo,Ee,E1,E0,E12)
      EPol = EPol + Eo + Ee
      Epas = Epas + E0
      ESym = ESym + Ee
      El1  = El1  + E1
      El12 = El12 + E12

      EPol = Efac * EPol * mes_oc%r(N_,2) / mes_oc%intr2dr(N_)
      Epas = Efac * Epas * mes_oc%r(N_,2) / mes_oc%intr2dr(N_)  
      ESym = Efac * ESym * mes_oc%r(N_,2) / mes_oc%intr2dr(N_)
      El1  = Efac * El1  * mes_oc%r(N_,2) / mes_oc%intr2dr(N_)
      El12  = Efac * El12  * mes_oc%r(N_,2) / mes_oc%intr2dr(N_)
      
      r  = 0
      nh_= 0
      _loop_lm_begin(Pol)
         if(l==1 .and. m_==0) then
            r  = mpi_rnk
            nh_= nh
         end if
      _loop_lm_end
      g10 = Pol%Re(N_,nh_)
#ifdef _MPI
      call mpi_allreduce(r, rt, 1, mpi_integer,  &
         mpi_sum, mpi_comm_world, mpi_er)
      call mpi_bcast(g10, 1, mpi_double_precision,  &
         rt, mpi_comm_world, mpi_er)
#endif
      
      if(i_Mp==1 .and. i_M1>0) then
         r  = 0
         nh_= 0
         _loop_lm_begin(Pol)
            if(l==1 .and. m_==1) then
               r  = mpi_rnk
               nh_= nh
            end if
         _loop_lm_end
         g11 = -2d0*Pol%Re(N_,nh_)
         h11 =  2d0*Pol%Im(N_,nh_)
#ifdef _MPI
         call mpi_allreduce(r, rt, 1, mpi_integer,  &
            mpi_sum, mpi_comm_world, mpi_er)
         call mpi_bcast(g11, 1, mpi_double_precision,  &
            rt, mpi_comm_world, mpi_er)
         call mpi_bcast(h11, 1, mpi_double_precision,  &
            rt, mpi_comm_world, mpi_er)
#endif            
         th = dsqrt(g10*g10 + g11*g11 + h11*h11)
         if(th/=0d0)  th = dacos( g10 / max(th,dabs(g10)) )
         if(dabs(g11)<1d-99) then
            if(h11>=0d0) ph = 0.5d0*d_Pi
            if(h11 <0d0) ph = 1.5d0*d_Pi
         else
            ph = datan( h11/g11 )
            if(g11<0d0) ph = ph + d_Pi
            if(ph <0d0) ph = ph + 2d0*d_Pi
         endif
      else
         if(g10>=0d0) th = 0d0
         if(g10< 0d0) th = d_Pi
         ph = 0d0
      end if

      ! saving values for writing non-dimensional outputs
      io_fdip = SQRT(El1/El12)
      io_thdip = th
          
      if(mpi_rnk/=0) return
      if (tim_step==0)  &
      write(ifile,'(8A20)') '#      Time         ','    Surf mag Energy  ', &
      '   Even mag Energy   ',' l=1 surf mag Energy ','   N Mag pole theta  ', &
                             '   N Mag pole phi    ',' Axisym Surf mag en.', &
                             '  l<=12 surf mag Energy '
      write(ifile,'(8e20.12)') tim_t, EPol, ESym, El1, th, ph,Epas, El12
! flushing
      CALL flush(ifile)
   end subroutine io_write_cmb


!--------------------------------------------------------------------------
!> Write nusselt number
!!
!!  write the Nusselt number:  actual flux / flux for no flow  at cmb
!!  Only used if b_boussinesq==.true., i_cod_bc=1, b_inhombc_c==.false.
!! Otherwise, definition of Nusselt number is unclear. For anelastic cases,
!! the reference state flux is best found in the matlab script that generates 
!! the reference state, and then if well-defined the Nusselt number can be 
!! found from the flux output in vel_energy.dat 
!--------------------------------------------------------------------------
   subroutine io_write_nussel(mes_oc,ifile, buoyancy_field, field_type)
      use dyn_mpi, only : mpi_rnk
      use parameters, only : i_N, d_rratio, i_cod_bc
      use timestep, only : tim_step, tim_t
      use meshs, only : rdom
      use variables, only : var_coll_meshmult

      type(rdom), intent(in) :: mes_oc !< outer core domain
      integer,     intent(in) :: ifile !< file unit identifier
      integer,     intent(in) :: field_type !< 1=codensity, 2=composition
      type (coll), intent(in) :: buoyancy_field !< codensity/composition
      double precision :: dCdr_ro, dCdr_ri, T_ro, T_ri
      double precision :: ri, ro
      double precision :: A, DT, DT_cond, nu_top, nu_bot, nu_flux, nu

      type (coll) :: dCdr

      integer :: N, j, n1

      dCdr_ro=0
      dCdr_ri=0

      ! on a unit sphere:
      ! int_0^2pi int_0^pi A sin theta dtheta dphi = 4pi A_0^0,
      ! this means that T is just cod_c%re(r_index,H=0), where h=0
      ! corresponds to l=0,m=0. We ignore r dependence and the factor of 
      ! 4pi as we want the surface average, not the integral.
      if (mpi_rnk==0) then
         T_ri = buoyancy_field%Re(1,0)
         T_ro = buoyancy_field%Re(i_N,0)
      end if
      
      ! get the codensity gradient in the dCdr coll variable
      call var_coll_meshmult(mes_oc%dr(1),buoyancy_field, dCdr)
      
      ! Use the H=0 (l=0, m=0) at r=ri, r=ro to find the flux
      ! on the inner and outer boundaries.
      if (mpi_rnk==0) then
         dCdr_ri = dCdr%Re(1,0)
         dCdr_ro = dCdr%Re(i_N,0)
      end if
      
      if(mpi_rnk/=0) return
      if(buoyancy_field%H%pH0>0) stop 'io_write_nussel'
            
      ! Find nusselt number for simple cases. This may not be correct
      ! if complex boundary conditions are in play
      ! This is a prototype, want to create a precompute subroutine
      ! that will calculated the conduction profile / Delta T for any
      ! given set of boundary conditions / source terms, and then
      ! use those precomputed conduction terms to compute the nusselt
      ! number from actual Delta T / dT/dr.
      
      ri = d_rratio/(1-d_rratio)
      ro = 1/(1-d_rratio)
      A = (ri*ro)/(ri-ro)
      Nu_top = 0
      Nu_bot = 0
      nu_flux = 0
      ! 
      if (i_cod_bc==1) then
      ! conduction profile for Delta T = 1 is dT/dr = A/r^2, where
      ! A = (ri*ro)/(ri-ro)
         nu_bot= dCdr_ri / (A/ri**2)
         nu_top = dCdr_ro / (A/ro**2)
         nu = nu_top 
      else if (i_cod_bc==4) then
      ! conduction profile when setting A=1 gives Delta T = 1/ri - 1/ro
         DT_cond =  1/(ri) - (1/ro)
         DT = T_ri - T_ro
         nu_flux = DT_cond/DT
         nu=nu_flux
      end if
      ! saving nusselt number to write out elsewhere
      if (field_type==1) io_Nus = nu
      if (field_type==2) io_Nus_comp = nu

      if (tim_step==0) & 
         write(ifile, '(9A20)') '#      Time         ','    Nusselt        ',&
         '    T_top    ','   T_bot   ', &
         '    dCdr_ri   ','   dCdr_ro   ',' Nu_top ',' Nu_bot ',' Nu_flux '

      write(ifile,'(9e20.12)') tim_t, nu, T_ro, T_ri,&
         dCdr_ri, dCdr_ro, nu_top, nu_bot, nu_flux
      CALL flush(ifile)
   end subroutine io_write_nussel


!--------------------------------------------------------------------------
!>  write to torque file
!--------------------------------------------------------------------------
   subroutine io_write_torque(ifile, rot_omega, rot_velTorq, rot_magTorq)
      use dyn_mpi, only : mpi_rnk
      use timestep, only : tim_t, tim_step

      integer, intent(in) :: ifile !< file unit identifier rot_torque.dat
      double precision, intent(in) :: rot_omega !< inner core rotation rate
      double precision, intent(in) :: rot_velTorq !< viscous torque
      double precision, intent(in) :: rot_magTorq !< magnetic torque

      if(mpi_rnk/=0) return
      if (tim_step==0)  &
!  K.K. added "hash" character to mark comment
      write(ifile,'(12A20)') '#       Time        ',' Inner core rotation',' &
&Viscous torque ','  Magnetic torque '

      write(ifile,'(4e20.12)') tim_t, rot_omega, rot_velTorq, rot_magTorq
      CALL flush(ifile)
   end subroutine io_write_torque


!--------------------------------------------------------------------------
!> Write output for 2011 Anelastic Benchmard
!!
!!  This routine is not currently used. It was used to generate outputs
!!  required in the 2011 anelastic benchmark paper. It is retained here
!!  as it could be easily modified to print out quantities that users
!!  might be interested in
!--------------------------------------------------------------------------
   subroutine io_write_bench(mes_oc,ifile, vel_pol, vel_r, vel_p, mag_t, cod_C)
      use meshs, only : rdom
      use variables, only : phys, spec, var_Th, var_coll2spec
      use parameters, only : i_Th
      use dyn_mpi
      use timestep, only : tim_t
      use transform, only : tra_spec2phys

      type(rdom), intent(in) :: mes_oc !< outer core domain
      integer, intent(in) :: ifile !< file unit identifier
      type(coll), intent(in) :: vel_pol !< poloidal component of velocity
      type(phys), intent(in) :: vel_r !< radial component of velocity
      type(phys), intent(in) :: vel_p !< phi component of velocity
      type(phys), intent(in) :: mag_t !< theta component of magnetic field
      type(coll), intent(in) :: cod_C !< codensity

      double precision :: dat(6)
      double precision, save :: t=-1d0,t_
      double complex,   save :: c_, c=0d0
      integer :: Nr,Nth, r,r_,rr,rth
      type (spec) :: sC
      type (phys) :: pC
      
      call var_coll2spec(cod_C, sC)
      call tra_spec2phys(sC, pC)
      call var_coll2spec(vel_Pol, sC)
      Nr  = (mes_oc%N+1)/2
      Nth = (i_Th+1)/2
      do r = 0, mpi_sze-1
         rr  = modulo(r,_Nr)
         rth = r/_Nr
         if(mes_oc%pNi_(rr)<=Nr .and.  &
            mes_oc%pNi_(rr)+mes_oc%pN_(rr)-1>=Nr .and.  &
            var_Th%pThi_(rth)<=Nth .and.  &
            var_Th%pThi_(rth)+var_Th%pTh_(rth)-1>=Nth )  &
            r_ = r
      end do
      t_   = t
      t    = tim_t
      if(r_==mpi_rnk) then
         Nr = Nr - mes_oc%pNi + 1
         Nth= Nth- var_Th%pThi + 1
         dat(1) = vel_r%Re(Nth,0,Nr)	! u_r
         dat(2) = vel_p%Re(Nth,0,Nr)	! u_ph
         dat(3) =    pC%Re(Nth,0,Nr)	! C
         dat(4) = mag_t%Re(Nth,0,Nr)	! B_th
!         c_   = c
!         c    = cmplx( sC%Re(i_L,Nr), sC%Im(i_L,Nr) )
!         if(abs(c_)/=0d0)  dat(5) =  real(log(c/c_))/(t-t_)  ! growth
!         if(abs(c_)/=0d0)  dat(6) = aimag(log(c/c_))/(t-t_)  ! drift*Mp
      end if
      if(t_==-1d0) return
#ifdef _MPI
      call mpi_bcast(dat, 6, mpi_double_precision,  &
                     r_, mpi_comm_world, mpi_er)
#endif
      if(mpi_rnk/=0) return
      write(ifile,'(7e20.12)') tim_t, (dat(r), r=1,4) !6)
   end subroutine io_write_bench


!--------------------------------------------------------------------------
!>  write to timestep file
!--------------------------------------------------------------------------
   subroutine io_write_timestep(ifile)
      use dyn_mpi, only : mpi_rnk
      use timestep, only : tim_t, tim_dt, tim_corr_dt, tim_cfl_dt
      integer, intent(in) :: ifile !< file unit identifier
      if(mpi_rnk/=0) return
      write(ifile,'(4e20.12)') tim_t, tim_dt, tim_corr_dt, tim_cfl_dt
   end subroutine io_write_timestep

!--------------------------------------------------------------------------
!>  write to non-dimensional outputs file
!!
!! quantities calculated from variables saved
!! in previously run io subroutines, such as io_write_vel_energy
!--------------------------------------------------------------------------
   subroutine io_write_dyn_outputs(ifile)
      use dyn_mpi, only : mpi_rnk
      use parameters, only : b_vel_tstep, b_mag_tstep, b_cod_tstep, b_comp_tstep,&
                             d_Pm, d_Pr, d_Sc, d_E
      use timestep, only : tim_step, tim_t

      integer, intent(in) :: ifile !< file unit identifier

      double precision :: Rm, Re, Pe, Pe_comp, El, Le, fohm, M, Rossby, ri, ro
      if(mpi_rnk/=0) return
      if(.not. b_vel_tstep) then
         Rm=0
         Re=0
         Pe=0
         Pe_comp=0
         Rossby = 0
      else
         Rm = SQRT((2*io_KinEn)/io_vol)
         Re = Rm / d_Pm
         Pe = Re * d_Pr
         Pe_comp = Re * d_Sc
         Rossby = Re * d_E
      endif
      if(.not. b_mag_tstep) then
         El = 0
         Le = 0
      else
         El = (d_E * io_MagEn)/(d_Pm * io_vol)
         Le = SQRT((2*d_E**2 * io_MagEn)/(d_Pm**2 * io_vol))
      endif
      if ( (b_vel_tstep) .and. (b_mag_tstep)) then
         M = io_MagEn / io_KinEn
         fohm = io_OhmDiss/(io_OhmDiss + io_ViscDiss)
      else
         M=0
         fohm = 0
      endif
      if (.not. b_cod_tstep) io_nus = 0
      if (.not. b_comp_tstep) io_Nus_comp = 0
      if (tim_step==0) &
      write(ifile,'(A1,A17,13A19)') '#', 'time', 'Rm', 'Reynolds', 'Pe', 'Pe_comp', 'Rossby', 'M_ratio', 'Elsasser', 'Lehnert', &
         'fdip', 'fohm', 'theta_dip', 'Nusselt', 'Sherwood'
      write(ifile,'(es18.11,13es19.11)') tim_t, Rm, Re, Pe, Pe_comp, Rossby, M, El, Le, io_fdip, fohm, io_thdip, io_Nus, io_Nus_comp
      CALL flush(ifile)
   end subroutine io_write_dyn_outputs


!--------------------------------------------------------------------------
!>  Save spectrum
!--------------------------------------------------------------------------
   subroutine io_save_spectrum(prefix,Tor,Pol,Efac)
      use parameters, only : i_KL, i_L1, i_M1, i_Mp
      use dyn_mpi, only : mpi_rnk
      use timestep, only : tim_t
      use variables, only : var_coll_TorPol2qst_mag, var_coll_normlm
      character(4) :: cnum

      character(3), intent(in) :: prefix !< 3 letter id for field, eg vel, mag
      type (coll), intent(in) :: Tor !< toroidal part of field
      type (coll), intent(in) :: Pol !< poloidal part of field
      double precision, intent(in) :: Efac !< prefactor for field energy

      double precision :: EQh(0:2*i_KL), EQl(0:i_L1), EQm(0:i_M1)
      double precision :: ETh(0:2*i_KL), ETl(0:i_L1), ETm(0:i_M1)
      double precision :: ESh(0:2*i_KL), ESl(0:i_L1), ESm(0:i_M1),lsum,msum
      integer :: i
   10 format(i4,1e20.12)
   11 format(1e20.12)

      write(cnum,'(I4)') io_save3
      if(io_save3<1000) cnum(1:1) = '0'
      if(io_save3< 100) cnum(2:2) = '0'
      if(io_save3<  10) cnum(3:3) = '0'
      if(mpi_rnk==0)  &
         open(11, status='unknown', file=prefix//'_spectrum_'//cnum//'.dat')

      call var_coll_TorPol2qst_mag(Tor,Pol, cq,cs,ct)

      if(mpi_rnk==0)  &
         write(11,*) '# t = ', tim_t

      call var_coll_normlm(ct, ETh,ETl,ETm)
      call var_coll_normlm(cq, EQh,EQl,EQm)
      call var_coll_normlm(cs, ESh,ESl,ESm)
      !Polidal part made up of q, s components, (from qst decompostion)
      !Toroidal part made up of t component.
      !EPol(l) = Efac * (EQl + ESl)
      !ETor(l) = Efac * (ETl)
      !EPol(m) = Efac * (EQm + ESm)
      !ETor(m) = Efac * (ETm)
     if(mpi_rnk==0) then
         write(11,*) '# l, E(l), Epol(l), Etor(l)'
         do i = 1, i_L1
            write(11,'(i4,3e20.12)') i, Efac*(ETl(i)+EQl(i)+ESl(i)), Efac*(EQl(i)+ESl(i)), Efac*(ETl(i))

         end do
         lsum=0d0
         do i = 1,i_L1
         lsum=lsum+Efac*(ETl(i)+EQl(i)+ESl(i))
         end do

         write(11,*) 'Sum over l'
         write(11,11) lsum
         ! Two blank lines here for gnuplot indices
         write(11,*)
         write(11,*)
         write(11,*) '# m, E(m), Epol(m), Etor(m)'

         do i = 0, i_M1
            write(11,'(i4,3e20.12)') i*i_Mp, Efac*(ETm(i)+EQm(i)+ESm(i)), Efac*(EQm(i)+ESm(i)), Efac*(ETm(i))

         end do
         msum=0d0
         do i = 0,i_M1
            msum=msum+Efac*(ETm(i)+EQm(i)+ESm(i))
         end do
         write(11,*) 'Sum over m'
         write(11,11) msum
         CALL flush(11)
         close(11)
      end if

   end subroutine io_save_spectrum


!--------------------------------------------------------------------------
!>  save spectrum
!!
!! Used for vector quantities, magnetic and velocity fields
!--------------------------------------------------------------------------
   subroutine io_save_spectrum_anelastic(prefix,Tor,Pol,Efac)
      use parameters, only : i_KL, i_L1, i_M1, i_Mp
      use dyn_mpi, only : mpi_rnk
      use timestep, only : tim_t
      use variables, only : var_coll_TorPol2qst, var_coll_normlm_comp
      character(4) :: cnum
      
      character(3), intent(in) :: prefix !< 3 letter id for field, eg vel, mag
      type (coll), intent(in) :: Tor !< toroidal part of field 
      type (coll), intent(in) :: Pol !< poloidal part of field
      double precision, intent(in) :: Efac  !< prefactor for field energy

      double precision :: ETh(0:2*i_KL), ETl(0:i_L1), ETm(0:i_M1)
      double precision :: EPh(0:2*i_KL), EPl(0:i_L1), EPm(0:i_M1), lsum,msum
      integer :: i
   10 format(i4,1e20.12)
   11 format(1e20.12)

      write(cnum,'(I4)') io_save3
!     CALL flush(cnum)
      if(io_save3<1000) cnum(1:1) = '0'
      if(io_save3< 100) cnum(2:2) = '0'
      if(io_save3<  10) cnum(3:3) = '0'
      if(mpi_rnk==0)  &
         open(11, status='unknown', file=prefix//'_spectrum_'//cnum//'.dat')

      call var_coll_TorPol2qst(Tor,Pol,cq,cs,ct)

      if(mpi_rnk==0)  &
         write(11,*) '# t = ', tim_t
! flushing
!     CALL flush(11)

      call var_coll_normlm_comp(ct, ETh,ETl,ETm)
      call var_coll_normlm_comp(cq, EPh,EPl,EPm)
      ETl = ETl + EPl
      ETm = ETm + EPm
      call var_coll_normlm_comp(cs, EPh,EPl,EPm)
      ETl = Efac * (ETl + EPl)
      ETm = Efac * (ETm + EPm)

      if(mpi_rnk==0) then
         write(11,*) '# l'
! flushing
!        CALL flush(11)
         do i = 1, i_L1
            write(11,10) i, ETl(i)
! flushing
!           CALL flush(11)

         end do
         lsum=0d0
         do i = 1,i_L1
            lsum=lsum+ETl(i)
         end do
         write(11,*) 'Sum over l'
         write(11,11) lsum
         write(11,*)
         write(11,*) '# m'
! flushing
!        CALL flush(11)

         do i = 0, i_M1
            write(11,10) i*i_Mp, ETm(i)
! flushing
!           CALL flush(11)

         end do
         msum=0d0
         do i = 0,i_M1
            msum=msum+ETm(i)
         end do
         write(11,*) 'Sum over m'
         write(11,11) msum
             CALL flush(11)
         close(11)
      end if

   end subroutine io_save_spectrum_anelastic


!--------------------------------------------------------------------------
!>  save spectrum
!!
!! Used for scalar quantities, codensity and composition fields
!--------------------------------------------------------------------------
   subroutine io_save_spectrum1(prefix,a,Efac)
      use parameters, only : i_Kl, i_L1, i_M1, i_Mp
      use dyn_mpi, only : mpi_rnk
      use timestep, only : tim_t
      use variables, only : var_coll_normlm
      character(4) :: cnum

      character(3), intent(in) :: prefix !< 3 letter id for field: cod/com
      type (coll),  intent(in) :: a !< scalar field to be analysed
      double precision, intent(in) :: Efac !< prefactor for energy calc

      double precision :: Eh(0:2*i_KL), El(0:i_L1), Em(0:i_M1)
      integer :: i
   10 format(i4,1e20.12)

      write(cnum,'(I4)') io_save3
!     CALL flush(cnum)
      if(io_save3<1000) cnum(1:1) = '0'
      if(io_save3< 100) cnum(2:2) = '0'
      if(io_save3<  10) cnum(3:3) = '0'
      if(mpi_rnk==0)  &
         open(11, status='unknown', file=prefix//'_spectrum_'//cnum//'.dat')

      if(mpi_rnk==0)  &
         write(11,*) '# t = ', tim_t
! flushing
!     CALL flush(11)

      call var_coll_normlm(a, Eh,El,Em)
      El = Efac * El
      Em = Efac * Em

      if(mpi_rnk==0) then
         write(11,*) '# l'
! flushing
!     CALL flush(11)
         do i = 1, i_L1
            write(11,10) i, El(i)
! flushing
!     CALL flush(11)
         end do
         write(11,*)
         write(11,*) '# m'
! flushing
!     CALL flush(11)
         do i = 0, i_M1
            write(11,10) i*i_Mp, Em(i)
! flushing
!     CALL flush(11)
         end do
         close(11)
      end if

   end subroutine io_save_spectrum1

!--------------------------------------------------------------------------
!> write out radial temp profile
!!
!! Write out radius, mean temperature
!! and mean temperature gradient.
!--------------------------------------------------------------------------
   subroutine io_save_profile(mes_oc, label,field)
      use meshs, only : rdom
      use parameters, only : i_N
      use dyn_mpi, only : mpi_rnk
      use variables, only : var_coll_meshmult
      type(rdom), intent(in) :: mes_oc !< outer core domain
      character(4) :: label !< prefix id for field: temp/comp
      type (coll), intent(in) :: field !< scalar field to be analysed
      character(4) :: cnum
      type (coll) :: dTdr
      double precision :: T_r(i_N), dTdr_r(i_N), r(i_N)
      integer :: i

      write(cnum,'(I4)') io_save3
!     CALL flush(cnum)
      if(io_save3<1000) cnum(1:1) = '0'
      if(io_save3< 100) cnum(2:2) = '0'
      if(io_save3<  10) cnum(3:3) = '0'
      if(mpi_rnk==0)  &
         open(11, status='unknown', file=label//'_profile_'//cnum//'.dat')

      call var_coll_meshmult(mes_oc%dr(1),field,dTdr)

      if (mpi_rnk==0) then
         T_r = field%Re(:,0)
         dTdr_r = dTdr%Re(:,0)
      
         r = mes_oc%r(:,1)


         write(11,'(3A10)') '# r     ','   C   ','   dC/dr   '
         do  i = 1, i_N
            write(11,'(3e20.12)') r(i), T_r(i), dTdr_r(i)
         end do
         close(11)
         call flush(11)
      end if
   end subroutine io_save_profile

!--------------------------------------------------------------------------
!> Write out radial kinetic and magnetic energy profile
!!
!! Writes out file with 3 letter prefix 'label', needs
!! to be given the toroidal and poloidal `coll` variables
!! for the given field, and the correct conversion factor
!! Efac.
!-------------------------------------------------------------------------- 
   subroutine io_rad_profile(label, Tor, Pol, Efac, mes_oc)
      use meshs, only : rdom
      use parameters, only : d_pi
      use dyn_mpi
      use variables, only : var_coll_TorPol2qst

      type(coll), intent(in) :: Tor !< toroidal component of field
      type(coll), intent(in) :: Pol !< poloidal component of field
      character(3), intent(in) :: label !< id for field: vel/mag
      double precision, intent(in) :: Efac !< prefactor for energy calc
      type(rdom), intent(in) :: mes_oc !< outer core domain

      double precision :: E(cq%D%N), Epol(cq%D%N), Etor(cq%D%N), w, x  
      double precision :: eq(cq%D%N), es(cq%D%N), et(cq%D%N), Etot
      character(4) :: cnum
      integer :: i, N_
      _loop_lm_vars
      
      N_ = cq%D%N
      
      write(cnum,'(i4)') io_save3
      if(io_save3<1000) cnum(1:1) = '0'
      if(io_save3<100)  cnum(2:2) = '0'
      if(io_save3<10)   cnum(3:3) = '0'
      
      call var_coll_TorPol2qst( Tor, Pol, cq, cs, ct)
      Epol=0d0
      Etor=0d0
      E = 0d0
      _loop_lm_begin(cq)
        w = 4d0*d_Pi / dble(2*l+1)
        if (m_/=0) w = 4d0*w
        eq= w*( cq%Re(1:N_,nh)*cq%Re(1:N_,nh) + cq%Im(1:N_,nh)*cq%Im(1:N_,nh))
        es= w*( cs%Re(1:N_,nh)*cs%Re(1:N_,nh) + cs%Im(1:N_,nh)*cs%Im(1:N_,nh))
        et= w*( ct%Re(1:N_,nh)*ct%Re(1:N_,nh) + ct%Im(1:N_,nh)*ct%Im(1:N_,nh))
        Epol=Epol+eq+es
        Etor=Etor + et
      _loop_lm_end
#ifdef _MPI
      do i=1, N_
        call mpi_allreduce ( Epol(i), x, 1, mpi_double_precision, mpi_sum, mpi_comm_world, mpi_er)
        Epol(i) = x
        call mpi_allreduce ( Etor(i), x, 1, mpi_double_precision, mpi_sum, mpi_comm_world, mpi_er)
        Etor(i) = x
      end do
#endif
      Etor = Efac * Etor
      Epol = Efac * Epol
      E = Etor + Epol
      Etot = dot_product(E,cq%D%intr2dr(1:N_))
      if (mpi_rnk==0) then
        open(12, status='unknown',file=label//'_profile_'//cnum//'.dat')
        write(12,'(5A10,1e20.12)') '# r     ','    E    ',' Epol  ','  Etor   ',' #total_energy:',Etot
        do i=1, N_
          write(12,'(4e20.12)') mes_oc%r(i,1), E(i), Epol(i), Etor(i)
        end do
        close(12)
        call flush(12)
      end if
   end subroutine io_rad_profile

!**************************************************************************
 end module io
!**************************************************************************

