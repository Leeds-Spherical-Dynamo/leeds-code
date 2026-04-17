!*************************************************************************
! main.f90 (executable)
!*************************************************************************
!
!
!*************************************************************************
#include "../parallel.h"
!> Driving routine for simulation.
!!
!! Runs the whole simulation. Calls initialise to set up file structure
!! and put RUNNING file in directory. Performs main timestepping loop,
!! then terminates the simulation using call to terminate, before
!! calling CLEANUP.
 PROGRAM MAIN
!*************************************************************************
   use meshs, only : rdom
   use variables, only : coll, phys
   use parameters, only : i_pH1, i_N, b_vel_tstep, b_cod_tstep, b_comp_tstep,&
                & b_mag_tstep, b_rot_tstep, d_timestep
   use timestep, only : tim_it, tim_new_dt, tim_t, tim_dt, tim_step, tim_check_cgce,&
                & tim_new_tstep
   use codensity, only : cod_transform, cod_matrices, cod_predictor, cod_corrector
   use composition, only : comp_transform, comp_matrices, comp_predictor, comp_corrector
   use velocity, only : vel_angmomclr, vel_transform, vel_matrices, vel_predictor, vel_corrector
   use magnetic, only : mag_transform, mag_matrices, mag_predictor, mag_corrector
   use nonlinear, only : non_rotation, non_codensity, non_composition, non_magnetic,&
                 & non_velocity, non_maxtstep
   use rotation, only : rot_predictor, rot_corrector
   use io, only : io_write2files

   implicit none
   double precision :: d_start !< wallclock start time of timestepping loop
   double precision :: d_stop !< used to check if wallclock stop time reached
!  _loop_lm_vars

   type(rdom) :: mes_oc !< outer core doamin
   type(rdom) :: mes_ic !< inner core domain

   type(coll) :: vel_tor !< toroidal component of velocity
   type(coll) :: vel_pol !< poloidal component of velocity
   type(coll) :: vel_ntor !< nonlinear term, toroidal component of velocity
   type(coll) :: vel_npol !< nonlinear term, poloidal component of velocity
   type(phys) :: vel_r !< radial componenet of velocity
   type(phys) :: vel_t !< theta compnent of velocity
   type(phys) :: vel_p !< phi component of velocity
   type(phys) :: vel_curlr !< radial component of curl(velocity)
   type(phys) :: vel_curlt !< theta component of curl(velocity)
   type(phys) :: vel_curlp !< phi component of curl(velocity)
   type(phys) :: vel_dur
   type(phys) :: vel_qv1
   type(phys) :: vel_dutr
   type(phys) :: vel_dupr
   type(phys) :: vel_qv2

   type(coll) :: mag_tor !< toroidal component of magnetic field
   type(coll) :: mag_pol !< poloidal component of magnetic field
   type(coll) :: mag_ntor !< nonlinear term, toroidal component of magnetic field
   type(coll) :: mag_npol !< nonlinear term, poloidal component of magnetic field
   type(coll) :: mag_icTor !< toroidal component of magnetic field, inner core
   type(coll) :: mag_icPol !< poloidal component of magnetic field, inner core
   type(coll) :: mag_icNTor !< nonlinear term, toroidal component of magnetic field, inner core
   type(coll) :: mag_icNPol !< nonlinear term, poloidal component of magnetic field, inner core
   type(phys) :: mag_r !< radial component of magnetic field
   type(phys) :: mag_t !< theta component of magnetic field
   type(phys) :: mag_p !< phi component of magnetic field
   type(phys) :: mag_curlr !< radial component of curl(magnetic field)
   type(phys) :: mag_curlt !< theta component of curl(magnetic field)
   type(phys) :: mag_curlp !< phi component of curl(magnetic field)
   double precision :: mag_bcRe(2,0:i_pH1) !< mag field toroidal Boundary Conditions, real part
   double precision :: mag_bcIm(2,0:i_pH1) !< mag field toroidal Boundary Conditions, imaginary part
   double precision :: mag_bc_PolRe(2,0:i_pH1) !< mag field poloidal Boundary Conditions, real part
   double precision :: mag_bc_PolIm(2,0:i_pH1) !< mag field toroidal Boundary Conditions, imaginary part

   type(coll) :: cod_C !< codensity (temperature in two-component convection)
   type(coll) :: cod_NC !< nonlinear term, codensity
   type(phys) :: cod_gradr !< radial component, grad(codensity)
   type(phys) :: cod_gradt !< theta component, grad(codensity)
   type(phys) :: cod_gradp !< phi component, grad(codensity)
   double precision :: cod_S(1:i_N) !< codensity source term, S(r)
   double precision :: cod_bcRe(2,0:i_pH1) !< codensity boundary conditions, real part
   double precision :: cod_bcIm(2,0:i_pH1) !< codensity boundary conditions, imaginary part
 
   type(coll) :: comp_C !< composition field
   type(coll) :: comp_NC !< nonlinear term, composition field
   type(phys) :: comp_gradr !< radial component, grad(composition)
   type(phys) :: comp_gradt !< theta component, grad(composition)
   type(phys) :: comp_gradp !< phi component, grad(composition)
   double precision :: comp_S(1:i_N) !< composition source term, S(r)
   double precision :: comp_bcRe(2,0:i_pH1) !< composition boundary conditions, real part
   double precision :: comp_bcIm(2,0:i_ph1) !< composition boundary conditions, imaginary part

   double precision :: rot_omega !< inner core rotation rate
   double precision :: rot_velTorq !< viscous torque on inner core
   double precision :: rot_magTorq !< magnetic torque on inner core
   double precision :: rot_inertia !< inertia of inner core

   call initialise()
   
      							! main loop
   do while(tim_it/=0 .or. .not.terminate())
      if(b_vel_tstep) call vel_angmomclr(mes_oc,vel_tor)
      if(b_vel_tstep) call vel_transform(mes_oc,vel_tor, vel_pol, vel_r, vel_t, vel_p,&
                                 vel_curlr, vel_curlt, vel_curlp,&
                                 vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2)
      if(b_cod_tstep) call cod_transform(mes_oc, cod_C, cod_gradr, cod_gradt, cod_gradp)
      if(b_comp_tstep) call comp_transform(mes_oc, comp_C, comp_gradr, comp_gradt, comp_gradp)
      if(b_mag_tstep) call mag_transform(mag_tor, mag_pol, mag_r, mag_t, mag_p, mag_curlr, mag_curlt, mag_curlp)
      if(b_rot_tstep) call non_rotation(mes_oc, mes_ic, vel_p, vel_curlt, mag_r, mag_p, rot_velTorq, rot_magTorq)
      if(b_vel_tstep) call non_velocity(mes_oc, vel_r, vel_t, vel_p,&
                                 vel_curlr, vel_curlt, vel_curlp, mag_r, mag_t, mag_p, mag_curlr, mag_curlt,&
                                  mag_curlp, cod_C, comp_C, vel_Ntor, vel_npol)
      if(b_cod_tstep) call non_codensity(mes_oc, vel_r,vel_t, vel_p,&
                                 vel_curlr, vel_curlt, vel_curlp,&
                                 vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2,&
                                 mag_curlr, mag_curlt, mag_curlp,cod_gradr,cod_gradt,cod_gradp, cod_S,cod_NC)
      if(b_comp_tstep) call non_composition(mes_oc, vel_r, vel_t, vel_p,&
                                 vel_curlr, vel_curlt, vel_curlp,&
                                 vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2,&
                                 mag_curlr, mag_curlt, mag_curlp,&
                                 comp_gradr, comp_gradt, comp_gradp, comp_S, comp_NC)
      if(b_mag_tstep) call non_magnetic(mes_oc, mes_ic, vel_r, vel_t, vel_p, mag_icTor, mag_icPol,&
                               mag_r, mag_t,mag_p, rot_omega, mag_bcRe, mag_bcIm,&
                               mag_NTor, mag_NPol, mag_icNTor, mag_icNPol)

      if(tim_it==0) then
         if(d_timestep<0d0) then
            call non_maxtstep(mes_oc, vel_r, vel_t, vel_p, mag_r, mag_t, mag_p)
            call tim_new_tstep()
            if(tim_new_dt) then
               call vel_matrices(mes_oc)
               call cod_matrices(mes_oc)
               call comp_matrices(mes_oc)
               call mag_matrices(mes_oc, mes_ic)
            endif  
         end if        
         if(b_rot_tstep) call rot_predictor(rot_inertia, rot_velTorq, rot_magTorq, rot_omega)
         if(b_vel_tstep) call vel_predictor(mes_oc, vel_ntor, vel_npol, rot_omega, vel_tor, vel_pol)

         if(b_cod_tstep) call cod_predictor(cod_Nc, cod_bcRe, cod_bcIm, cod_C)
         if(b_comp_tstep) call comp_predictor(comp_NC, comp_bcRe, comp_bcIm, comp_C)

         if(b_mag_tstep) call mag_predictor(mes_oc, mes_ic, mag_ntor, mag_npol, mag_icNTor, mag_icNpol,&
                                  mag_bcRe, mag_bcIm, mag_bc_PolRe, mag_bc_PolIm,&
                                  mag_tor, mag_pol, mag_icTor, mag_icPol)

         tim_it = 1
      else
         if(b_rot_tstep) call rot_corrector(rot_inertia, rot_velTorq, rot_magTorq, rot_omega)
         if(b_vel_tstep) call vel_corrector(mes_oc, vel_ntor,vel_npol,rot_omega, vel_tor, vel_pol)
         if(b_cod_tstep) call cod_corrector(cod_NC, cod_C)
         if(b_comp_tstep) call comp_corrector(comp_NC, comp_C)
         if(b_mag_tstep) call mag_corrector(mes_oc, mes_ic, mag_ntor, mag_npol, mag_icNTor, mag_icNPol,&
                                 mag_bcRe, mag_bcIm,&
                                 mag_tor, mag_pol, mag_icTor, mag_icPol)
         call tim_check_cgce()
      end if

      if(tim_it==0) then
         tim_t    = tim_t    + tim_dt
         tim_step = tim_step + 1
         call io_write2files(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol,mag_icTor, mag_icPol,&
                               rot_omega, rot_velTorq, rot_magTorq, cod_C, comp_C)
      end if
   end do
      							! end main loop   

   call cleanup()

 contains


!-------------------------------------------------------------------------
!  Termaination conditions
!-------------------------------------------------------------------------
   !>  Termination check
   !!   
   !! Sets terminate=.true. if i_maxstep is reached,
   !! if wall time exceeds d_cpuhours, or if 
   !! file "RUNNING" no longer exists
   logical function terminate()
      use dyn_mpi
      use parameters, only : i_maxtstep, i_save_rate2, d_cpuhours

      logical :: file_exist
            
      if(mpi_rnk==0) then
         terminate = .false.
      
         if(tim_step==i_maxtstep) then
            terminate = .true.
            print*, 'maxtstep reached!'
         end if

         call get_time(d_stop)
         if(dble(d_stop-d_start)/36d2 >= d_cpuhours) then
            terminate = .true.
            print*, 'cpuhours reached!'
         end if

         if( modulo(tim_step,i_save_rate2)==0) then
            inquire(file='RUNNING', exist=file_exist)
            if(.not. file_exist) then
               terminate = .true.
               print*, 'RUNNING deleted !'
            end if
         end if
      end if

#ifdef _MPI
      call mpi_bcast(terminate,1,mpi_logical, 0,mpi_comm_world,mpi_er)
#endif

   end function terminate


!-------------------------------------------------------------------------
!  Initialisation
!-------------------------------------------------------------------------
   !> Initialise Data for simulation
   subroutine initialise()

      use dyn_mpi
      use io, only : io_precompute, code_version, io_load_codbc, io_load_codStxt, io_load_compBC, &
                     & io_comp_flayer, io_openfiles, io_load_state, io_load_radial_txtfiles
      use parameters, only : par_precompute, i_inhombc_C, i_source_load, i_inhombc_Comp,&
                       & i_source_comp_load, b_inhom_magPol, b_fixed_at_inner, b_fixed_at_outer, i_Nic
      use meshs, only : mes_precompute
      use legendre, only : leg_precompute
      use variables, only : var_precompute
      use transform, only : tra_precompute
      use timestep, only : tim_precompute
      use rotation, only : rot_precompute
      use velocity, only : vel_precompute
      use codensity, only : cod_precompute
      use composition, only : comp_precompute
      use magnetic, only : mag_precompute
      use nonlinear, only : non_precompute

      logical :: file_exist

      call mpi_precompute()
      if(mpi_rnk==0) then
         call system('touch PRECOMPUTING')
         call system('echo $HOSTNAME > HOST')
         open(66, file='LOG_copy.dat')
         write(66,*) 'This file duplicates some of the output to the console.'
         write(66,*) 'Leeds-Spherical-Dynamo version: ',code_version
      end if

!     if(mpi_rnk==0)  print*, 'loading txt files...'
      if(mpi_rnk==0)  print*, 'Leeds-Spherical-Dynamo version: ',code_version
      if(mpi_rnk==0)  print*, 'precomputing function requisites...'
#ifdef _SHTNS
      if(mpi_rnk==0)  print*, 'using SHTNS for transforms'
#else
      if(mpi_rnk==0)  print*, 'using native code for transforms'
#endif
      call io_precompute(mes_ic) 
      call par_precompute()
      call io_load_radial_txtfiles(mes_oc, mes_ic)
      call mes_precompute(mes_oc, mes_ic)
      call leg_precompute()
      call var_precompute()
      call tra_precompute(1)
      call tim_precompute()
      call rot_precompute(mes_oc, rot_omega, rot_inertia, rot_velTorq, rot_magTorq)
      call vel_precompute(mes_oc, vel_Tor, vel_Pol)
      call cod_precompute(mes_oc, cod_C, cod_S, cod_bcRe, cod_bcIm)
      call comp_precompute(mes_oc, comp_C, comp_NC, comp_S, comp_bcRe, comp_bcIm)
      call mag_precompute(mes_oc, mes_ic, mag_Tor, mag_Pol, mag_icTor, mag_icPol,&
                  mag_bcRe, mag_bcIm, mag_bc_PolRe, mag_bc_PolIm)
      call non_precompute()
!     call io_precompute()
   
      ! load codensity boundary conditions
      if (i_inhombc_C==2) then
         call io_load_codBC(mes_oc, cod_C, cod_S, cod_bcRe, cod_bcIm)
         if (mpi_rnk==0) print*, ' loaded codBC cdf file'
      else if (i_inhombc_C==3) then
         call io_load_codBCtxt(cod_C, cod_bcRe, cod_bcIm)
         if (mpi_rnk==0) print*, ' loaded codBCtxt file'
      end if

      ! load codensity source from codBC.cdf.in
      if (i_source_load == 2) then
         if (i_inhombc_C .ne. 2) then
            stop 'specified loading source from codBC.cdf but didnt load codBC.cdf ' 
            return
         end if
      end if

      ! load codensity source from codS.txt.in
      if (i_source_load == 3) then
         if (i_inhombc_C == 2) then
            if (mpi_rnk==0) print*, ' warning, overwriting cod_S with codS.txt.in'
         end if
         call io_load_codStxt(mes_oc, cod_S)
      end if
         
      ! load composition boundary conditions
      if (i_inhombc_Comp==2) then
        call io_load_compBC(mes_oc, comp_C, comp_S, comp_bcRe, comp_bcIm)
        if (mpi_rnk==0) print*, ' loaded compBC cdf file'
      else if (i_inhombc_Comp==3) then
         call io_load_compBCtxt(comp_C, comp_bcRe, comp_bcIm)
         if (mpi_rnk==0) print*, ' loaded compBCtxt file'
      end if

      ! load composition source from compBC.cdf.in
      if (i_source_comp_load == 2) then
         if (i_inhombc_comp .ne. 2) then
            stop 'specified loading source from compBC.cdf but didnt load compBC.cdf'
            return
         end if
      end if
      
      ! load composition source from compS.txt.in
      if (i_source_comp_load == 3) then
         if (i_inhombc_comp == 2) then
            if (mpi_rnk==0) print*, ' warning, overwriting comp_S with compS.txt.in'
         end if
         call io_load_compStxt(mes_oc, comp_S)
      end if

      ! create flayer for composition
      if (i_inhombc_Comp==4) then
         if (i_source_comp_load .ne. 4) then
            stop 'flayer params not consistent'
         endif
         call io_comp_flayer(mes_oc, comp_C, comp_S, comp_bcRe, comp_bcim)
         if (mpi_rnk==0) print *, ' created flayer using composition '
      end if

!     If inhomogeneous magnetic poloidal condition imposed (magnetoconvection)
!
      if (b_inhom_magPol) then
        call io_load_magPolBC(mag_pol, mag_bc_PolRe, mag_bc_PolIm)
        if (mpi_rnk==0) print*, ' loaded magPolBC file'
!     end if

         if ((.not.b_fixed_at_inner).and.(.not.b_fixed_at_outer)) then
            if (mpi_rnk==0) print*,' One of b_fixed_at_inner and b_fixed_at_outer must be true'
            stop ' One of b_fixed_at_inner and b_fixed_at_outer must be true'
         end if
         if (b_fixed_at_inner.and.i_Nic/=1) then
            if (mpi_rnk==0) print*,' Must have i_Nic=1 if b_fixed_at_inner = .true.'
            stop ' Must have i_Nic=1 if b_fixed_at_inner = .true.'
         end if

      end if

      if(mpi_rnk==0)  print*,    'loading state...'
      call io_load_state(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol, mag_icTor, mag_icPol,&
                cod_C, comp_C, rot_omega)
      if(mpi_rnk==0)  print*,    'state loaded ...' 
      call vel_angmomclr(mes_oc, vel_tor)
      call vel_transform(mes_oc, vel_tor, vel_pol, vel_r, vel_t, vel_p,&
                  vel_curlr, vel_curlt, vel_curlp,&
                  vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2)
      if(mpi_rnk==0)  print*,    'vel_trans fin ..'    

      call cod_transform(mes_oc,cod_C,cod_gradr, cod_gradt, cod_gradp)
      if(mpi_rnk==0)  print*,    'cod_trans fin ..'  

      call comp_transform(mes_oc, comp_C, comp_gradr, comp_gradt, comp_gradp)
      if(mpi_rnk==0)  print*,    'comp_trans fin ..'

      call mag_transform(mag_tor, mag_pol, mag_r, mag_t, mag_p, mag_curlr, mag_curlt, mag_curlp)
      if(mpi_rnk==0)  print*,    'mag_trans fin ..'  


      call non_rotation(mes_oc, mes_ic, vel_p, vel_curlt, mag_r, mag_p, rot_velTorq, rot_magTorq)
      if(mpi_rnk==0)  print*,    'non_rotat fin ..'   
      call non_maxtstep(mes_oc, vel_r, vel_t, vel_p,mag_r, mag_t, mag_p)

      if(mpi_rnk==0)  print*, 'computing timestepping matrices...'
      call vel_matrices(mes_oc)
      call cod_matrices(mes_oc)
      call comp_matrices(mes_oc)
      call mag_matrices(mes_oc, mes_ic)

      if(mpi_rnk==0)  print*, 'initialising output files...'
      call io_openfiles()
      call io_write2files(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol, mag_icTor, mag_icPol,&
                      rot_omega, rot_velTorq, rot_magTorq, cod_C, comp_C)

      if(mpi_rnk==0) then
         open(99, file='PRECOMPUTING')
         close(99, status='delete')

         open(99, file='RUNNING')
         write(99,*) 'This file indicates a job running in this directory.'
         write(99,*) 'Delete this file to cleanly terminate the process.'
         close(99)
         print*, 'timestepping.....'
      end if
      
      call get_time(d_start)
   
   end subroutine initialise

!-------------------------------------------------------------------------
!  Program closure
!-------------------------------------------------------------------------
   !> Cleanup of modules, remove 'RUNNING' file, finalize mpi. 
   subroutine cleanup()
      use dyn_mpi
      use io, only : io_save_state, io_closefiles
      logical :: file_exist
   
      if(mpi_rnk==0) then
         print*, 'cleanup...'
         call get_time(d_stop)
         print*, ' sec/step = ', (d_stop-d_start)/real(tim_step)
         print*, ' cpu time = ', int((d_stop-d_start)/6d1), ' mins.'
      end if
      
      call io_save_state(mes_oc, mes_ic, vel_tor, vel_pol, mag_tor, mag_pol, mag_icTor, mag_icPol,&
                   cod_C, comp_C, rot_omega)
      call io_closefiles()

      if(mpi_rnk==0) then
         inquire(file='RUNNING', exist=file_exist)
         if(file_exist) open(99, file='RUNNING')
         if(file_exist) close(99, status='delete')
      end if      

#ifdef _MPI
      call mpi_barrier(mpi_comm_world, mpi_er)
      call mpi_finalize(mpi_er)
#endif
      if(mpi_rnk==0) print*, '...done!'

   end subroutine cleanup

        
!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************


