#include "../parallel.h"
!**************************************************************************
!> Non-linear terms of the several quantities and their evolution.
!!*************************************************************************
 module nonlinear
!*************************************************************************
   use variables, only : coll, spec, phys
   implicit none
   save

   type (coll), private :: cq,cs,ct
   type (spec), private :: sq,ss,st
   type (phys), private :: r,th,ph

 contains


!------------------------------------------------------------------------
!>  Initialise nonlinear stuff
!------------------------------------------------------------------------
   subroutine non_precompute()
   end subroutine non_precompute


!-------------------------------------------------------------------------
!>  nonlinear terms for velocity.  q,s->NPol, t->NTor
!------------------------------------------------------------------------
   subroutine non_velocity(mes_oc,vel_r, vel_t,vel_p,vel_curlr, vel_curlt, vel_curlp,&
                    mag_r, mag_t, mag_p, mag_curlr, mag_curlt, mag_curlp,&
                    cod_C, comp_C,&
                    vel_NTor, vel_npol)
      use meshs, only : rdom
      use parameters, only : i_pTh, i_Ph, b_boussinesq, d_pm, d_E, d_Ra, d_Pr, d_Sc, d_Ra_comp
      use variables, only : var_H, var_Th, thdom, harm, var_spec_init, var_spec2coll,&
                           & var_coll_qstllcurlr, var_coll_qstllcurlcurlr
      use legendre, only : leg_sin, leg_gp
      use transform, only : tra_rtp2qst
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(phys), intent(in) :: vel_r !< radial component of velocity
      type(phys), intent(in) :: vel_t !< theta component of velocity
      type(phys), intent(in) :: vel_p !< phi component of velocity
      type(phys), intent(in) :: vel_curlr !< radial component of curled velocity
      type(phys), intent(in) :: vel_curlt !< theta component of curled velocity
      type(phys), intent(in) :: vel_curlp !< phi component of curled velocity
      type(phys), intent(in) :: mag_r !< radial component of magnetic field
      type(phys), intent(in) :: mag_t !< theta component of magnetic field
      type(phys), intent(in) :: mag_p !< phi component of magnetic field
      type(phys), intent(in) :: mag_curlr !< radial component of curled magnetic field
      type(phys), intent(in) :: mag_curlt !< theta component of curled magnetic field
      type(phys), intent(in) :: mag_curlp !< phi component of curled magnetic field
      type(coll), intent(in) :: cod_C !< codensity
      type(coll), intent(in) :: comp_C !< composition 
      type(coll), intent(inout) :: vel_NTor !< nonlinear toroidal velocity
      type(coll), intent(inout) :: vel_npol !< nonlinear poloidal velocity
      double precision :: d1, d2, xi0_,xi1_
      double precision :: rm1_,Cmag,d_Cor,c1,Dz0

      double precision :: lsin(i_pTh), lgp(i_pTh)
      integer :: j,n,m,nh

      call var_spec_init(mes_oc,var_H, sq)
      call var_spec_init(mes_oc,var_H, ss)
      call var_spec_init(mes_oc,var_H, st)
      lsin = leg_sin(var_Th%pThi:var_Th%pThi+i_pTh-1)
      lgp  = leg_gp (var_Th%pThi:var_Th%pThi+i_pTh-1)
      
         			! advection   u x curlu
      j = var_Th%pTh
      do n = 1, mes_oc%pN
         do m = 0, i_Ph-1
            r%Re (:j,m,n) = (vel_t%Re(:j,m,n)*vel_curlp%Re(:j,m,n)  &
                                - vel_p%Re(:j,m,n)*vel_curlt%Re(:j,m,n) )
            th%Re(:j,m,n) = (vel_p%Re(:j,m,n)*vel_curlr%Re(:j,m,n)  &
                                - vel_r%Re(:j,m,n)*vel_curlp%Re(:j,m,n) )
            ph%Re(:j,m,n) = (vel_r%Re(:j,m,n)*vel_curlt%Re(:j,m,n)  &
                                - vel_t%Re(:j,m,n)*vel_curlr%Re(:j,m,n) )
         end do
      end do


! Compressible part of viscosity added in here
!  Fhat = (-2 xi u_r / r - d xi /dr u_r -1/3 xi^2 u_r) rhat
!         - {xi zeta_phi + u_theta (xi^2 +d xi /dr +2 xi/r)} thetahat
!         + {xi zeta_theta - (xi^2 +d xi /dr +2 xi/r)} phihat
!   Here zeta = curl u
!
!  ## Piotr ## 
!  due to switching back to the magnetic time scale this term gets 
!  multiplied by d_Pm (see equations)

      if (.not.b_boussinesq) then
        j = var_Th%pTh
        do n = 1, mes_oc%pN
           xi0_= mes_oc%rho(n+mes_oc%pNi-1,1)
           xi1_= mes_oc%rho(n+mes_oc%pNi-1,2)
           rm1_= mes_oc%r(n+mes_oc%pNi-1,-1)
           do m = 0, i_Ph-1
             r%Re(:j,m,n)  =  r%Re (:j,m,n) + d_Pm*vel_r%Re(:j,m,n) *      &
             (- 2.0d0*xi0_*rm1_- xi1_  - (xi0_**2)/3.0d0)
             d1 = 2.0d0*xi0_ * rm1_ + xi1_+ xi0_**2
             th%Re(:j,m,n) = th%Re(:j,m,n) - d_Pm*(xi0_*vel_curlp%Re(:j,m,n) &
             + d1*vel_t%Re(:j,m,n))
             ph%Re(:j,m,n) = ph%Re(:j,m,n) + d_Pm*(xi0_*vel_curlt%Re(:j,m,n) &
             - d1*vel_p%Re(:j,m,n))
           end do
        end do
      end if

!  ## Piotr ##           
!  the Coriolis force has now 2 Pm/E factor

         			! coriolis  - z x u,  z=(cos,-sin,0)

      do n = 1, mes_oc%pN
         do m = 0, i_Ph-1
            d_Cor=2.0d0*d_Pm/d_E
            r%Re (:j,m,n) = r%Re (:j,m,n) + d_Cor*lsin*vel_p%Re(:j,m,n)
            th%Re(:j,m,n) = th%Re(:j,m,n) + d_Cor*lgp *vel_p%Re(:j,m,n)
            ph%Re(:j,m,n) = ph%Re(:j,m,n) - d_Cor*(lgp*vel_t%Re(:j,m,n)  &
                                          + lsin*vel_r%Re(:j,m,n))
         end do
      end do

!  ## Piotr ##
! #CAJ Mar 2012 #  Lorentz force has now Pm/E*rho(:,0) factor
        			! lorentz
      do n = 1, mes_oc%pN
         if (.not.b_boussinesq) then
           Cmag=d_Pm/(d_E*mes_oc%rho(n+mes_oc%pNi-1,0))
         else
           Cmag=2 * d_Pm/d_E
         end if
         do m = 0, i_Ph-1
            r%Re(:j,m,n)  = r%Re(:j,m,n)  &
               +  Cmag*(mag_curlt%Re(:j,m,n)*mag_p%Re(:j,m,n)  &
               -  mag_curlp%Re(:j,m,n)*mag_t%Re(:j,m,n)) 
            th%Re(:j,m,n) = th%Re(:j,m,n)  &
               +  Cmag*(mag_curlp%Re(:j,m,n)*mag_r%Re(:j,m,n)  &
               -  mag_curlr%Re(:j,m,n)*mag_p%Re(:j,m,n))
            ph%Re(:j,m,n) = ph%Re(:j,m,n)  &
               +  Cmag*(mag_curlr%Re(:j,m,n)*mag_t%Re(:j,m,n)  &
               -  mag_curlt%Re(:j,m,n)*mag_r%Re(:j,m,n))
         end do
      end do

      call tra_rtp2qst(r,th,ph, sq,ss,st)

      call var_spec2coll(sq,cq, s2=ss,c2=cs, s3=st,c3=ct)


         			! codensity

!     Now set buoyancy in terms of basic state temperature gradient
!           d2 = d_Ra_now*mes_oc%r(n,i_buoyancy_exponent)*(d_Pm**2)/d_Pr

      do nh = 0, cq%H%pH1
         do n = 1, mes_oc%N
         if (.not.b_boussinesq) then
            d2 = -d_Ra*mes_oc%Temp0(n,1)*(d_Pm**2)/(d_Pr)
         else 
            d2 = d_Ra*mes_oc%r(n,1)*(d_Pm**2)/d_Pr
         end if            
            cq%Re(n,nh) = cq%Re(n,nh) + d2*cod_C%Re(n,nh)
            cq%Im(n,nh) = cq%Im(n,nh) + d2*cod_C%Im(n,nh)
         end do
      end do

               			! composition

!     Now set buoyancy in terms of basic state composition gradient
!           d2 = d_Ra_now*mes_oc%r(n,i_buoyancy_exponent)*(d_Pm**2)/d_Sc !!!!CHECK THIS!!!!!

      do nh = 0, cq%H%pH1
         do n = 1, mes_oc%N
         if (.not.b_boussinesq) then
            d2 = -d_Ra_comp*mes_oc%Temp0(n,1)*(d_Pm**2)/(d_Sc) ! Check this
         else
            d2 = d_Ra_comp*mes_oc%r(n,1)*(d_Pm**2)/d_Sc ! Check this
         end if
            cq%Re(n,nh) = cq%Re(n,nh) + d2*comp_C%Re(n,nh)
            cq%Im(n,nh) = cq%Im(n,nh) + d2*comp_C%Im(n,nh)
         end do
      end do

                                ! r/(l(l+1)) \hat{r} . curl N
      call var_coll_qstllcurlr(ct, vel_NTor)
                                ! r/(l(l+1)) \hat{r} . curl curl N
      call var_coll_qstllcurlcurlr(cq,cs, vel_NPol)

      if(mpi_rnk==0) then
         vel_NTor%Re(:,0) = 0d0
         vel_NTor%Im(:,0) = 0d0
         vel_NPol%Re(:,0) = 0d0
         vel_NPol%Im(:,0) = 0d0
      end if
      
   end subroutine non_velocity
   

!---------------------------------------------------------------------
!>  nonlinear terms for the codensity.  C->NPol
!------------------------------------------------------------------------
   subroutine non_codensity(mes_oc, vel_r, vel_t, vel_p, vel_curlr, vel_curlt, vel_curlp,&
                       vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2,&
                       mag_curlr,mag_curlt,mag_curlp,&
                       cod_gradr, cod_gradt, cod_gradp, cod_S, cod_NC)
      use meshs, only : rdom
      use parameters, only : d_Pr, d_Pm, d_Ra, d_E, b_boussinesq, i_ph, i_source_load, d_source
      use variables, only : var_H, var_Th, var_spec2coll, var_spec_init
      use transform, only : tra_phys2spec
      use dyn_mpi, only : mpi_rnk

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(phys), intent(in) :: vel_r !< radial component of velocity
      type(phys), intent(in) :: vel_t !< theta component of velocity
      type(phys), intent(in) :: vel_p !< phi component of velocity
      type(phys), intent(in) :: vel_curlr !< radial component of curled velocity
      type(phys), intent(in) :: vel_curlt !< theta component of curled velocity
      type(phys), intent(in) :: vel_curlp !< phi component of curled velocity
      type(phys), intent(in) :: vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2
      type(phys), intent(in) :: mag_curlr !< radial component of curled magnetic field
      type(phys), intent(in) :: mag_curlt !< theta component of curled magnetic field
      type(phys), intent(in) :: mag_curlp !< phi component of curled magnetic field
      type(phys), intent(in) :: cod_gradr !< radial component of grad(codensity)
      type(phys), intent(in) :: cod_gradt !< theta component of grad(codensity)
      type(phys), intent(in) :: cod_gradp !< phi component of grad(codensity)
      double precision, intent(in) :: cod_S(:) !< source term for codensity
      type(coll), intent(out) :: cod_NC !< nonlinear codensity term
      integer :: j, n, m
      double precision :: xi0_,Diss,rrr,d_qm1
      double precision :: dTemp0dr_,Temp0_,rho_ 
!     _loop_lm_vars
     
      d_qm1=d_Pr/d_Pm
      Diss=d_qm1**2/d_Ra       
      call var_spec_init(mes_oc,var_H, ss)
! - u . gradC -> -Pr u . grad C and - xi u_r term added in
      j = var_Th%pTh
      do n = 1, mes_oc%pN
!  Add in the dissipation terms if not Boussinesq
         if (.not.b_boussinesq) then
           xi0_= mes_oc%rho(n+mes_oc%pNi-1,1)
           Temp0_ = mes_oc%Temp0(n+mes_oc%pNi-1,0)
           dTemp0dr_= mes_oc%Temp0(n+mes_oc%pNi-1,1)
           rho_= mes_oc%rho(n+mes_oc%pNi-1,0)
           rrr= mes_oc%r(n+mes_oc%pNi-1,1)
         end if 
           do m = 0, i_Ph-1
             if (.not.b_boussinesq) then  
               r%Re(:j,m,n) =  &
               - d_qm1*vel_r%Re(:j,m,n) * cod_gradr%Re(:j,m,n)  &
               + (dTemp0dr_/Temp0_) * cod_gradr%Re(:j,m,n)  &
               - d_qm1*vel_t%Re(:j,m,n) * cod_gradt%Re(:j,m,n)  &
               - d_qm1*vel_p%Re(:j,m,n) * cod_gradp%Re(:j,m,n) &
               + (Diss/Temp0_)*(2.0d0*(vel_dur%Re(:j,m,n)**2) &
               + 2.0d0*((xi0_+1.0d0/rrr)*vel_r%Re(:j,m,n) + vel_dur%Re(:j,m,n) &
               + vel_qv1%Re(:j,m,n))**2 - 2.0d0*((xi0_*vel_r%Re(:j,m,n))**2)/3.d0 &
               + 2.0d0*((vel_qv1%Re(:j,m,n)+vel_r%Re(:j,m,n)/rrr)**2)  &
               + (2.0d0*vel_dutr%Re(:j,m,n) - vel_curlp%Re(:j,m,n))**2 &
               + (2.0d0*vel_dupr%Re(:j,m,n) + vel_curlt%Re(:j,m,n))**2 &
               + (2.0d0*vel_qv2%Re(:j,m,n) + vel_curlr%Re(:j,m,n))**2  &
!  Variable conductivity adds factor eta to ohmic dissipation May 3rd 2011.
               + mes_oc%eta(n+mes_oc%pNi-1,0)*((mag_curlr%Re(:j,m,n)**2+mag_curlp%Re(:j,m,n)**2 &
               + mag_curlt%Re(:j,m,n)**2)/(mes_oc%rho(n+mes_oc%pNi-1,0)*d_E)))
             else
!  Boussinesq terms with no dissipation
               r%Re(:j,m,n) =  &
               - d_qm1*vel_r%Re(:j,m,n) * cod_gradr%Re(:j,m,n)  &
               - d_qm1*vel_t%Re(:j,m,n) * cod_gradt%Re(:j,m,n)  &
               - d_qm1*vel_p%Re(:j,m,n) * cod_gradp%Re(:j,m,n)
             end if
           end do
      end do
      call tra_phys2spec(r, ss)
      call var_spec2coll(ss, cod_NC)
      ! internal source. This is taken from d_source in parameters.F90
      ! for simple cases (i_inhombc_C = 1)
      ! When i_inhombc_C = 2, 3 , it will use cod_S, which is loaded
      ! from either codBC.cdf.in or codS.txt.in 
      if(mpi_rnk==0) then  
         if( (i_source_load == 2) .or. (i_source_load == 3) ) then
            cod_NC%Re(:,0) = cod_NC%Re(:,0) + cod_S
         else if ( i_source_load ==1) then
            cod_NC%Re(:,0) = cod_NC%Re(:,0) + d_source
         end if
      end if         
   end subroutine non_codensity

!---------------------------------------------------------------------
!  nonlinear terms for the composition.  C->NPol
!------------------------------------------------------------------------
   subroutine non_composition(mes_oc, vel_r, vel_t, vel_p, vel_curlr, vel_curlt, vel_curlp,&
                      vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2,mag_curlr,mag_curlt,mag_curlp,&
                      comp_gradr, comp_gradt, comp_gradp, comp_S, comp_NC)

      use meshs, only : rdom
      use parameters, only : d_Sc, d_Pm, d_Ra_comp, b_boussinesq, i_Ph, i_source_comp_load, d_source_comp
      use variables, only : var_spec_init, var_H, var_Th, var_spec2coll
      use transform, only : tra_phys2spec
      use dyn_mpi, only : mpi_rnk
      
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(phys), intent(in) :: vel_r !< radial component of velocity
      type(phys), intent(in) :: vel_t !< theta component of velocity
      type(phys), intent(in) :: vel_p !< phi component of velocity
      type(phys), intent(in) :: vel_curlr !< radial component of curled velocity
      type(phys), intent(in) :: vel_curlt !< theta component of curled velocity
      type(phys), intent(in) :: vel_curlp !< phi component of curled velocity
      type(phys), intent(in) :: vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2
      type(phys), intent(in) :: mag_curlr !< radial component of curled magnetic field
      type(phys), intent(in) :: mag_curlt !< theta component of curled magnetic field
      type(phys), intent(in) :: mag_curlp !< phi component of curled magnetic field
      type(phys), intent(in) :: comp_gradr !< radial component of grad(composition)
      type(phys), intent(in) :: comp_gradt !< theta component of grad(composition)
      type(phys), intent(in) :: comp_gradp !< phi component of grad(composition)
      double precision, intent(in) :: comp_S(:) !< source term for composition
      type(coll), intent(out) :: comp_NC !< nonlinear composition term
      integer :: j, n, m
      double precision :: xi0_,Diss,rrr,d_qm1
      double precision :: dTemp0dr_,Temp0_,rho_
!     _loop_lm_vars

      d_qm1=d_Sc/d_Pm
      Diss=d_qm1**2/d_Ra_comp
      call var_spec_init(mes_oc,var_H, ss)
! - u . gradC -> -Pr u . grad C and - xi u_r term added in
      j = var_Th%pTh
      do n = 1, mes_oc%pN
!  Add in the dissipation terms if not Boussinesq
         if (.not.b_boussinesq) then
           xi0_= mes_oc%rho(n+mes_oc%pNi-1,1)
           Temp0_ = mes_oc%Temp0(n+mes_oc%pNi-1,0)
           dTemp0dr_= mes_oc%Temp0(n+mes_oc%pNi-1,1)
           rho_= mes_oc%rho(n+mes_oc%pNi-1,0)
           rrr= mes_oc%r(n+mes_oc%pNi-1,1)
         end if
           do m = 0, i_Ph-1
               !  Boussinesq terms with no dissipation
               r%Re(:j,m,n) =  &
               - d_qm1*vel_r%Re(:j,m,n) * comp_gradr%Re(:j,m,n)  &
               - d_qm1*vel_t%Re(:j,m,n) * comp_gradt%Re(:j,m,n)  &
               - d_qm1*vel_p%Re(:j,m,n) * comp_gradp%Re(:j,m,n)
           end do
      end do
      call tra_phys2spec(r, ss)
      call var_spec2coll(ss, comp_NC)
      ! internal source, fixed at d_source
      if(mpi_rnk==0) then
         if( (i_source_comp_load==2) .or. (i_source_comp_load==3) .or. (i_source_comp_load==4 )) then
            comp_NC%Re(:,0) = comp_NC%Re(:,0) + comp_S
         else if ( i_source_comp_load==1) then
            comp_NC%Re(:,0) = comp_NC%Re(:,0) + d_source_comp
         end if
      end if
   end subroutine non_composition
   

!-------------------------------------------------------------------------
!>  nonlinear terms for magnetic field.  q,s->NTor, t->NPol
!-------------------------------------------------------------------------
   subroutine non_magnetic(mes_oc, mes_ic, vel_r, vel_t, vel_p,&
                      mag_icTor, mag_icPol, mag_r, mag_t, mag_p,&
                      rot_omega, mag_bcRe, mag_bcIm,&
                      mag_NTor, mag_NPol, mag_icNTor, mag_icNPol)

      use meshs, only : rdom
      use parameters, only : i_H1, i_pTh, i_vel_bc, i_Ph
      use variables, only : var_spec_init, var_H, var_th, var_spec2coll, &
                          & var_coll_qstllcurlr, var_coll_qstllcurlcurlr, &
                          & var_coll_copy, var_coll_dph
      use legendre, only : leg_sin
      use dyn_mpi, only : mpi_rnk, r_rnk
      use transform, only : tra_rtp2qst

      type(rdom), intent(inout) :: mes_oc !< outer core domain
      type(rdom), intent(inout) :: mes_ic !< inner core domain
      type(phys), intent(in) :: vel_r !< radial component of velocity
      type(phys), intent(in) :: vel_t !< theta component of velocity
      type(phys), intent(in) :: vel_p !< phi component of velocity
      type(coll), intent(in) :: mag_icTor, mag_icPol
      type(phys), intent(in) :: mag_r !< radial component of magnetic field
      type(phys), intent(in) :: mag_t !< theta component of magnetic field
      type(phys), intent(in) :: mag_p !< phi component of magnetic field
      double precision, intent(in) :: rot_omega !< inner core rotation speed
      double precision, intent(inout) :: mag_bcRe(:,0:) !< real part mag BCs
      double precision, intent(inout) :: mag_bcIm(:,0:) !< imaginary part mag BCs
      type(coll), intent(inout) :: mag_NTor !< nonlinear part toroidal magnetic field
      type(coll), intent(inout) :: mag_NPol !< nonlinear part poloidal magnetic field
      type(coll), intent(inout) :: mag_icNTor !< nonlinear part toroidal magnetic field (inner core)
      type(coll), intent(inout) :: mag_icNPol !< nonlinear part poloidal magnetic field (inner core)
      
      integer :: nhm, j, n, m
      logical :: magbc
      double precision :: bcRe(0:i_H1), bcIm(0:i_H1), velph(i_pTh)

      call var_spec_init(mes_oc,var_H, sq)
      call var_spec_init(mes_oc,var_H, ss)
      call var_spec_init(mes_oc,var_H, st)

      magbc = (mes_ic%N/=1 .and. i_vel_bc>2)
         			! jump for stress-free cond inner core
      if(magbc) then
         r%Re (:,:,1) = 0d0
         th%Re(:,:,1) = mag_r%Re(:,:,1) * vel_t%Re(:,:,1)
         velph = rot_omega * mes_oc%r(1,1)  &
            * leg_sin(var_Th%pThi:var_Th%pThi+i_pTh-1)
         do m = 0, i_Ph-1 
            ph%Re(:,m,1) = mag_r%Re(:,m,1) * (vel_p%Re(:,m,1) - velph)
         end do
         n = mes_oc%N
         mes_oc%N = 1
         call tra_rtp2qst(r,th,ph, sq,ss,st)
         mes_oc%N = n
         bcRe = st%Re(:,1)
         bcIm = st%Im(:,1)
      end if
                                ! E = u x B
      j = var_Th%pTh
      do n = 1, mes_oc%pN
         do m = 0, i_Ph-1
            r%Re (:j,m,n) = vel_t%Re(:j,m,n) * mag_p%Re(:j,m,n)  &
                          - vel_p%Re(:j,m,n) * mag_t%Re(:j,m,n)
            th%Re(:j,m,n) = vel_p%Re(:j,m,n) * mag_r%Re(:j,m,n)  &
                          - vel_r%Re(:j,m,n) * mag_p%Re(:j,m,n)
            ph%Re(:j,m,n) = vel_r%Re(:j,m,n) * mag_t%Re(:j,m,n)  &
                          - vel_t%Re(:j,m,n) * mag_r%Re(:j,m,n)
         end do
      end do
      call tra_rtp2qst(r,th,ph, sq,ss,st)
    
      if(magbc .and. r_rnk==0) then  !NB no dr t in curl, curlcurl E
         st%Re(:,1) = bcRe
         st%Im(:,1) = bcIm
      end if
            
      call var_spec2coll(sq,cq, s2=ss,c2=cs, s3=st,c3=ct)
      
                                ! r/(l(l+1) \hat{r} . curl E
      call var_coll_qstllcurlr(ct, mag_NPol)
                                ! r/(l(l+1) \hat{r} . curl curl E
      call var_coll_qstllcurlcurlr(cq,cs, mag_NTor)
      
      if(magbc) then
         mag_bcRe(1,:) = - mag_NPol%Re(1,:)
         mag_bcIm(1,:) = - mag_NPol%Im(1,:)
      end if

               			! inner core
      call var_coll_copy(mag_icTor, ct)
      call var_coll_copy(mag_icPol, cq)
      n = mes_ic%N
      nhm = var_H%pH1
      call var_coll_dph(ct, mag_icNTor)
      call var_coll_dph(cq, mag_icNPol)
      mag_icNTor%Re(1:n,0:nhm) = - rot_omega * mag_icNTor%Re(1:n,0:nhm)
      mag_icNTor%Im(1:n,0:nhm) = - rot_omega * mag_icNTor%Im(1:n,0:nhm)
      mag_icNPol%Re(1:n,0:nhm) = - rot_omega * mag_icNPol%Re(1:n,0:nhm)
      mag_icNPol%Im(1:n,0:nhm) = - rot_omega * mag_icNPol%Im(1:n,0:nhm)

      if(mpi_rnk==0) then
         mag_NTor%Re(:,0) = 0d0
         mag_NTor%Im(:,0) = 0d0
         mag_NPol%Re(:,0) = 0d0
         mag_NPol%Im(:,0) = 0d0
         mag_icNTor%Re(:,0) = 0d0
         mag_icNTor%Im(:,0) = 0d0
         mag_icNPol%Re(:,0) = 0d0
         mag_icNPol%Im(:,0) = 0d0
      end if

   end subroutine non_magnetic


!-------------------------------------------------------------------------
!>  inner core rotation, torques
!------------------------------------------------------------------------
   subroutine non_rotation(mes_oc, mes_ic, vel_p, vel_curlt,mag_r, mag_p, rot_velTorq, rot_magTorq)

      use meshs, only : rdom
      use parameters, only : i_pth, i_vel_bc, i_Ph, d_Pm, d_E
      use dyn_mpi
      use legendre, only : leg_sin
      use variables, only : var_th
      use transform, only : tra_surfintegral

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(rdom), intent(in) :: mes_ic !< inner core domain
      type(phys), intent(in) :: vel_p !< phi component velocity
      type(phys), intent(in) :: vel_curlt !< theta component curl of velocity
      type(phys), intent(in) :: mag_r !< radial component of magnetic field
      type(phys), intent(in) :: mag_p !< phi component of magnetic field
      double precision, intent(out) :: rot_velTorq, rot_magTorq
      
      double precision :: S, vfac, mfac, lsin(i_pTh), d(2)
      integer :: m, n
      
      if(r_rnk==0) then
         lsin = leg_sin(var_Th%pThi:var_Th%pThi+i_pTh-1)
 
         vfac = mes_oc%r(1,3)*d_Pm
         mfac = mes_oc%r(1,3)*2*d_Pm/d_E

         			! viscous torque
         if(i_vel_bc>2) then
            rot_velTorq = 0d0
         else
            S = - 2d0/mes_oc%r(1,1)
            do m = 0, i_Ph-1
               r%Re(:,m,1) =  &
                  (S*vel_p%Re(:,m,1)-vel_curlt%Re(:,m,1)) * lsin
            end do
            call tra_surfintegral(r%Re(1,0,1), S)
            rot_velTorq = S * vfac
         end if
         			! lorentz torque
         if(mes_ic%N==1) then
            rot_magTorq = 0d0
         else
            do m = 0, i_Ph-1
               r%Re(:,m,1) = mag_r%Re(:,m,1) * mag_p%Re(:,m,1) * lsin
            end do
            call tra_surfintegral(r%Re(1,0,1), S)
            rot_magTorq = S * mfac
         end if
      end if

#ifdef _MPI
      if(mpi_rnk==0) d(1) = rot_velTorq
      if(mpi_rnk==0) d(2) = rot_magTorq      
      call mpi_bcast(d, 2, mpi_double_precision,  &
         0, mpi_comm_world, mpi_er)
      rot_velTorq = d(1)
      rot_magTorq = d(2)
#endif      
      
   end subroutine non_rotation
   

!------------------------------------------------------------------------
!>  get cfl max dt
!!
!! Saved to tim_cfl_dt from timestep.F90
!------------------------------------------------------------------------
   subroutine non_maxtstep(mes_oc, vel_r, vel_t, vel_p, mag_r,mag_t,mag_p)

      use meshs, only : rdom
      use parameters, only : i_pN, i_pTh, i_Ph, d_E, d_Pm, i_L, i_L1
      use variables, only : var_Th
      use timestep, only : tim_cfl_dt
      use dyn_mpi

      type(rdom), intent(in) :: mes_oc
      type(phys), intent(in) :: vel_r !< radial component of velocity
      type(phys), intent(in) :: vel_t !< theta component of velocity
      type(phys), intent(in) :: vel_p !< phi component of velocity
      type(phys), intent(in) :: mag_r !< radial component of magnetic field
      type(phys), intent(in) :: mag_t !< theta component of magnetic field
      type(phys), intent(in) :: mag_p !< phi component of magnetic field

      double precision :: dr, mx, r_(0:i_pN+1)
      double precision :: d, p(i_pTh, 0:i_Ph-1)
      integer :: N_,n, n__, NTh

      N_ = mes_oc%pN
      do n = 0, N_+1
         n__ = n+mes_oc%pNi-1
         if(n__>=1 .and. n__<=mes_oc%N)  r_(n) = mes_oc%r(n__,1)
      end do
      
      Nth = var_Th%pTh
      tim_cfl_dt = 1d99
         				! local
      do n = 1, N_

         if(n==1 .and. r_rnk==0) then
            dr = r_(2) - r_(1)
         else if(n==N_ .and. r_rnk==_Nr-1) then
            dr = r_(N_) - r_(N_-1)
         else
            dr = min( r_(n)-r_(n-1), r_(n+1)-r_(n) )
         end if
!  d_Ro -> d_E/d_Pm

         d = (d_E+d_E/d_Pm)/(2d0*dr)
         d = d*d
         p = mag_r%Re(:,:,n) * mag_r%Re(:,:,n)
         mx = maxval( p/dsqrt(p*(d_E/d_Pm)+d) + dabs(vel_r%Re(:,:,n)) )
         if(mx/=0d0) tim_cfl_dt = min( tim_cfl_dt, dr/mx )

         dr = r_(n)/dsqrt(dble(i_L*i_L1))
         d = (d_E+d_E/d_Pm)/(2d0*dr)
         d = (d_E)/(2d0*dr)
         d = d*d

         p =   mag_t%Re(:,:,n)*mag_t%Re(:,:,n)  &
            +  mag_p%Re(:,:,n)*mag_p%Re(:,:,n)
         mx = maxval( p/dsqrt(p*d_E/d_Pm+d)  &
            +  dsqrt( vel_t%Re(:,:,n)*vel_t%Re(:,:,n)  &
            +         vel_p%Re(:,:,n)*vel_p%Re(:,:,n) ) )
         if(mx/=0d0) tim_cfl_dt = min( tim_cfl_dt, dr/mx )
      end do
                                        ! global
      tim_cfl_dt = min(tim_cfl_dt, d_E/d_Pm)
      tim_cfl_dt = min(tim_cfl_dt, dsqrt(d_E))

#ifdef _MPI
      call mpi_allreduce(tim_cfl_dt, mx, 1, mpi_double_precision,  &
         mpi_min, mpi_comm_world, mpi_er)
      tim_cfl_dt = mx
#endif
  

   end subroutine non_maxtstep


!*************************************************************************
 end module nonlinear
!*************************************************************************

