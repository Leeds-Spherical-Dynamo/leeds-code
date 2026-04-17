#include "../parallel.h"
!*************************************************************************
!> Provides the infrastructure to deal with the velocity field and its boundary conditions.
!!*************************************************************************
 module velocity
!*************************************************************************
   use variables, only : coll, spec
   use meshs, only : lumesh, mesh
   use parameters, only : i_L1, i_N
   
   implicit none
   save

   type (coll), private :: Tor !< temporary coll for computation
   type (coll), private :: Pol !< temporary coll for computation
   type (coll), private :: NTor !< temporary coll for computation
   type (coll), private :: NPol !< temporary coll for computation

   type (lumesh), private :: XTor(0:i_L1) !< XTor lumesh for computation
   type (lumesh), private :: XPol(0:i_L1) !< XPol lumesh for computation
   type (lumesh), private :: XGre(0:i_L1) !< XGre lumesh for computation
   type (mesh),   private :: YTor(0:i_L1) !< YTor mesh for computation
   type (mesh),   private :: YPol(0:i_L1) !< Ypol mesh for computation
   double precision, private :: Gre (i_N,2,i_L1)
   double precision, private :: invG(  2,2,i_L1)

   type (coll), private :: cq !< coll for computation (qst, q)
   type (coll), private :: cs !< coll for computation (qst, s)
   type (coll), private :: ct !< coll for computation (qst, t)
   type (coll), private :: cqr !< coll for computation (qst, qr)
   type (coll), private :: csr !< coll for computation (qst, sr)
   type (coll), private :: ctr !< coll for computation (qst, tr)
   type (spec), private :: sq !< spec for computation (qst, q)
   type (spec), private :: ss !< spec for computation (qst, s)
   type (spec), private :: st !< spec for computation (qst, t)
   type (spec), private :: sqr !< spec for computation (qst, qr)
   type (spec), private :: ssr !< spec for computation (qst, sr)
   type (spec), private :: str !< spec for computation (qst, tr)
   
 contains

!------------------------------------------------------------------------
!>  initialise velocity/pressure field
!------------------------------------------------------------------------
   subroutine vel_precompute(mes_oc, vel_Tor, vel_pol)
      use meshs, only : rdom
      use variables, only : var_H, var_coll_init

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(inout) :: vel_Tor !< toroidal part of velocity
      type(coll), intent(inout) :: vel_pol !< poloidal part of velocity

      call var_coll_init(mes_oc,var_H, vel_Tor)
      call var_coll_init(mes_oc,var_H, vel_Pol)
   end subroutine vel_precompute


!------------------------------------------------------------------------
!>  precomputation of velocity/pressure field matrices
!! 
!! Sets up XTor, YTor, XPol, YPol, XGre, invG
!------------------------------------------------------------------------
   subroutine vel_matrices(mes_oc)
      use meshs, only : rdom, mes_lu_find, mes_mat_invert
      use variables, only : harm, var_coll_init
      use parameters, only : d_Pm, b_boussinesq, smart_hd_vel, ldiff_vel, qh_vel,&
                           & d_implicit
      use timestep, only : tim_lumesh_vel_X, tim_lumesh_flayer_X, tim_mesh_vel_Y, tim_mesh_flayer_Y,&
                          & tim_invX, tim_zerobc, tim_lumesh_gre_X

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type (harm) :: H
      double precision :: c1, c2, c2_default, hyp_factor
      integer :: N_, l

      c1 = 1d0
      c2_default = d_Pm

      do l = 1, i_L1
         ! smart hyperdiffusion for toroidal velocity field
         if (b_boussinesq .and. smart_hd_vel) then
            if (l>=ldiff_vel) then
               hyp_factor = qh_vel**(l-ldiff_vel)
               c2 = hyp_factor * c2_default
            else
               c2 = c2_default
            end if
         else
            c2 = c2_default
         end if

         if (.not.b_boussinesq) then
           call tim_lumesh_vel_X      ( mes_oc, c1, c2, l, XTor(l) )
         else
          call tim_lumesh_flayer_X      ( mes_oc, c1, c2, l, XTor(l) )
         end if  
         call vel_bc_Tor        ( mes_oc, XTor(l), l)
         call mes_lu_find       ( mes_oc, XTor(l) )
         if (.not.b_boussinesq) then
           call tim_mesh_vel_Y        ( mes_oc, c1, c2, l, YTor(l) )
         else
           call tim_mesh_flayer_Y        ( mes_oc, c1, c2, l, YTor(l) )
         end if
      end do


      do l = 1, i_L1
         ! smart hyperdiffusion for poloidal velocity field
         if (b_boussinesq .and. smart_hd_vel) then
            if (l>=ldiff_vel) then
               hyp_factor = qh_vel**(l-ldiff_vel)
               c2 = hyp_factor * c2_default
            else
               c2 = c2_default
            end if
         else
            c2 = c2_default
         end if

         if (.not.b_boussinesq) then
           call tim_lumesh_vel_X      ( mes_oc, c1, c2, l, XPol(l) )
         else
           call tim_lumesh_flayer_X      ( mes_oc, c1, c2, l, XPol(l) )
         end if 
         call vel_bc_Pol        ( mes_oc, XPol(l), l, mes_oc)
         call mes_lu_find       ( mes_oc, XPol(l) )
         if (.not.b_boussinesq) then
           call tim_mesh_vel_Y        ( mes_oc, c1, c2, l, YPol(l) )
         else
           call tim_mesh_flayer_Y        ( mes_oc, c1, c2, l, YPol(l) )
         end if
      end do

                                ! greens functions & influence matrix
      c1 = 0d0
      c2 = 1d0/d_implicit
      N_ = mes_oc%N
      do l = 1, i_L1
         if (.not.b_boussinesq) then
           call tim_lumesh_Gre_X      ( mes_oc, c1, c2, l, XGre(l))
         else
           call tim_lumesh_flayer_X      ( mes_oc, c1, c2, l, XGre(l))
         end if
         call vel_bc_Gre        ( mes_oc, XGre(l), l)
         call mes_lu_find       ( mes_oc, XGre(l) )
      end do

      call var_coll_init(mes_oc,H, Pol)
      do l = 1, i_L1
         Pol%H%l(0) = l
         Pol%H%m(0) = 0
         Pol%H%pH0 = 0
         Pol%H%pH1 = 0
         Pol%Re( :,0) = 0d0
         Pol%Im( :,0) = 0d0
         Pol%Re( 1,0) = 1d0
         Pol%Im(N_,0) = 1d0
         call tim_invX(.false.,XGre, Pol)
         call tim_zerobc(Pol)
         call tim_invX(.false.,XPol, Pol)
         Gre(:,1,l) = Pol%Re(:,0)
         Gre(:,2,l) = Pol%Im(:,0)
         invG(1,1,l) = Gre( 1,1,l)
         invG(1,2,l) = Gre( 1,2,l)
         invG(2,1,l) = Gre(N_,1,l)
         invG(2,2,l) = Gre(N_,2,l)
         call mes_mat_invert(2,invG(1,1,l),2)
      end do

   end subroutine vel_matrices


!---------------------------------------------------------
!>  Velocity toroidal boundary condtions
!!
!--------------------------------------------------------
   subroutine vel_bc_Tor(D,A,l)
      use meshs, only : rdom, lumesh
      use parameters, only : i_vel_bc, i_Kl

      type (rdom),   intent(in)  :: D !< outer core domain
      type (lumesh), intent(out) :: A !< lumesh for l
      integer,       intent(in)  :: l !< spherical harmonic degree

      integer :: j

      if(i_vel_bc==1 .or. i_vel_bc==2) then
         A%M(2*i_KL+1,1)   = 1d0
      else
         do j = 1, 1+i_KL
            A%M(2*i_KL+1+1-j,j) = D%dr(1)%M(i_KL+1+1-j,j)
         end do
         A%M(2*i_KL+1,1) = A%M(2*i_KL+1,1) - D%r(1,-1)
      end if
      
      if(i_vel_bc==1 .or. i_vel_bc==3) then
         A%M(2*i_KL+1,D%N) = 1d0
      else
         do j = D%N-i_KL, D%N
            A%M(2*i_KL+1+D%N-j,j) = D%dr(1)%M(i_KL+1+D%N-j,j)
         end do
         A%M(2*i_KL+1,D%N) = A%M(2*i_KL+1,D%N) - D%r(D%N,-1)
      end if

   end subroutine vel_bc_Tor
   
!---------------------------------------------------------
!>  Velocity Poloidal boundary condtions
!!
!--------------------------------------------------------
   subroutine vel_bc_Pol(D,A,l, mes_oc)
      use meshs, only : rdom, lumesh
      use parameters, only : i_vel_bc, i_KL, b_boussinesq

      type(rdom), intent(in) :: mes_oc
      type (rdom),   intent(in)  :: D !< outer core domain
      type (lumesh), intent(out) :: A !< lu mesh
      integer,       intent(in)  :: l !< spherical harmonic degree

      integer :: j

      if(i_vel_bc==1 .or. i_vel_bc==2) then
         do j = 1, 1+i_KL
            A%M(2*i_KL+1+1-j,j) = D%dr(1)%M(i_KL+1+1-j,j)
         end do
      else
         if (b_boussinesq) then
            do j = 1, 1+i_KL
               A%M(2*i_KL+1+1-j,j) = D%dr(2)%M(i_KL+1+1-j,j)
            end do
         else ! (if anelastic, add compressible part)
            do j = 1, 1+i_KL
               A%M(2*i_KL+1+1-j,j) = D%dr(2)%M(i_KL+1+1-j,j) + &
    !#CAJ Compressible part of the bc
               mes_oc%rho(1,1)*D%dr(1)%M(i_KL+1+1-j,j)
            end do
         end if
      end if
      
      if(i_vel_bc==1 .or. i_vel_bc==3) then
         do j = D%N-i_KL, D%N
            A%M(2*i_KL+1+D%N-j,j) = D%dr(1)%M(i_KL+1+D%N-j,j)
         end do
      else
         if (.not.b_boussinesq) then
           do j = D%N-i_KL, D%N
              A%M(2*i_KL+1+D%N-j,j) = D%dr(2)%M(i_KL+1+D%N-j,j) + &
! #CAJ Compressible part of the bc
              mes_oc%rho(D%N,1)*D%dr(1)%M(i_KL+1+D%N-j,j)    
           end do
         else
           do j = D%N-i_KL, D%N
             A%M(2*i_KL+1+D%N-j,j) = D%dr(2)%M(i_KL+1+D%N-j,j)
           end do
         end if 
      end if

   end subroutine vel_bc_Pol


!---------------------------------------------------------
!> vel_bc_Gre
!!
!--------------------------------------------------------
   subroutine vel_bc_Gre(D,A,l)
      use meshs, only : rdom, lumesh
      use parameters, only : i_Kl
      type (rdom),   intent(in)  :: D !< outer core mesh
      type (lumesh), intent(out) :: A !< greens lu mesh
      integer,       intent(in)  :: l
      A%M(2*i_KL+1,1)   = 1d0
      A%M(2*i_KL+1,D%N) = 1d0
   end subroutine vel_bc_Gre


!------------------------------------------------------------------------
!>  set boundary condition rotating core on T10 mode.
!------------------------------------------------------------------------
   subroutine vel_setbc_Tor(mes_oc,vel_tor, rot_omega, a)
      use meshs, only : rdom
      use parameters, only : i_vel_bc
      use timestep, only : tim_it

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(in) :: vel_tor !< toroidal part of velocity
      double precision,  intent(in) :: rot_omega !< inner core rotation
      type (coll), intent(out) :: a !< field to apply rotating core to

      integer :: N_
      _loop_lm_vars

      N_ = a%D%N
      _loop_lm_begin(a)
         a%Re(1 ,nh) = 0d0
         a%Re(N_,nh) = 0d0
         a%Im(1 ,nh) = 0d0
         a%Im(N_,nh) = 0d0
         if(i_vel_bc>2) cycle
         if(l/=1 .or. m_/=0) cycle
         a%Re(1,nh) = rot_omega * mes_oc%r(1,nh)
         if(tim_it/=0) a%Re(1,nh) = a%Re(1,nh) - vel_Tor%Re(1,nh)
      _loop_lm_end

   end subroutine vel_setbc_Tor


!------------------------------------------------------------------------
!>  Evaluate in physical space  u  and  curl u
!------------------------------------------------------------------------
   subroutine vel_transform(mes_oc, vel_Tor, vel_pol,vel_r,vel_t,vel_p,&
                             vel_curlr,vel_curlt,vel_curlp,&
                             vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2)
      use meshs, only : rdom
      use variables, only : phys, var_coll_TorPol2qstderiv, var_coll2spec, var_coll_qstcurl,&
                           &var_coll_TorPol2qst_mag
      use transform, only : tra_qst2rtp_vel, tra_qst2rtp
      use parameters, only : b_boussinesq

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(in) :: vel_Tor !< toroidal part of velocity
      type(coll), intent(in) :: vel_pol !< poloidal part of velocity
      type(phys), intent(inout) :: vel_r !< radial component of velocity
      type(phys), intent(inout) :: vel_t !< theta component of velocity
      type(phys), intent(inout) :: vel_p !< phi component of velocity
      type(phys), intent(inout) :: vel_curlr !< radial component of curled field
      type(phys), intent(inout) :: vel_curlt !< theta component of curled field
      type(phys), intent(inout) :: vel_curlp !< phi component of curled field
      type(phys), intent(inout) :: vel_dur
      type(phys), intent(inout) :: vel_qv1
      type(phys), intent(inout) :: vel_dutr
      type(phys), intent(inout) :: vel_dupr
      type(phys), intent(inout) :: vel_qv2

      if (.not.b_boussinesq) then
        call var_coll_TorPol2qstderiv(vel_Tor,vel_Pol, cq,cs,ct,cqr,csr,ctr)
        call var_coll2spec(cq,sq, c2=cs,s2=ss, c3=ct,s3=st)
        call var_coll2spec(cqr,sqr, c2=csr,s2=ssr, c3=ctr,s3=str)
        call tra_qst2rtp_vel(mes_oc,sq,ss,st,sqr,ssr,str,vel_r,vel_t,vel_p, &
          vel_dur,vel_qv1,vel_dutr,vel_dupr,vel_qv2)
      else
        call var_coll_TorPol2qst_mag(vel_Tor,vel_Pol, cq,cs,ct)
        call var_coll2spec(cq,sq, c2=cs,s2=ss, c3=ct,s3=st)
        call tra_qst2rtp(sq,ss,st, vel_r,vel_t,vel_p)
      end if

      call var_coll_qstcurl(cq,cs,ct, cq,cs,ct)
      call var_coll2spec(cq,sq, c2=cs,s2=ss, c3=ct,s3=st)
      call tra_qst2rtp(sq,ss,st, vel_curlr,vel_curlt,vel_curlp)

   end subroutine vel_transform


!------------------------------------------------------------------------
!> Velocity Predictor
!!
!! Tor:
!! - N  := vN,            save N at time t
!! - T  := Y vV + vN,     rhs
!! - vT := invX T,   	get prediction B*
!! - vT := T              copy prediction
!! Pol:
!! - N  := vN,            save N at time t
!! - Np := vNp
!! - P  := Y vV + vN      rhs
!! - Pr := vNp
!! - P,Pr -> P            solve -> prediction
!! - vP := P              copy pred
!------------------------------------------------------------------------
   subroutine vel_predictor(mes_oc, vel_NTor, vel_Npol, rot_omega, vel_Tor, vel_pol)
      use meshs, only : rdom
      use variables, only : var_coll_copy
      use timestep, only : tim_multY, tim_invX, tim_zerobc
      
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(in) :: vel_Ntor !< toroidal field nonlinear term
      type(coll), intent(in) :: vel_npol !< poloidal field nonlinear term
      double precision, intent(in) :: rot_omega !< rotation rate
      type(coll), intent(inout) :: vel_Tor !< toroidal part of field
      type(coll), intent(inout) :: vel_pol !< poloidal part of field
      
      call var_coll_copy (vel_NTor, NTor)
      call tim_multY     (.false.,YTor,vel_Tor,vel_NTor, Tor)
      call vel_setbc_Tor (mes_oc,vel_tor,rot_omega,Tor)
      call tim_invX      (.false.,XTor, Tor)
      call var_coll_copy (Tor, vel_Tor)
      call var_coll_copy (vel_NPol, NPol)
      call var_coll_copy (vel_NPol, Pol)
      call tim_zerobc    (Pol)
      call tim_invX      (.false.,XGre, Pol)
      call tim_multY     (.false.,YPol,vel_Pol,Pol, Pol)
      call tim_zerobc    (Pol)
      call tim_invX      (.false.,XPol, Pol)
      call vel_gre_adjust()
      call var_coll_copy (Pol, vel_Pol)
   end subroutine vel_predictor
   
   
!------------------------------------------------------------------------
!> Velocity Corrector
!!
!! Tor:
!! - T  := c (vN - N),   	using N* get nlin correction
!! - N  := vN
!! - T  := invX T,   	get correction to V (includes factor c)
!! - vT := vT + T		update correction
!! Pol:
!! - P  := c (vN - N),   	using N* get nlin correction
!! - N  := vN
!! - Pr := c (vNp- Np)
!! - P,Pr -> P            solve -> correction
!! - vP := vP + P		update correction
!------------------------------------------------------------------------
   subroutine vel_corrector(mes_oc, vel_NTor, vel_npol, rot_omega,vel_Tor, vel_pol)
      use meshs, only : rdom
      use timestep, only : tim_nlincorr, tim_invX, tim_addcorr, tim_zerobc
      use variables, only : var_coll_copy

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(in) :: vel_Ntor !< toroidal field nonlinear term
      type(coll), intent(in) :: vel_npol !< poloidal field nonlinear term
      double precision, intent(in) :: rot_omega !< rotation rate
      type(coll), intent(inout) :: vel_Tor !< toroidal part of field
      type(coll), intent(inout) :: vel_pol !< poloidal part of field

      call tim_nlincorr  (vel_NTor,NTor, Tor)
      call var_coll_copy (vel_NTor,NTor)
      call vel_setbc_Tor (mes_oc,vel_tor,rot_omega,Tor)
      call tim_invX      (.false.,XTor, Tor)
      call tim_addcorr   (Tor, vel_Tor)
      call tim_nlincorr  (vel_NPol,NPol, Pol)
      call var_coll_copy (vel_NPol,NPol)
      call tim_zerobc    (Pol)
      call tim_invX      (.false.,XGre, Pol)
      call tim_zerobc    (Pol)
      call tim_invX      (.false.,XPol, Pol)
      call vel_gre_adjust()
      call tim_addcorr   (Pol, vel_Pol)

   end subroutine vel_corrector
   

!------------------------------------------------------------------------
!> No information
!------------------------------------------------------------------------
   subroutine vel_gre_adjust()
      double precision :: a, b, BCi,BCo
      integer :: N_
      _loop_lm_vars

      N_ = Pol%D%N
      _loop_lm_begin(Pol)
         if(l==0) cycle
         BCi = Pol%Re( 1,nh)            
         BCo = Pol%Re(N_,nh)            
         a = invG(1,1,l)*BCi + invG(1,2,l)*BCo
         b = invG(2,1,l)*BCi + invG(2,2,l)*BCo
         Pol%Re(:,nh) =  Pol%Re(:,nh) - a*Gre(:,1,l) - b*Gre(:,2,l)
         BCi = Pol%Im( 1,nh)
         BCo = Pol%Im(N_,nh)
         a = invG(1,1,l)*BCi + invG(1,2,l)*BCo
         b = invG(2,1,l)*BCi + invG(2,2,l)*BCo
         Pol%Im(:,nh) =  Pol%Im(:,nh) - a*Gre(:,1,l) - b*Gre(:,2,l)
      _loop_lm_end

   end subroutine vel_gre_adjust


!------------------------------------------------------------------------
!> Removes angular momentum for stress free
!!   
!!   Subroutine vel_angmomclr removes the angular momentum in the stress-free
!!    case by subtracting a small rigid body rotation omz, omy, omx
!!   CAJ 6th Nov 2010
!------------------------------------------------------------------------
   subroutine vel_angmomclr(mes_oc,vel_tor)
      use meshs, only : rdom
      use parameters, only : b_boussinesq, i_vel_bc, d_Pi
      
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type(coll), intent(inout) :: vel_tor !< toroidal part of velocity

      integer :: N_,i,No
      double precision :: w,f1(vel_Tor%D%N),f2(vel_Tor%D%N),Lz,Lz2,omz
      double precision :: Lx,Lx2,omx,Ly,Ly2,omy
      _loop_lm_vars

      if (b_boussinesq) then
     
         if(i_vel_bc/=4) return	!apply to stress-free case only
     
         N_  = vel_Tor%D%N
         f2  = mes_oc%r(1:N_,2)  ! r^2
         w   = 8d0/3d0 * d_PI
         
         _loop_lm_begin(vel_Tor)
                ! Lz(l=1,m=0) = 8*pi/3 * int(r^3*Re(T(r)),dr)
         if (l==1 .and. m_==0) then 
            f1  = vel_Tor%Re(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Re(T(r))
            Lz  = w*dot_product(f1,vel_Tor%D%intr2dr(1:N_)) ! int(r^2*f,dr)
            Lz2 = w*dot_product(f2,vel_Tor%D%intr2dr(1:N_)) ! int(r^2*f2,dr)
            omz = Lz/Lz2
            vel_Tor%Re(1:N_,nh)=vel_Tor%Re(1:N_,nh)-omz*mes_oc%r(1:N_,1)
               ! Lx(l=1,m=1) =  16*pi/3 * int(r^3*Re(T(r)),dr)
         else if(l==1 .and. m_==1) then           
            f1  =  vel_Tor%Re(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Re(T(r))
            Lx  =  2d0*w*dot_product(f1,vel_Tor%D%intr2dr(1:N_)) ! int(r^2*f,dr)
            Lx2 =  2d0*w*dot_product(f2,vel_Tor%D%intr2dr(1:N_)) ! int(r^2*f2,dr)
            f1  =  vel_Tor%Im(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Im(T(r))
            Ly  = -2d0*w*dot_product(f1,vel_Tor%D%intr2dr(1:N_)) ! int(r^2*f,dr)
            Ly2 = -2d0*w*dot_product(f2,vel_Tor%D%intr2dr(1:N_)) ! int(r^2*f2,dr)
            omx = Lx/Lx2
            omy = Ly/Ly2
            vel_Tor%Re(1:N_,nh)=vel_Tor%Re(1:N_,nh)-omx*mes_oc%r(1:N_,1)
            vel_Tor%Im(1:N_,nh)=vel_Tor%Im(1:N_,nh)-omy*mes_oc%r(1:N_,1)
         end if
         _loop_lm_end

      else

         if (i_vel_bc/=4) return  ! only apply if bothh bcs stress-free
         N_  = vel_Tor%D%N
         Lz =0.0
         Lz2=0.0
         Lx =0.0
         Lx2=0.0
         Ly =0.0
         Ly2=0.0

         if (i_vel_bc==4) then
            _loop_lm_begin(vel_Tor)
            if ((l == 1) .and. (m_ <= 1)) then
               w = 8d0/3d0 * d_PI
               if (m_ == 0) then ! Lz(l=1,m=0) = 8*pi/3 * int(r^3*zeta^n*Re(T(r)),dr)
                  f1 =  vel_Tor%Re(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Re(T(r))
                  Lz  = w*dot_product(f1,vel_Tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f,dr)
                  f2 =  mes_oc%r(1:N_,1)*mes_oc%r(1:N_,1)  ! r*r
                  Lz2 = w*dot_product(f2,vel_Tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f2,dr) 
               else ! Lx(l=1,m=1) =  16*pi/3 * int(r^3*zeta^n*Re(T(r)),dr)
                  f1 =  vel_Tor%Re(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Re(T(r))
                  Lx  = 2d0*w*dot_product(f1,vel_Tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f,dr)
                  f2 =  mes_oc%r(1:N_,1)*mes_oc%r(1:N_,1)  ! r*r
                  Lx2 = 2d0*w*dot_product(f2,vel_Tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f2,dr)
                  f1 =  vel_Tor%Im(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Im(T(r))
                  Ly  = -2d0*w*dot_product(f1,vel_Tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f,dr)
                  f2 =  mes_oc%r(1:N_,1)*mes_oc%r(1:N_,1)  ! r*r
                  Ly2 = -2d0*w*dot_product(f2,vel_Tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f2,dr)

               end if
            end if
             
            _loop_lm_end

            omz = Lz/Lz2
            omx = Lx/Lx2
            omy = Ly/Ly2 
            _loop_lm_begin(vel_Tor)
            if ((l == 1) .and. (m_ <= 1)) then
               w = 8d0/3d0 * d_PI
               if (m_ == 0) then
                  do i=1,N_
                     vel_Tor%Re(i,nh)=vel_Tor%Re(i,nh)-omz*mes_oc%r(i,1)
                  enddo
               else 
                  do i=1,N_ 
                     vel_Tor%Re(i,nh)=vel_Tor%Re(i,nh)-omx*mes_oc%r(i,1)
                     vel_Tor%Im(i,nh)=vel_Tor%Im(i,nh)-omy*mes_oc%r(i,1)
                  enddo
               end if
            end if
            _loop_lm_end
         end if
      end if

   end subroutine vel_angmomclr
   
!*************************************************************************
 end module velocity
!*************************************************************************
