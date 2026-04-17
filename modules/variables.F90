#include "../parallel.h"
!*************************************************************************
!>  Variables in spectral,collocated and physical space.
!!  Simple functions on variables.
!!
!!~~~~
!!     r_n    n in [1,N]
!!     th_j   j in [1,Th]
!!     ph_i   i in [0,Ph)
!!~~~~
!!
!!  The harmonics are ordered such that the usual loop is defined as:
!!
!!~~~~
!!     nh = -1
!!     do m = 0, M-1
!!        m_ = m*Mp
!!        do l = m_, L-1
!!           nh = nh + 1
!!           ...
!!~~~~
!! - nh  harmonic index in the array
!! - m   m index
!! - m_  m of P_l^m
!! - l   l of P_l^m
!!
!! See parallel.h
!!*************************************************************************
 module variables
!*************************************************************************
   use parameters, only : i_H1, i_pN, i_M, i_N, i_pH1, i_PTh, i_Ph, i_L1
   use meshs, only : rdom
   implicit none
   save
   
   !> Harmonic data type
   !!
   !!The "harmonic" data type stores the maximum degree and order
   !! in the decomposition of the physical quantities.\n
   !! It also stores information about the processor in which each harmonic is treated.
   type harm
      integer              :: l(0:i_H1) !< spherical harmonic degree
      integer              :: m(0:i_H1) !< spherical harmonic order
      integer              :: s(0:i_H1) !< odd/even symmetry
      integer              :: pH0 !< Index of the first harmonic being treated by the current CPU.
      integer              :: pH1 !< The last harmonic being treated by the current CPU is nh = pH0 + pH1.
      integer              :: pH0_(0: _Nr -1) !< Index of the first harmonic being treated by each CPU.
      integer              :: pH1_(0:_Nr-1) !< The last harmonic being treated by CPU k is nh = pH0_(k) + pH1_(k).
   end type harm

   !> Decomposition information for theta parallelisation
   type thdom
      integer              :: m(0:2* _pHc -1,0:_Nth-1)
      integer              :: pM !< num M points treated by cpu
      integer              :: pM_(0:_Nth-1) !< num M points treated by each th cpu group
      integer              :: pThi !< first theta point treated by processor
      integer              :: pTh !< num theta points treated by processor
      integer              :: pThi_(0:_Nth-1) !< first theta point treated by each th cpu group
      integer              :: pTh_(0:_Nth-1) !< num theta points treated by each th cpu group
   end type thdom

   !> Variable in the spectral domain.\n
   !! Re and Im refer to real and imaginary parts of the field.
   !! First index refers to harmonic index and second to radial shell.
   type spec
      type (rdom), pointer :: D !< The radial domain.
      type (harm), pointer :: H !< The harmonic truncations and CPU distribution.
      double precision     :: Re(0:i_H1, i_pN) !< Real coefficients of the field.
      double precision     :: Im(0:i_H1, i_pN) !< Imaginary coefficients of the field.
   end type spec

   !> A coll (or co-located) variable represents a scalar field or a component of a vector field.
   !! Re and Im are the real and imaginary parts of the field.\n
   !! The first index represents the radial shell and varies between 1 and resolution::i_n .\n
   !! The second index represents the harmonic coefficient and varies between 0 and i_pH1
   type coll
      type (rdom), pointer :: D !< The radial domain.
      type (harm), pointer :: H !< The harmonic truncations and CPU distribution.
      double precision     :: Re(i_N, 0:i_pH1) !< Real coefficients of the field
      double precision     :: Im(i_N, 0:i_pH1) !< Imaginary coefficients of the field
   end type coll

   !> Variables in physical space.
   type phys
      double precision     :: Re(i_pTh, 0:i_Ph-1, i_pN) !< indexes are ordered as (theta, phi, r)
   end type phys
   
   type (harm)               :: var_H
   type (thdom)              :: var_Th
   double precision, private :: ad_ll1(0:i_L1) !< l*(l+1)
   double precision, private :: ad_ll11(0:i_L1) !< 1/[ l*(l+1) ]
   double precision, private :: ad_rtll1(0:i_L1) !< sqrt[ l*(l+1) ]
   double precision, private :: ad_rtll11(0:i_L1) !< 1/( sqrt[l*(l+1)] )
   type (coll),      private :: coll_


 contains


!------------------------------------------------------------------------
!>  Initialise stuff for variables
!------------------------------------------------------------------------
   subroutine var_precompute()
      integer :: l
      ad_ll1(0)    = 0d0
      ad_ll11(0)   = 0d0
      ad_rtll1(0)  = 0d0
      ad_rtll11(0) = 0d0
      do l = 1, i_L1
         ad_ll1(l)    = dble(l*(l+1))
         ad_ll11(l)   = 1d0 / ad_ll1(l)
         ad_rtll1(l)  = dsqrt(ad_ll1(l))
         ad_rtll11(l) = 1d0 / ad_rtll1(l)
      end do
      call var_harm_init(var_H)
      call var_thdom_init(var_Th)
   end subroutine var_precompute


!-------------------------------------------------------------------------
!>  initialise harmonic indices
!-------------------------------------------------------------------------

   subroutine var_harm_init(H)
      use dyn_mpi, only : th_rnk, r_rnk
      use parameters, only : i_H, i_L, i_Mp
      type (harm), intent(out) :: H !< harmonic to initialise
      integer :: n,r,p
     _loop_setlm_vars

      H%pH1_ = -1
      do n = 0, i_H1
         r = _Nr - modulo(n,_Nr) - 1
         H%pH1_(r) = H%pH1_(r) + 1
      end do

      H%pH0_(0) = 0
      do r = 1, _Nr-1
         H%pH0_(r) = H%pH0_(r-1) + H%pH1_(r-1) + 1
      end do 
      p = r_rnk
      H%pH0 = H%pH0_(p)
      H%pH1 = H%pH1_(p)

      ! storing l,m in for parallel groups
      ! s is new to optimise for odd/even symmetry
      _loop_setlm_begin
         H%l(nh) = l
         H%m(nh) = m_
         H%s(nh) = modulo(l+m_,2)
      _loop_setlm_end
   end subroutine var_harm_init



!-------------------------------------------------------------------------
!>  initialise theta domains
!-------------------------------------------------------------------------
   subroutine var_thdom_init(Th)
      use parameters, only : i_Th, i_Mp
      use dyn_mpi

      type (thdom), intent(out) :: Th !< theta domain to initialise
      integer :: n,p
      _loop_sumth_vars

      Th%pThi_(0) = 1
      n = modulo(i_Th,i_pTh)
      if(n==0) Th%pTh_(0) = i_pTh
      if(n/=0) Th%pTh_(0) = n
      do p = 1, _Nth-1
         Th%pThi_(p) = Th%pThi_(p-1) + Th%pTh_(p-1)
         Th%pTh_ (p) = i_pTh
      end do
      p = th_rnk
      Th%pThi = Th%pThi_(p)
      Th%pTh  = Th%pTh_ (p)

      _loop_sumth_begin
         if(m_==m__) cycle
         Th%m(m,p) = m_/i_Mp
      _loop_sumth_end
      Th%pM_(p) = m+1
      Th%pM = m+1
#ifdef _MPI
      do p = 0, _Nth-1
         call mpi_bcast(Th%m(0,p), 2*_pHc, mpi_integer,  &
            p*_Nr, mpi_comm_world, mpi_er)
         call mpi_bcast(Th%pM_(p), 1, mpi_integer,  &
            p*_Nr, mpi_comm_world, mpi_er)
      end do
#endif
      
   end subroutine var_thdom_init


!-------------------------------------------------------------------------
!>  initialise collocated variable
!-------------------------------------------------------------------------
   subroutine var_coll_init(D,H, a)
      type (rdom), target, intent(in) :: D !< outer/inner core domain
      type (harm), target, intent(in) :: H !< harmonic domain
      type (coll),         intent(out) :: a !< scalar field in spectral space (coll)
      a%D => D
      a%H => H
      a%Re = 0d0
      a%Im = 0d0
   end subroutine var_coll_init


!-------------------------------------------------------------------------
!>  initialise spectral variable
!-------------------------------------------------------------------------
   subroutine var_spec_init(D,H, a)
      type (rdom), target, intent(in) :: D !< outer/inner core domain
      type (harm), target, intent(in) :: H !< harmonic domain
      type (spec),         intent(out) :: a !< scalar field in spectral space
      a%D => D
      a%H => H
   end subroutine var_spec_init


!-------------------------------------------------------------------------
!>  Copy a collocated variable
!-------------------------------------------------------------------------
   subroutine var_coll_copy(in, out)
      type (coll), intent(in)  :: in
      type (coll), intent(out) :: out
      integer :: n, nhm
      out%D => in%D
      out%H => in%H
      n   = in%D%N
      nhm = in%H%pH1
      out%Re(1:n,0:nhm) = in%Re(1:n,0:nhm)
      out%Im(1:n,0:nhm) = in%Im(1:n,0:nhm)
   end subroutine var_coll_copy


!-------------------------------------------------------------------------
!>  convert collocated -> spectral
!!
!! Transpose data (serial)
!! Parallel transpose over radial mpi group (mpi)
!-------------------------------------------------------------------------
#ifndef _MPI
   subroutine var_coll2spec(c,s, c2,s2, c3,s3)
      type (coll), intent(in)  :: c
      type (spec), intent(out) :: s
      type (coll), intent(in),  optional :: c2,c3
      type (spec), intent(out), optional :: s2,s3
      integer :: n, nhm
      s%D => c%D 
      s%H => c%H
      n   = s%D%N
      nhm = s%H%pH1
      s%Re(0:nhm,1:n) = transpose(c%Re(1:n,0:nhm))
      s%Im(0:nhm,1:n) = transpose(c%Im(1:n,0:nhm))
      if(.not. present(c2)) return
      s2%D => c2%D 
      s2%H => c2%H
      s2%Re(0:nhm,1:n) = transpose(c2%Re(1:n,0:nhm))
      s2%Im(0:nhm,1:n) = transpose(c2%Im(1:n,0:nhm))
      if(.not. present(c3)) return
      s3%D => c3%D 
      s3%H => c3%H
      s3%Re(0:nhm,1:n) = transpose(c3%Re(1:n,0:nhm))
      s3%Im(0:nhm,1:n) = transpose(c3%Im(1:n,0:nhm))      
   end subroutine var_coll2spec

#else
   subroutine var_coll2spec(c,s, c2,s2, c3,s3)
      use dyn_mpi
      type (coll), intent(in)  :: c
      type (spec), intent(out) :: s
      type (coll), intent(in),  optional :: c2,c3
      type (spec), intent(out), optional :: s2,s3
      double precision :: bsend(2*i_pN*(i_pH1+1)*3,0:_Nr-1)
      double precision :: brecv(2*i_pN*(i_pH1+1)*3,0:_Nr-1)
      integer :: stp, dst,src, n,nh,l, nc
      integer :: ic
      nc = 1

      s%D => c%D 
      s%H => c%H

      if(present(c2)) nc = 2
      if(nc==2) s2%D => c2%D
      if(nc==2) s2%H => c2%H
      if(present(c3)) nc = 3
      if(nc==3) s3%D => c3%D
      if(nc==3) s3%H => c3%H

      do stp = 0, _Nr-1
         src  = modulo(_Nr-stp+r_rnk, _Nr)         
         mpi_tg = stp
         call mpi_irecv( brecv(1,stp), 2*c%D%pN*(c%H%pH1_(src)+1)*nc,  &
            mpi_double_precision, src, mpi_tg, r_comm,  &
            mpi_rq(stp), mpi_er)
      end do


      do stp = 0, _Nr-1
         dst  = modulo(stp+r_rnk, _Nr)
         l = 1
         do n = c%D%pNi_(dst), c%D%pNi_(dst)+c%D%pN_(dst)-1
            do nh = 0, c%H%pH1
               bsend(l,  stp) = c%Re(n,nh)
               bsend(l+1,stp) = c%Im(n,nh)
               l = l + 2
               if(nc<2) cycle
               bsend(l,  stp) = c2%Re(n,nh)
               bsend(l+1,stp) = c2%Im(n,nh)
               l = l + 2
               if(nc<3) cycle
               bsend(l,  stp) = c3%Re(n,nh)
               bsend(l+1,stp) = c3%Im(n,nh)
               l = l + 2
            end do
         end do
         mpi_tg = stp

         call mpi_isend( bsend(1,stp), 2*c%D%pN_(dst)*(c%H%pH1+1)*nc,  &
            mpi_double_precision, dst, mpi_tg, r_comm,  &
            mpi_rq(_Nr+stp), mpi_er)
      end do

      do stp = 0, _Nr-1
         src  = modulo(_Nr-stp+r_rnk, _Nr)      
     
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
       do n = 1, s%D%pN
          do nh = s%H%pH0_(src), s%H%pH0_(src)+s%H%pH1_(src)
               s%Re(nh,n) = brecv(l,  stp)
               s%Im(nh,n) = brecv(l+1,stp)
               l = l + 2
               if(nc<2) cycle
               s2%Re(nh,n) = brecv(l,  stp)
               s2%Im(nh,n) = brecv(l+1,stp)
               l = l + 2
               if(nc<3) cycle
               s3%Re(nh,n) = brecv(l,  stp)
               s3%Im(nh,n) = brecv(l+1,stp)
               l = l + 2
            end do
         end do

      end do      

      do stp = 0, _Nr-1
         call mpi_wait( mpi_rq(_Nr+stp), mpi_st, mpi_er)
      end do

   end subroutine var_coll2spec
#endif
   
!-------------------------------------------------------------------------
!>  convert spectral -> collocated
!!
!! Transpose data (serial)
!! Parallel transpose over radial mpi group (mpi)
!! s -> c, s2 -> c2 (optional), s3 -> c3 (optional) 
!-------------------------------------------------------------------------
#ifndef _MPI
   subroutine var_spec2coll(s,c, s2,c2, s3,c3)
      type (spec), intent(in)  :: s !< s -> c
      type (coll), intent(out) :: c !< s -> c
      type (spec), intent(in),  optional :: s2 !< s2 -> c2
      type (coll), intent(out), optional :: c2 !< s2 -> c2
      type (spec), intent(in),  optional :: s3 !< s3 -> c3
      type (coll), intent(out), optional :: c3 !< s3 -> c3
      integer :: n, nhm
      c%D => s%D
      c%H => s%H
      n   = s%D%N
      nhm = s%H%pH1
      c%Re(1:n,0:nhm) = transpose(s%Re(0:nhm,1:n))
      c%Im(1:n,0:nhm) = transpose(s%Im(0:nhm,1:n))
      if(.not. present(s2)) return
      c2%D => s2%D
      c2%H => s2%H
      c2%Re(1:n,0:nhm) = transpose(s2%Re(0:nhm,1:n))
      c2%Im(1:n,0:nhm) = transpose(s2%Im(0:nhm,1:n))
      if(.not. present(s3)) return
      c3%D => s3%D
      c3%H => s3%H
      c3%Re(1:n,0:nhm) = transpose(s3%Re(0:nhm,1:n))
      c3%Im(1:n,0:nhm) = transpose(s3%Im(0:nhm,1:n))
   end subroutine var_spec2coll

#else
   subroutine var_spec2coll(s,c, s2,c2, s3,c3)
      use dyn_mpi
      type (spec), intent(in)  :: s
      type (coll), intent(out) :: c
      type (spec), intent(in),  optional :: s2,s3
      type (coll), intent(out), optional :: c2,c3
      double precision :: bsend(2*(i_pH1+1)*i_pN*3,0:_Nr-1)
      double precision :: brecv(2*(i_pH1+1)*i_pN*3,0:_Nr-1)
      integer :: stp, dst,src, n,nh,l, ns

      ns = 1
      c%D => s%D 
      c%H => s%H
      if(present(s2)) ns = 2
      if(ns==2) c2%D => s2%D
      if(ns==2) c2%H => s2%H
      if(present(s3)) ns = 3
      if(ns==3) c3%D => s3%D
      if(ns==3) c3%H => s3%H
      
      do stp = 0, _Nr-1
         src  = modulo(_Nr-stp+r_rnk, _Nr)
         mpi_tg = stp
         call mpi_irecv( brecv(1,stp), 2*(s%H%pH1+1)*s%D%pN_(src)*ns,  &
            mpi_double_precision, src, mpi_tg, r_comm,  &
            mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Nr-1
         dst  = modulo(stp+r_rnk, _Nr)
         l = 1
         do nh = s%H%pH0_(dst), s%H%pH0_(dst)+s%H%pH1_(dst)
            do n = 1, s%D%pN
               bsend(l,  stp) = s%Re(nh,n)
               bsend(l+1,stp) = s%Im(nh,n)
               l = l + 2
               if(ns<2) cycle
               bsend(l,  stp) = s2%Re(nh,n)
               bsend(l+1,stp) = s2%Im(nh,n)
               l = l + 2
               if(ns<3) cycle
               bsend(l,  stp) = s3%Re(nh,n)
               bsend(l+1,stp) = s3%Im(nh,n)
               l = l + 2
            end do
         end do
         mpi_tg = stp
         call mpi_isend( bsend(1,stp), 2*(s%H%pH1_(dst)+1)*s%D%pN*ns,  &
            mpi_double_precision, dst, mpi_tg, r_comm,  &
            mpi_rq(_Nr+stp), mpi_er)
      end do
      
      do stp = 0, _Nr-1
         src  = modulo(_Nr-stp+r_rnk, _Nr)
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do nh = 0, c%H%pH1
            do n = c%D%pNi_(src), c%D%pNi_(src)+c%D%pN_(src)-1
               c%Re(n,nh) = brecv(l,  stp)
               c%Im(n,nh) = brecv(l+1,stp)
               l = l + 2
               if(ns<2) cycle
               c2%Re(n,nh) = brecv(l,  stp)
               c2%Im(n,nh) = brecv(l+1,stp)
               l = l + 2
               if(ns<3) cycle
               c3%Re(n,nh) = brecv(l,  stp)
               c3%Im(n,nh) = brecv(l+1,stp)
               l = l + 2
            end do
         end do
      end do

      do stp = 0, _Nr-1
         call mpi_wait( mpi_rq(_Nr+stp), mpi_st, mpi_er)
      end do

   end subroutine var_spec2coll
#endif

!------------------------------------------------------------------------
!>  Multiply a collocated variable by a mesh-type
!!     out = A.in
!------------------------------------------------------------------------
   subroutine var_coll_meshmult(A,in, out)
      use meshs, only : mesh
      use parameters, only : i_Kl

      type (mesh), intent(in)  :: A !< mesh, A in "out=A.in"
      type (coll), intent(in)  :: in !< coll, in in "out=A.in"
      type (coll), intent(out) :: out !< coll, out in "out=A.in"
      double precision :: re(i_N), im(i_N)
      integer :: N_,nh, n,j,l,r

      out%D => in%D
      out%H => in%H
      N_  = in%D%N

      do nh = 0, in%H%pH1
         re = 0d0
         im = 0d0
         do j = 1, N_
            l = max(1,j-i_KL)
            r = min(j+i_KL,N_)
            do n = l, r
               re(n) = re(n) + A%M(i_KL+1+n-j,j) * in%Re(j,nh)
               im(n) = im(n) + A%M(i_KL+1+n-j,j) * in%Im(j,nh)
            end do
         end do
         out%Re(1:N_,nh) = re(1:N_)
         out%Im(1:N_,nh) = im(1:N_)
      end do

   end subroutine var_coll_meshmult


!------------------------------------------------------------------------
!>  change to qst form  LM: Pol->q,s, Tor->t
!------------------------------------------------------------------------
   subroutine var_coll_TorPol2qst(Tor,Pol, q,s,t)
!
!   This is now the compressible version of var_coll_TorPol2qst. The Boussinesq
!   version is now var_coll_TorPol2qst_mag below.
! 
      type (coll), intent(in)  :: Tor !< toroidal part of field, spectral space
      type (coll), intent(in)  :: Pol !< poloidal part of field, spctral space
      type (coll), intent(out) :: q !< q, part of (q,s,t) decomposition
      type (coll), intent(out) :: s !< s, part of (q,s,t) decomposition
      type (coll), intent(out) :: t !< t, part of (q,s,t) decomposition
      double precision :: a
      integer :: N_,n
      _loop_lm_vars
         			! q = l(l+1)/r Pol
      q%D =>Pol%D
      q%H =>Pol%H
      N_  = Pol%D%N
      _loop_lm_begin(q)
         do n = 1, N_
            a = ad_ll1(l) * q%D%r(n,-1)
            q%Re(n,nh) = Pol%Re(n,nh) * a
            q%Im(n,nh) = Pol%Im(n,nh) * a
         end do
      _loop_lm_end
! Modified for compressiblity here
!
!  s = sqrt(l(l+1)) (1/r+xi+dr) Pol
!
      call var_coll_meshmult(Pol%D%compr1dr,Pol, s)
      _loop_lm_begin(s)
         do n = 1, N_
            s%Re(n,nh) = s%Re(n,nh) * ad_rtll1(l)
            s%Im(n,nh) = s%Im(n,nh) * ad_rtll1(l)
         end do
      _loop_lm_end
         			! t = - sqrt(l(l+1)) Tor
      t%D =>Tor%D
      t%H =>Tor%H
      _loop_lm_begin(t)
         do n = 1, N_
            t%Re(n,nh) = - Tor%Re(n,nh) * ad_rtll1(l)
            t%Im(n,nh) = - Tor%Im(n,nh) * ad_rtll1(l)
         end do
      _loop_lm_end

   end subroutine var_coll_TorPol2qst

! ## Piotr ## (14/07/2009)
! This is original Ashley's version of var_coll_TorPol2qst
! (using definition F = curl(rT)+curl(curl(rP))
! to be used for magnetic field
!------------------------------------------------------------------------
!>  change to qst form  LM (for magnetic field): Pol->q,s, Tor->t
!------------------------------------------------------------------------
   subroutine var_coll_TorPol2qst_mag(Tor,Pol, q,s,t)


      type (coll), intent(in)  :: Tor !< toroidal part of field, spectral space
      type (coll), intent(in)  :: Pol !< poloidal part of field, spctral space
      type (coll), intent(out) :: q !< q, part of (q,s,t) decomposition
      type (coll), intent(out) :: s !< s, part of (q,s,t) decomposition
      type (coll), intent(out) :: t !< t, part of (q,s,t) decomposition
      double precision :: a
      integer :: N_,n
      _loop_lm_vars
                                ! q = l(l+1)/r Pol
      q%D =>Pol%D
      q%H =>Pol%H
      N_  = Pol%D%N
      _loop_lm_begin(q)
         do n = 1, N_
            a = ad_ll1(l) * q%D%r(n,-1)
            q%Re(n,nh) = Pol%Re(n,nh) * a
            q%Im(n,nh) = Pol%Im(n,nh) * a
         end do
      _loop_lm_end
                                ! s = sqrt(l(l+1)) (1/r+dr) Pol
      call var_coll_meshmult(Pol%D%r1dr,Pol, s)
      _loop_lm_begin(s)
         do n = 1, N_
            s%Re(n,nh) = s%Re(n,nh) * ad_rtll1(l)
            s%Im(n,nh) = s%Im(n,nh) * ad_rtll1(l)
         end do
      _loop_lm_end
                                ! t = - sqrt(l(l+1)) Tor
      t%D =>Tor%D
      t%H =>Tor%H
      _loop_lm_begin(t)
         do n = 1, N_
            t%Re(n,nh) = - Tor%Re(n,nh) * ad_rtll1(l)
            t%Im(n,nh) = - Tor%Im(n,nh) * ad_rtll1(l)
         end do
      _loop_lm_end

   end subroutine var_coll_TorPol2qst_mag

!-----------------------------------------------------------
!> Derivatives used for compressible case
!-----------------------------------------------------------
   subroutine var_coll_TorPol2qstderiv(Tor,Pol, q,s,t,qr,sr,tr)
!
!  This is only used in the compressible case 

      type (coll), intent(in)  :: Tor !< toroidal part of field, spectral space
      type (coll), intent(in)  :: Pol !< poloidal part of field, spctral space
      type (coll), intent(out) :: q !< q, part of (q,s,t) decomposition
      type (coll), intent(out) :: s !< s, part of (q,s,t) decomposition
      type (coll), intent(out) :: t !< t, part of (q,s,t) decomposition
      type (coll), intent(out) :: qr !< dq/dr, q: part of (q,s,t) decomposition
      type (coll), intent(out) :: sr !< ds/dr, s: part of (q,s,t) decomposition
      type (coll), intent(out) :: tr !< dt/dr, t: part of (q,s,t) decomposition
! calculates dq/dr, ds/dr and dt/dr as well as q,s and t
      double precision :: a
      integer :: N_,n
      _loop_lm_vars
                                ! q = l(l+1)/r Pol, qr=dq/dr
      q%D =>Pol%D
      q%H =>Pol%H
      N_  = Pol%D%N
      _loop_lm_begin(q)
         do n = 1, N_
            a = ad_ll1(l) * q%D%r(n,-1)
            q%Re(n,nh) = Pol%Re(n,nh) * a
            q%Im(n,nh) = Pol%Im(n,nh) * a
         end do
      _loop_lm_end
      call var_coll_meshmult(q%D%dr(1),q,qr)
                                ! s = sqrt(l(l+1)) (1/r+dr) Pol
 
      call var_coll_meshmult(Pol%D%compr1dr,Pol, s)
      _loop_lm_begin(s)
         do n = 1, N_
            s%Re(n,nh) = s%Re(n,nh) * ad_rtll1(l)
            s%Im(n,nh) = s%Im(n,nh) * ad_rtll1(l)
         end do
      _loop_lm_end
      call var_coll_meshmult(s%D%dr(1),s,sr)
 
                                ! t = - sqrt(l(l+1)) Tor
      t%D =>Tor%D
      t%H =>Tor%H
      _loop_lm_begin(t)
         do n = 1, N_
            t%Re(n,nh) = - Tor%Re(n,nh) * ad_rtll1(l)
            t%Im(n,nh) = - Tor%Im(n,nh) * ad_rtll1(l)
         end do
      _loop_lm_end
      call var_coll_meshmult(t%D%dr(1),t,tr) 

   end subroutine var_coll_TorPol2qstderiv



!------------------------------------------------------------------------
!>  change qst to Tor,Pol form  LM: q->Pol, t->Tor
!!
!!  Only really valid for solenoindal vectors, only
!!  q,t required
!------------------------------------------------------------------------
   subroutine var_coll_qst2TorPol(q,t, Tor,Pol)
      type (coll), intent(in)  :: q !< q, part of (q,s,t) decomposition
      type (coll), intent(in)  :: t !< t, part of (q,s,t) decomposition
      type (coll), intent(out) :: Tor !< toroidal part of field
      type (coll), intent(out) :: Pol !< poloidal part of field
      double precision :: a
      integer :: N_,n
      _loop_lm_vars

         			! q = l(l+1)/r Pol
      Pol%D =>q%D
      Pol%H =>q%H
      N_    = q%D%N
      _loop_lm_begin(q)
         do n = 1, N_
            a = q%D%r(n,1) * ad_ll11(l)
            Pol%Re(n,nh) = q%Re(n,nh) * a
            Pol%Im(n,nh) = q%Im(n,nh) * a
         end do
      _loop_lm_end
         			! t = - sqrt(l(l+1)) Tor
      Tor%D =>t%D
      Tor%H =>t%H
      _loop_lm_begin(t)
         do n = 1, N_
            Tor%Re(n,nh) = - t%Re(n,nh) * ad_rtll11(l)
            Tor%Im(n,nh) = - t%Im(n,nh) * ad_rtll11(l)
         end do
      _loop_lm_end

   end subroutine var_coll_qst2TorPol


!------------------------------------------------------------------------
!>  take the curl of a vector in qst form;
!!
!!  curl (q, s, t)
!!   = (  -sqrt(l(l+1))t/r , 
!!        -(1/r) * d_r(rt) ,
!!          (1/r) * d_r(rs) - sqrt(l(l+1))q/r
!!      )
!!  change to qst form  LM:  q,s->ot,  t->oq,os
!------------------------------------------------------------------------
   subroutine var_coll_qstcurl(q,s,t, oq,os,ot)
      type (coll), intent(in)  :: q !< q part of (q, s, t) decomposition
      type (coll), intent(in)  :: s !< s part of (q, s, t) decomposition
      type (coll), intent(in)  :: t !< t part of (q, s, t) decomposition
      type (coll), intent(out) :: oq !< -sqrt(l(l+1))t/r
      type (coll), intent(out) :: os !< -(1/r) * d_r(rt) 
      type (coll), intent(out) :: ot !< (1/r) * d_r(rs) - sqrt(l(l+1))q/r
      double precision :: a
      integer :: N_,n
      _loop_lm_vars

      call var_coll_copy(t,coll_)
         			! t := (1/r+dr)s - sqrt(l(l+1))/r q
      call var_coll_meshmult(s%D%r1dr,s, ot)
      N_  = ot%D%N
      _loop_lm_begin(ot)
         do n = 1, N_
            a = ad_rtll1(l) * ot%D%r(n,-1)
            ot%Re(n,nh) = ot%Re(n,nh) - a*q%Re(n,nh)
            ot%Im(n,nh) = ot%Im(n,nh) - a*q%Im(n,nh)
         end do
      _loop_lm_end
         			! s := - (1/r+dr)t
      call var_coll_meshmult(coll_%D%r1dr,coll_, os)
      _loop_lm_begin(os)
         do n = 1, N_
            os%Re(n,nh) =  - os%Re(n,nh)
            os%Im(n,nh) =  - os%Im(n,nh)
         end do
      _loop_lm_end
         			! q := - sqrt(l(l+1))/r t
      oq%D =>coll_%D
      oq%H =>coll_%H
      _loop_lm_begin(oq)
         do n = 1, N_
            a = - ad_rtll1(l) * oq%D%r(n,-1)
            oq%Re(n,nh) = coll_%Re(n,nh)*a
            oq%Im(n,nh) = coll_%Im(n,nh)*a
         end do
      _loop_lm_end

   end subroutine var_coll_qstcurl


!------------------------------------------------------------------------
!>  \f$r/(l(l+1)) \hat{r} \cdot curl A\f$  where \f$A\f$ is in qst form;
!!
!! only the toroidal (t) part of \f$A\f$ is required.
!------------------------------------------------------------------------
   subroutine var_coll_qstllcurlr(t, r)
      type (coll), intent(in)  :: t !< toroidal part of field
      type (coll), intent(out) :: r !<  r/(l(l+1)) * curl(A)
      double precision :: a
      integer :: N_,n
      _loop_lm_vars
         			!  - (r/(l(l+1))) * sqrt(l(l+1))/r t
      r%D =>t%D
      r%H =>t%H
      N_  = t%D%N
      _loop_lm_begin(t)
         do n = 1, N_
            r%Re(n,nh) =  - t%Re(n,nh) * ad_rtll11(l)
            r%Im(n,nh) =  - t%Im(n,nh) * ad_rtll11(l)
         end do
      _loop_lm_end

   end subroutine var_coll_qstllcurlr


!------------------------------------------------------------------------
!>  \f$r/(l(l+1)) \hat{r} \cdot \nabla\times\nabla\times \vec{A}\f$
!!  where A is in qst form;  only q,s required
!------------------------------------------------------------------------
   subroutine var_coll_qstllcurlcurlr(q,s, r)
      type (coll), intent(in)  :: q !< q part of A(q, s, t) decomposition
      type (coll), intent(in)  :: s !< s part of A(q, s, t) decomposition
      type (coll), intent(out) :: r !< r(l(l+1))r * curl curl A(q,s,t)
      integer :: N_,n
      _loop_lm_vars
         			! r/(l(l+1)) * ...
         			! l(l+1)/r^2 q - sqrt(l(l+1))/r (1/r+dr)s
      call var_coll_meshmult(s%D%r1dr,s, r)
      N_  = q%D%N
      _loop_lm_begin(q)
         do n = 1, N_
            r%Re(n,nh) = r%Re(n,nh)*ad_rtll11(l)
            r%Re(n,nh) = q%Re(n,nh)*r%D%r(n,-1) - r%Re(n,nh)
            r%Im(n,nh) = r%Im(n,nh)*ad_rtll11(l)
            r%Im(n,nh) = q%Im(n,nh)*r%D%r(n,-1) - r%Im(n,nh)
         end do
      _loop_lm_end

   end subroutine var_coll_qstllcurlcurlr


!------------------------------------------------------------------------
!>  phi-derivative
!------------------------------------------------------------------------
   subroutine var_coll_dph(in, out)
      type (coll), intent(in)  :: in !< scalar field in spectral space
      type (coll), intent(out) :: out !< d (in) / d phi
      double precision :: a
      integer :: N_,n
      _loop_lm_vars

      out%D =>in%D
      out%H =>in%H
      N_  = in%D%N
      _loop_lm_begin(out)
         do n = 1, N_
            a = in%Re(n,nh)
            out%Re(n,nh) =  - in%Im(n,nh) * m_
            out%Im(n,nh) =    a * m_
         end do
      _loop_lm_end

   end subroutine var_coll_dph


!------------------------------------------------------------------------
!>  norm  \int a.a zeta^n dV
!------------------------------------------------------------------------
! WD 25/09/2013, i replaced the El1 (l=1, energy) with the equatorial
! symmetric and axisymmetric (m=0,l=even)
   subroutine var_coll_norm_comp(a, Eodd,Eeven,Em0le,Em0)
      use parameters, only : d_Pi
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field in spectral space
      double precision, intent(out) :: Eodd !< int(a.a*zeta^n) dV (for (l-m)=odd)
      double precision, intent(out) :: Eeven !< int(a.a*zeta^n) dV (for (l-m)=even)
      double precision, intent(out) :: Em0le !< int(a.a*zeta^n) dV (for m=0, l:even)
      double precision, intent(out) :: Em0 !< int(a.a*zeta^n) dV (for m=0)
      double precision :: w, f(a%D%N), Fev(a%D%N), Fod(a%D%N)
      double precision :: Fm0le(a%D%N), Fm0(a%D%N)
      integer :: N_,n
      _loop_lm_vars

      Fod   = 0d0
      Fev   = 0d0
      Fm0le = 0d0
      Fm0   = 0d0
      N_    = a%D%N

      ! Integral over sphere
      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         if(m_/=0) w = 4d0*w
         f =  w*( a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         if(m_==0)  &
            Fm0 = Fm0 + f
         if((m_==0) .and. (modulo(l,2)==0))  &
            Fm0le = Fm0le + f
         if(modulo(l-m_,2)==1) then
            Fod = Fod + f
         else
            Fev = Fev + f
         end if
      _loop_lm_end

      ! Integral over radius
      Eodd  = dot_product(Fod,a%D%intr2dr_comp(1:N_))
      Eeven = dot_product(Fev,a%D%intr2dr_comp(1:N_))
      Em0le = dot_product(Fm0le,a%D%intr2dr_comp(1:N_))
      Em0   = dot_product(Fm0,a%D%intr2dr_comp(1:N_))
#ifdef _MPI
      call mpi_allreduce( Eodd,  w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eodd  = w
      call mpi_allreduce( Eeven, w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eeven = w
      call mpi_allreduce( Em0le,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em0le   = w
      call mpi_allreduce( Em0,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em0   = w
#endif

   end subroutine var_coll_norm_comp



!------------------------------------------------------------------------
! WD 25/09/2013, i replaced the El1 (l=1, energy) with the equatorially
! symmetric and axisymmetric (m=0,l=even)
!------------------------------------------------------------------------
!>  Computes \f$ \int a.a dV \f$
!------------------------------------------------------------------------
   subroutine var_coll_norm(a, Eodd,Eeven,Em0le,Em0)
      use parameters, only : d_Pi
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field in spectral space
      double precision, intent(out) :: Eodd !< int( a.a ) dV (for (l-m)=odd)
      double precision, intent(out) :: Eeven !< int( a.a ) dV (for (l-m)=even)
      double precision, intent(out) :: Em0le !< int( a.a) dV (for m=0, l:even)
      double precision, intent(out) :: Em0 !< int( a.a ) dV (for m=0)
      double precision :: w, f(a%D%N), Fev(a%D%N), Fod(a%D%N)
      double precision :: Fm0le(a%D%N), Fm0(a%D%N)
      integer :: N_,n
      _loop_lm_vars
      
      Fod   = 0d0
      Fev   = 0d0
      Fm0le = 0d0
      Fm0   = 0d0
      N_    = a%D%N
      
      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         if(m_/=0) w = 4d0*w
         f =  w*( a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         if(m_==0)  &
            Fm0 = Fm0 + f
         if((m_==0) .and. (modulo(l,2)==0))  &
            Fm0le = Fm0le + f
         if(modulo(l-m_,2)==1) then
            Fod = Fod + f
         else
            Fev = Fev + f
         end if
      _loop_lm_end

      Eodd  = dot_product(Fod,a%D%intr2dr(1:N_))
      Eeven = dot_product(Fev,a%D%intr2dr(1:N_))
      Em0le = dot_product(Fm0le,a%D%intr2dr(1:N_))
      Em0   = dot_product(Fm0,a%D%intr2dr(1:N_))
#ifdef _MPI
      call mpi_allreduce( Eodd,  w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eodd  = w
      call mpi_allreduce( Eeven, w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eeven = w
      call mpi_allreduce( Em0le,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em0le   = w
      call mpi_allreduce( Em0,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em0   = w
#endif

   end subroutine var_coll_norm
   
!------------------------------------------------------------------------
!>  norm  \int a.a dV
!!
!! Original version of var_coll_norm from LSD version of code.
!! Outputs are integrals for  odd, even, l=1, and m=0.
!! Used to compute mag_cmb data.
!------------------------------------------------------------------------
   subroutine var_coll_norm_orig(a, Eodd, Eeven, El1, Em0, El12)
      use parameters, only : d_Pi
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field in spectral space
      double precision, intent(out) :: Eodd !< int( a.a ) dV (for (l-m)=odd)
      double precision, intent(out) :: Eeven !< int( a.a ) dV (for (l-m)=even)
      double precision, intent(out) :: El1 !< int( a.a ) dV for (for l=1)
      double precision, intent(out) :: Em0 !< int( a.a ) dV for (m=0)
      double precision, intent(out) :: El12 !< int( a.a ) dV (for 1 <= l <= 12)
      double precision :: w, f(a%D%N), Fev(a%D%N), Fod(a%D%N)
      double precision :: Fl1(a%D%N), Fm0(a%D%N), Fl12(a%D%N)
      integer :: N_,n
      _loop_lm_vars
      
      Fod = 0d0
      Fev = 0d0
      Fl1 = 0d0
      Fm0 = 0d0
      Fl12 = 0d0
      N_  = a%D%N
      
      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         if(m_/=0) w = 4d0*w
         f =  w*( a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         if(m_==0)  &
            Fm0 = Fm0 + f
         if(l>0 .and. l<=12) &
            Fl12 = Fl12 + f
         if(l==1)  &
            Fl1 = Fl1 + f
         if(modulo(l-m_,2)==1) then
            Fod = Fod + f
         else
            Fev = Fev + f
         end if
      _loop_lm_end

      Eodd  = dot_product(Fod,a%D%intr2dr(1:N_))
      Eeven = dot_product(Fev,a%D%intr2dr(1:N_))
      El1   = dot_product(Fl1,a%D%intr2dr(1:N_))
      Em0   = dot_product(Fm0,a%D%intr2dr(1:N_))
      El12   = dot_product(Fl12,a%D%intr2dr(1:N_))
#ifdef _MPI
      call mpi_allreduce( Eodd,  w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eodd  = w
      call mpi_allreduce( Eeven, w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eeven = w
      call mpi_allreduce( El1,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      El1   = w
      call mpi_allreduce( El12,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      El12   = w
      call mpi_allreduce( Em0,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em0   = w
#endif

   end subroutine var_coll_norm_orig

!-----------------------------------------------------------------
!> Computes the extra term in the viscous dissipation integral when
!> non-Boussinesq: int_v |q|^2 2 rho ( d xi/dr - (1/3) xi^2) dv 
!------------------------------------------------------------------
   subroutine var_coll_norm_veldiss(a, Ecomp)
      use parameters, only : d_Pi
      use dyn_mpi
      type (coll), intent(in) :: a !< q part of (qst) decomposition of velocity
      double precision, intent(out) :: Ecomp !< int q^2 rho dV
      double precision :: w, f(a%D%N), Fcomp(a%D%N)
      integer :: N_,n
      _loop_lm_vars

      Fcomp   = 0d0
      N_    = a%D%N

      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         if(m_/=0) w = 4d0*w
         f =  w*( a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         Fcomp = Fcomp + f
      _loop_lm_end

      Ecomp = dot_product(Fcomp,a%D%intr2dr_veldiss(1:N_))
#ifdef _MPI
      call mpi_allreduce( Ecomp,  w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Ecomp  = w
#endif

   end subroutine var_coll_norm_veldiss

!---------------------------------------------------------
!> int( a.a r^2 eta) dV 
!!
!! Calculation of magnetic dissipation
!---------------------------------------------------------

    subroutine var_coll_norm_magdiss(a, Eodd,Eeven,Em0le,Em0)
      use parameters, only : d_Pi
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field in spectral space
      double precision, intent(out) :: Eodd !< int( a.a r^2 eta ) dV (for (l-m)=odd)
      double precision, intent(out) :: Eeven !< int( a.a r^2 eta ) dV (for (l-m)=even)
      double precision, intent(out) :: Em0le !< int( a.a r^2 eta ) dV (for m=0, l:even)
      double precision, intent(out) :: Em0 !< int( a.a r^2 eta ) dV (for m=0)
      double precision :: w, f(a%D%N), Fev(a%D%N), Fod(a%D%N)
      double precision :: Fm0le(a%D%N), Fm0(a%D%N)
      integer :: N_,n
      _loop_lm_vars

      Fod   = 0d0
      Fev   = 0d0
      Fm0le = 0d0
      Fm0   = 0d0
      N_    = a%D%N

      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         if(m_/=0) w = 4d0*w
         f =  w*( a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         if(m_==0)  &
            Fm0 = Fm0 + f
         if((m_==0) .and. (modulo(l,2)==0))  &
            Fm0le = Fm0le + f
         if(modulo(l-m_,2)==1) then
            Fod = Fod + f
         else
            Fev = Fev + f
         end if
      _loop_lm_end

      Eodd  = dot_product(Fod,a%D%intr2dr_magdiss(1:N_))
      Eeven = dot_product(Fev,a%D%intr2dr_magdiss(1:N_))
      Em0le = dot_product(Fm0le,a%D%intr2dr_magdiss(1:N_))
      Em0   = dot_product(Fm0,a%D%intr2dr_magdiss(1:N_))
#ifdef _MPI
      call mpi_allreduce( Eodd,  w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eodd  = w
      call mpi_allreduce( Eeven, w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eeven = w
      call mpi_allreduce( Em0le,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em0le   = w
      call mpi_allreduce( Em0,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em0   = w
#endif

   end subroutine var_coll_norm_magdiss

!-------------------------------------------------------------------
!!> var_coll_norm_surf 
!!
!! evaluates
!! int 0_2pi int 0_pi rho r |a|^2  sin theta dtheta dphi on the
!! inner and outer surfaces
!-------------------------------------------------------------------
   subroutine var_coll_norm_surf(a, ESin, ESout)
      use parameters, only : d_Pi, b_boussinesq
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field, spectral space
      double precision, intent(out) :: ESin !< int rho r a.a dS, inner boundary
      double precision, intent(out) :: ESout !< int rho r a.a dS, outer boundary
      double precision :: w, fin, fout, FSin, FSout
      integer :: N_,n
      _loop_lm_vars

      FSin  = 0d0
      FSout = 0d0
      N_    = a%D%N

      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         if (m_/=0) w = 4d0*w
         fin =  w*( a%Re(1,nh)*a%Re(1,nh)  &
                + a%Im(1,nh)*a%Im(1,nh) )
         fout =  w*( a%Re(N_,nh)*a%Re(N_,nh)  &
                + a%Im(N_,nh)*a%Im(N_,nh) )
         FSin  = FSin + fin
         FSout = FSout + fout
      _loop_lm_end

      if (.not.b_boussinesq) then
        FSin = FSin * a%D%rho(1,0)*a%D%r(1,1)
        FSout= FSout * a%D%rho(N_,0)*a%D%r(N_,1)
      else
        FSin = FSin * a%D%r(1,1)
        FSout= FSout * a%D%r(N_,1)
      end if 
#ifdef _MPI
      call mpi_allreduce( FSin,  w, 1, mpi_double_precision,  &
        mpi_sum, mpi_comm_world, mpi_er)
       ESin  = w
       call mpi_allreduce( FSout, w, 1, mpi_double_precision,  &
          mpi_sum, mpi_comm_world, mpi_er)
       ESout = w
#endif

    end subroutine var_coll_norm_surf


   ! ## Piotr ## (14/07/2009)
! This is compressible version of angular momentum
! (to be used with velocity toroidal potential)
!------------------------------------------------------------------------
!> Compressible version of angular momentum
!!
!! This is compressible version of angular momentum
!! (to be used with velocity toroidal potential)
!!
!!  Lx = \int r (u_th cos(phi) - u_phi cos(theta)sin(phi) ) zeta^n  dV
!!     = 16pi/3 \int r^3 Tr zeta^n dr    ;       l=1, m=1
!!  Ly = \int r (-u_th sin(phi) - u_phi cos(theta)cos(phi) ) zeta^n  dV
!!     = -16pi/3 \int r^3 Ti zeta^n dr   ;       l=1, m=1
!!  Lx = \int r u_phi sin(theta) zeta^n  dV = 8pi/3 \int r^3 Tr zeta^n dr
!!      ;l=1, m=0
!!___
!------------------------------------------------------------------------
   subroutine var_coll_mom_comp(mes_oc, tor, Lx, Ly, Lz)
      use parameters, only : d_Pi
      use dyn_mpi
      type(rdom), intent(in) :: mes_oc !< outer core domain
      type (coll), intent(in) :: tor !< toroidal part of field
      double precision, intent(out) :: Lx !< int r (u_th cos(phi) - u_phi cos(theta)sin(phi) ) zeta^n  dV
      double precision, intent(out) :: Ly !< int r (-u_th sin(phi) - u_phi cos(theta)cos(phi) ) zeta^n  dV
      double precision, intent(out) :: Lz !< int r u_phi sin(theta) zeta^n  dV = 8pi/3 \int r^3 Tr zeta^n dr
      double precision :: w, f(tor%D%N)
      integer :: N_,n
      _loop_lm_vars
      Lx=0.0
      Ly=0.0
      Lz=0.0
      N_  = tor%D%N
      _loop_lm_begin(tor)
         if ((l == 1) .and. (m_ <= 1)) then
                w = 8d0/3d0 * d_PI
                if (m_ == 0) then ! Lz(l=1,m=0) = 8*pi/3 * int(r^3*zeta^n*Re(T(r)),dr)
                        f =  tor%Re(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Re(T(r))
                        Lz  = w*dot_product(f,tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f,dr)
                else ! Lx(l=1,m=1) =  16*pi/3 * int(r^3*zeta^n*Re(T(r)),dr)
                         ! Ly(l=1,m=1) = -16*pi/3 * int(r^3*zeta^n*Im(T(r)),dr)
                        f =  tor%Re(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Re(T(r))
                    Lx  = 2d0*w*dot_product(f,tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f,dr)
                    f =  tor%Im(1:N_,nh)*mes_oc%r(1:N_,1) ! r*Im(T(r))
                    Ly  = -2d0*w*dot_product(f,tor%D%intr2dr_comp(1:N_)) ! int(r^2*zeta^n*f,dr)
                end if
         end if
      _loop_lm_end

#ifdef _MPI
      call mpi_allreduce( Lx,  w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Lx  = w
      call mpi_allreduce( Ly, w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Ly = w
      call mpi_allreduce( Lz,   w, 1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Lz   = w
#endif

   end subroutine var_coll_mom_comp

!------------------------------------------------------------------------
!>  norm  \int a zeta^(n+1) dV
!------------------------------------------------------------------------
   subroutine var_coll_normS(a, E)
      use parameters, only : d_Pi
      use dyn_mpi, only : mpi_rnk
      type (coll), intent(in) :: a !< scalar field, spectral space
      double precision, intent(out) :: E !< int a zeta^(n+1) dV
      integer :: N_,n,j


      if (mpi_rnk==0) then

       N_  = a%D%N
       E = 0d0
       do j=1,N_
         E = E + 4d0 * d_PI * a%Re(j,0) * a%D%intr2dr_compS(j)
       enddo

      end if

   end subroutine var_coll_normS


!------------------------------------------------------------------------
!>  norm   \int a.a dV  for each l,m
!!
!! Used to produce energy spectra over l, m.
!------------------------------------------------------------------------
   subroutine var_coll_normlm(a, Eh,El,Em)
      use parameters, only : d_Pi, i_Kl, i_M1, i_Mp
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field, spectral space
      double precision, intent(out) :: Eh(0:2*i_KL)
      double precision, intent(out) :: El(0:i_L1) !< int(a.a)dV for each l
      double precision, intent(out) :: Em(0:i_M1) !< int(a.a)dV for each m
#ifdef _MPI
      double precision :: Eh_(0:2*i_KL)
      double precision :: El_(0:i_L1), Em_(0:i_M1)
#endif
      double precision :: f(a%D%N)
      double precision :: w, e
      integer :: N_,n, h, m
      _loop_lm_vars
      
      Eh = 0d0
      El = 0d0
      Em = 0d0
      N_  = a%D%N
      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         m = m_/i_Mp
         if (m /= 0) w = 4d0*w 
         f =  w * ( a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                  + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         e = dot_product(f,a%D%intr2dr(1:N_))
         El(l) = El(l) + e
         Em(m) = Em(m) + e
         do h = 0, 2*i_KL
            Eh(h) = Eh(h) + dot_product(f,a%D%intr2drh(1:N_,h))
         end do
      _loop_lm_end

#ifdef _MPI
      call mpi_allreduce( Eh, Eh_, 1+2*i_KL, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eh = Eh_
      call mpi_allreduce( El, El_, 1+i_L1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      El = El_
      call mpi_allreduce( Em, Em_, 1+i_M1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em = Em_
#endif
   
   end subroutine var_coll_normlm

!------------------------------------------------------------------------
!>  compressible norm   \int a.a rho dV  for each l,m
!------------------------------------------------------------------------

   subroutine var_coll_normlm_comp(a, Eh,El,Em)
      use parameters, only : d_Pi, i_Kl, i_M1, i_Mp
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field, spectral space
      double precision, intent(out) :: Eh(0:2*i_KL)
      double precision, intent(out) :: El(0:i_L1) !< int(a.a rho)dV for each l
      double precision, intent(out) :: Em(0:i_M1) !< int(a.a rho)dV for each m
#ifdef _MPI
      double precision :: Eh_(0:2*i_KL)
      double precision :: El_(0:i_L1), Em_(0:i_M1)
#endif
      double precision :: f(a%D%N)
      double precision :: w, e
      integer :: N_,n, h, m
      _loop_lm_vars
      Eh = 0d0
      El = 0d0
      Em = 0d0
      N_  = a%D%N
      _loop_lm_begin(a)
         w = 4d0 * d_PI / dble(2*l+1)
         m = m_/i_Mp
         if(m /=0 ) w = 4d0*w
         f =  w * ( a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                  + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         e = dot_product(f,a%D%intr2dr_comp(1:N_))
         El(l) = El(l) + e
         Em(m) = Em(m) + e
         do h = 0, 2*i_KL
            Eh(h) = Eh(h) + dot_product(f,a%D%intr2drh(1:N_,h))
         end do
      _loop_lm_end

#ifdef _MPI
      call mpi_allreduce( Eh, Eh_, 1+2*i_KL, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Eh = Eh_
      call mpi_allreduce( El, El_, 1+i_L1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      El = El_
      call mpi_allreduce( Em, Em_, 1+i_M1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      Em = Em_
#endif

   end subroutine var_coll_normlm_comp

!------------------------------------------------------------------------
!>  norm  0.5  \int ubar_phi^2 dV where ubar_phi is the zonal average of
!!  u_phi. 
!------------------------------------------------------------------------
   subroutine var_coll_norm_zon(a, E0)
      use parameters, only : d_PI
      use dyn_mpi
      type (coll), intent(in) :: a !< scalar field, spectral space
      double precision, intent(out) :: E0 !< Total E, int ubar_phi^2 dv
      double precision :: El(0:i_L1) !< total E for each l
#ifdef _MPI
      double precision :: El_(0:i_L1) !< local (mpi) E for each l
#endif
      double precision :: f(a%D%N)
      double precision :: w, e
      integer :: N_,n, h
      _loop_lm_vars

      El = 0d0
      N_  = a%D%N
      _loop_lm_begin(a)
         w = 2d0 * d_PI * dble(l)*dble(l+1) / dble(2*l+1)
         if(m_==0) then
         f =  w *  (a%Re(1:N_,nh)*a%Re(1:N_,nh)  &
                  + a%Im(1:N_,nh)*a%Im(1:N_,nh) )
         e = dot_product(f,a%D%intr2dr_comp(1:N_))
         El(l) = El(l) + e
         end if
      _loop_lm_end

#ifdef _MPI
      call mpi_allreduce( El, El_, 1+i_L1, mpi_double_precision,  &
         mpi_sum, mpi_comm_world, mpi_er)
      El = El_
#endif
      E0=0d0
      do l=0,i_L1
        E0 = El(l)+E0
      end do

   end subroutine var_coll_norm_zon



!*************************************************************************
 end module variables
!*************************************************************************

