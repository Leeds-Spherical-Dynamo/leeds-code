#include "../parallel.h"
!*************************************************************************
!> Transformation using the FFT, and Gauss Pts.
!!
!! De-aliased.\n
!! Written for FFTW version 3 -- documentation at  www.fftw.org\n
!! See FFTW reference, sections 2.3, 4.3.3, 4.4.1-2
!!
!! \f$ \sum_{-Ma}^{Ma} == \sum_0^{2Ma-1}\f$\n
!! logical size of DFT n==2Ma=Ph, \n
!! fftw complex array has n/2+1 elements => Ma+1 elements, or index 0:Ma \n
!! fftw real array 2*(n/2+1) elts, only first n accessed => index 0:Ph-1\n
!!
!! fft at each theta => howmany=Th\n
!! theta along first index of X,Y => dist=1, stride=Th
!!*************************************************************************
 module transform
!*************************************************************************
   use parameters, only : i_ph, i_Ma, i_M1, i_pTh, i_Th, i_pN, i_H, i_M
   use iso_c_binding
   use iso_fortran_env, only: real64
   implicit none
#ifdef _SHTNS
   include 'shtns.f03'
#endif
   save
      					! see fftw3.f
   integer, parameter, private :: FFTW_ESTIMATE=64
   integer, parameter, private :: FFTW_PATIENT =32
   double precision,   private :: a2d_X(i_pTh,0:i_Ph-1) !< fft temp var:  PHYS, real
   double complex,     private :: a2c_Y(i_pTh,0:i_Ma) !< fft temp var: SPEC, complex
   double precision,   private :: a3d_X8(i_pTh,0:i_Ph-1,  i_pN,8)
   double complex,     private :: a3c_Y8(i_Th, 0:2*_pHc-1,i_pN,8)
   double complex,     private :: a3c_T8(i_pTh,0:i_M1,    i_pN,8)
   double precision,   private :: a3d_X3(i_pTh,0:i_Ph-1,  i_pN,8) !< fft temp var:  PHYS, real
   double complex,     private :: a3c_Y3(i_Th, 0:2*_pHc-1,i_pN,3) !< fft temp var:  SPEC, cpmplex
   double complex,     private :: a3c_T3(i_pTh,0:i_M1,    i_pN,3)
   integer*8,          private :: plan_s2r, plan_r2s

   integer, parameter,        private :: dpi=real64

#ifdef _SHTNS
   type(shtns_info), pointer, private :: shtns
   type(c_ptr),               private :: shtns_c
   complex(dpi), allocatable,  private :: shtns_Slm(:), shtns_tra(:)
   complex(dpi), allocatable, private :: shtns_qst_q(:), shtns_qst_s(:), shtns_qst_t(:)
   complex(dpi), allocatable, private :: shtns_tra_r(:), shtns_tra_th(:), shtns_tra_ph(:)
   double precision, parameter :: one_div_sqrt_two = 1.0_dp/sqrt(2.0_dp)
   double precision, private :: shtns_sqrt_ll1(0:i_L) !< normalisation of qst harms for m=0
   double precision, private :: shtns_sqrt_2_ll1(0:i_L) !< normalisation of qst harms for m>0
#endif   
   integer :: m_start(0:2*_pHc) 
   integer :: m_array(0:2*_pHc-1)
   integer :: m_pos_array(0:2*_pHc-1)
   integer :: length(0:2*_pHc-1)

 contains


!------------------------------------------------------------------------
!> Setup plans for multiple transforms.
!!  -  priority = 0   Quick plan, slower ffts -- for only a few ffts
!!  -  priority = 1   Slow plan, quicker ffts -- for many ffts.
!------------------------------------------------------------------------
   subroutine tra_precompute(priority)

      use variables, only : var_th

      integer, intent(in) :: priority !< 0: quick fft plan, 1: slow fft plan

      integer :: flag, n(1), howmany, inembed(1), onembed(1)
      
      if(i_M>1 .and. modulo(i_M,2)==1) stop 'if(M>1) set M even'
      
      if(priority.eq.0) flag = FFTW_ESTIMATE
      if(priority.eq.1) flag = FFTW_PATIENT

      n = (/i_Ph/)
      howmany = var_Th%pTh
      inembed = (/i_pTh*(i_Ma+1)/)
      onembed = (/i_pTh*i_Ph/)
      call dfftw_plan_many_dft_c2r(plan_s2r, 1, n, howmany,  &
                                   a2c_Y, inembed, i_pTh, 1,  &
                                   a2d_X, onembed, i_pTh, 1,  &
                                   flag)
      n = (/i_Ph/)
      howmany = var_Th%pTh
      inembed = (/i_pTh*i_Ph/)
      onembed = (/i_pTh*(i_Ma+1)/)
      call dfftw_plan_many_dft_r2c(plan_r2s, 1, n, howmany,  &
                                   a2d_X, inembed, i_pTh, 1,  &
                                   a2c_Y, onembed, i_pTh, 1,  &
                                   flag)

#ifdef _SHTNS
      call tra_shtns_pre()
#endif
   end subroutine tra_precompute


#ifdef _SHTNS
   subroutine tra_shtns_pre()
      integer :: norm, layout, nthreads
      real(dp) :: eps_polar = 0.0d0
      integer :: m_before, m_after ,i, l, ll
      double precision :: ldbl
      
      _loop_sumth_vars

      ! allocate shtns memory
      !---------------------------
      
      norm = SHT_SCHMIDT + SHT_NO_CS_PHASE 
      layout = SHT_GAUSS + SHT_THETA_CONTIGUOUS + SHT_SOUTH_POLE_FIRST

      if (mpi_rnk==0) then 
         call shtns_verbose(1)
      else 
         call shtns_verbose(0)
      end if
      nthreads = shtns_use_threads(0)
      shtns_c=shtns_create(i_L-1,i_M-1,i_Mp,norm)
      call shtns_set_grid(shtns_c, layout, eps_polar, i_Th, i_Ph)
      call c_f_pointer(cptr=shtns_c, fptr=shtns)
      allocate( shtns_Slm(0:i_L))
      allocate(shtns_tra(1:i_Th) )
      !!!! Allocate the other shtns arrays here!!!!
      allocate(shtns_qst_q(0:i_L))
      allocate(shtns_qst_s(0:i_L))
      allocate(shtns_qst_t(0:i_L))
      allocate(shtns_tra_r(1:i_Th))
      allocate(shtns_tra_th(1:i_Th))
      allocate(shtns_tra_ph(1:i_Th))

      shtns_slm=0.0_dpi
      shtns_tra=0.0_dpi
      shtns_qst_q=0.0_dpi
      shtns_qst_s=0.0_dpi
      shtns_qst_t=0.0_dpi
      shtns_tra_r=0.0_dpi
      shtns_tra_th=0.0_dpi
      shtns_tra_ph=0.0_dpi

      ! set up normalisation array for leeds qst to shtns qst
      do ll=0, i_L-1
         ldbl=dble(ll)
         shtns_sqrt_ll1(ll) = sqrt(1.0_dp/(ldbl*(ldbl+1))) 
         shtns_sqrt_2_ll1(ll) = sqrt(2.0_dp/(ldbl*(ldbl+1))) 
      end do

      ! get start and end of each l.
      m_before=0
      i=0
      _loop_sumth_begin
         m_after = m

         if ( (nh==0) .or. (m_before /= m_after)) then 
            m_start(i) = nh 
            m_array(i) = m_ 
            m_pos_array(i) = m
            i = i+1
            m_before=m_after
         end if
      _loop_sumth_end
      do i = 0, 2*_pHc-1
         length(i)=m_start(i+1) - m_start(i)
      end do
      m_start(2*_pHc) = i_H1
      length(2*_pHc-1) = i_H1-m_start(2*_pHc-1)
   end subroutine tra_shtns_pre
#endif

!------------------------------------------------------------------------
!>  transpose theta_j and m-modes then fft
!!
!! Input data should be prepared in a3c_Y8 before calling tra_s2p_vel.
!! Output data should be found in a3d_X8.
!! For each radial point:
!!    a3c_Y8 is put into a2c_Y. FFT then transforms, outputting
!!    to a2d_X. a2d_X is then put into relevant part of a3d_X8.  
!------------------------------------------------------------------------
   subroutine tra_s2p_vel(s,ns)
      use variables, only : spec, var_th
      use dyn_mpi
      type (spec), intent(in) :: s !< used to determine size of data
      integer,     intent(in) :: ns
      integer :: n, m, mp, i
#ifndef _MPI
      do i = 1, ns
         do n = 1, s%D%pN
            do m = 0, var_Th%pM-1
               mp = var_Th%m(m,0)
               a2c_Y(:,mp) = a3c_Y8(:,m,n,i)
            end do
            a2c_Y(:,i_M:) = 0d0
            call dfftw_execute(plan_s2r)
            a3d_X8(:,:,n,i) = a2d_X

         end do
      end do
#else
      double precision  :: bsend8(2*i_pN*i_pTh*(2*_pHc)*8,0:_Nth-1)
      double precision  :: brecv8(2*i_pN*i_pTh*(2*_pHc)*8,0:_Nth-1)
      integer :: stp, dst,src, l,j
         					! transpose
      
      do stp = 0, _Nth-1
         src  = modulo(_Nth-stp+th_rnk, _Nth)         
         mpi_tg =  stp*_Nr
         n = 2*s%D%pN*var_Th%pTh*var_Th%pM_(src)*ns
         call mpi_irecv( brecv8(1,stp), n, mpi_double_precision,  &
            src, mpi_tg, th_comm, mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Nth-1
         dst  = modulo(stp+th_rnk, _Nth)
         src  = modulo(_Nth-stp+th_rnk, _Nth)         
         l = 1
         do i = 1, ns
           do n = 1, s%D%pN
             do m = 0, var_Th%pM-1
               do j = var_Th%pThi_(dst), var_Th%pThi_(dst)+var_Th%pTh_(dst)-1
                  bsend8(l,  stp) =  dble(a3c_Y8(j,m,n,i))
                  bsend8(l+1,stp) = dimag(a3c_Y8(j,m,n,i))
                  l = l + 2
               end do
             end do
           end do
         end do
         mpi_tg = stp*_Nr
         n = 2*s%D%pN*var_Th%pTh_(dst)*var_Th%pM*ns
         call mpi_isend( bsend8(1,stp), n, mpi_double_precision,  &
            dst, mpi_tg, th_comm, mpi_rq(_Nth+stp), mpi_er)
      end do
      
      do stp = 0, _Nth-1
         src  = modulo(_Nth-stp+th_rnk, _Nth)      
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do i = 1, ns
           do n = 1, s%D%pN
             do m = 0, var_Th%pM_(src)-1
               mp = var_Th%m(m,src)
               do j = 1, var_Th%pTh
                  a3c_T8(j,mp,n,i) = dcmplx(brecv8(l,stp),brecv8(l+1,stp))
                  l = l + 2
               end do
             end do
           end do
         end do
      end do

      do stp = 0, _Nth-1
         call mpi_wait( mpi_rq(_Nth+stp), mpi_st, mpi_er)
      end do
         					! FFT
      do i = 1, ns
         do n = 1, s%D%pN
            a2c_Y(:,:i_M1) = a3c_T8(:,:,n,i)
            a2c_Y(:,i_M:)  = 0d0
            call dfftw_execute(plan_s2r)
            a3d_X8(:,:,n,i) = a2d_X
         end do
      end do
#endif      

   end subroutine tra_s2p_vel



!------------------------------------------------------------------------
!>  transpose theta_j and m-modes then fft
!!
!! Input data should be prepared in a3c_Y3 before calling tra_s2p.
!! Output data should be found in a3d_X3.
!! For each radial point:
!!    a3c_Y3 is put into a2c_Y. FFT then transforms, outputting
!!    to a2d_X. a2d_X is then put into relevant part of a3d_X3.  
!------------------------------------------------------------------------
   subroutine tra_s2p(s,ns)
      use variables, only : spec, var_th
      use dyn_mpi
      type (spec), intent(in) :: s !< used to determine size of data
      integer,     intent(in) :: ns !< number of fields to transform
      integer :: n, m, mp, i
#ifndef _MPI
      do i = 1, ns
         do n = 1, s%D%pN
            do m = 0, var_Th%pM-1
               mp = var_Th%m(m,0)
               a2c_Y(:,mp) = a3c_Y3(:,m,n,i)
            end do
            a2c_Y(:,i_M:) = 0d0
            call dfftw_execute(plan_s2r)
            a3d_X3(:,:,n,i) = a2d_X
         end do
      end do
#else
      double precision :: bsend(2*i_pN*i_pTh*(2*_pHc)*3,0:_Nth-1)
      double precision :: brecv(2*i_pN*i_pTh*(2*_pHc)*3,0:_Nth-1)
      integer :: stp, dst,src, l,j
         					! transpose
      
      do stp = 0, _Nth-1
         src  = modulo(_Nth-stp+th_rnk, _Nth)         
         mpi_tg = stp*_Nr
         n = 2*s%D%pN*var_Th%pTh*var_Th%pM_(src)*ns
         call mpi_irecv( brecv(1,stp), n, mpi_double_precision,  &
            src, mpi_tg, th_comm, mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Nth-1
         dst  = modulo(stp+th_rnk, _Nth)
         src  = modulo(_Nth-stp+th_rnk, _Nth)         
         l = 1
         do i = 1, ns
           do n = 1, s%D%pN
             do m = 0, var_Th%pM-1
               do j = var_Th%pThi_(dst), var_Th%pThi_(dst)+var_Th%pTh_(dst)-1
                  bsend(l,  stp) =  dreal(a3c_Y3(j,m,n,i))
                  bsend(l+1,stp) = dimag(a3c_Y3(j,m,n,i))
                  l = l + 2
               end do
             end do
           end do
         end do
         mpi_tg = stp*_Nr
         n = 2*s%D%pN*var_Th%pTh_(dst)*var_Th%pM*ns
         call mpi_isend( bsend(1,stp), n, mpi_double_precision,  &
            dst, mpi_tg, th_comm, mpi_rq(_Nth+stp), mpi_er)
      end do
      
      do stp = 0, _Nth-1
         src  = modulo(_Nth-stp+th_rnk, _Nth)      
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do i = 1, ns
           do n = 1, s%D%pN
             do m = 0, var_Th%pM_(src)-1
               mp = var_Th%m(m,src)
               do j = 1, var_Th%pTh
                  a3c_T3(j,mp,n,i) = dcmplx(brecv(l,stp),brecv(l+1,stp))
                  l = l + 2
               end do
             end do
           end do
         end do
      end do

      do stp = 0, _Nth-1
         call mpi_wait( mpi_rq(_Nth+stp), mpi_st, mpi_er)
      end do
         					! FFT
      do i = 1, ns
         do n = 1, s%D%pN
            a2c_Y(:,:i_M1) = a3c_T3(:,:,n,i)
            a2c_Y(:,i_M:)  = 0d0
            call dfftw_execute(plan_s2r)
            a3d_X3(:,:,n,i) = a2d_X
         end do
      end do
#endif      

   end subroutine tra_s2p



!------------------------------------------------------------------------
!>  fft then transpose theta_j and m-modes
!!
!! Data should be prepared in a3d_X3 before calling tra_p2s.
!! Output data found in a3c_Y3.
!! For each radial gridpoint:
!!    data is moved from
!!    a3d_X3 to a2d_X. FFTW then transforms,
!!    taking input from a2d_X and outputting in a2c_Y.
!!    a2c_Y output then stored in a3c_Y3
!------------------------------------------------------------------------
   subroutine tra_p2s(s,ns)
      use variables, only : spec, var_th
      use dyn_mpi
      type (spec), intent(in) :: s !< used to determine size of data
      integer,     intent(in) :: ns !< number of fields to transform
      double precision :: scale
      integer :: n, m, mp, i
         			! scale, FFTW section 4.7.2
#ifndef _MPI
      scale = 1d0 / dble(i_Ph)
      do i = 1, ns
         do n = 1, s%D%pN
            a2d_X = a3d_X3(:,:,n,i) * scale
            call dfftw_execute(plan_r2s)

            do m = 0, var_Th%pM-1
               mp = var_Th%m(m,0)
               a3c_Y3(:,m,n,i) = a2c_Y(:,mp)            
            end do
         end do
      end do
#else
      double precision :: bsend(2*i_pN*i_pTh*(2*_pHc)*3,0:_Nth-1)
      double precision :: brecv(2*i_pN*i_pTh*(2*_pHc)*3,0:_Nth-1)
      integer :: stp, dst,src, l,j
         					! FFT
      scale = 1d0 / dble(i_Ph)
      do i = 1, ns
         do n = 1, s%D%pN
            a2d_X = a3d_X3(:,:,n,i) * scale
            call dfftw_execute(plan_r2s)

            a3c_T3(:,:,n,i) = a2c_Y(:,:i_M1)
         end do
      end do
         					! transpose
      
      do stp = 0, _Nth-1
         src  = modulo(_Nth-stp+th_rnk, _Nth)         
         mpi_tg = stp*_Nr
         n = 2*s%D%pN*var_Th%pTh_(src)*var_Th%pM*ns
         call mpi_irecv( brecv(1,stp), n, mpi_double_precision,  &
            src, mpi_tg, th_comm, mpi_rq(stp), mpi_er)
      end do

      do stp = 0, _Nth-1
         dst  = modulo(stp+th_rnk, _Nth)
         l = 1
         do i = 1, ns
           do n = 1, s%D%pN
             do m = 0, var_Th%pM_(dst)-1
               mp = var_Th%m(m,dst)
               do j = 1, var_Th%pTh
                  bsend(l,  stp) =  dble(a3c_T3(j,mp,n,i))
                  bsend(l+1,stp) = dimag(a3c_T3(j,mp,n,i))
                  l = l + 2
               end do
             end do
           end do
         end do
         mpi_tg = stp*_Nr
         n = 2*s%D%pN*var_Th%pTh*var_Th%pM_(dst)*ns
         call mpi_isend( bsend(1,stp), n, mpi_double_precision,  &
            dst, mpi_tg, th_comm, mpi_rq(_Nth+stp), mpi_er)
      end do
      
      do stp = 0, _Nth-1
         src  = modulo(_Nth-stp+th_rnk, _Nth)      
         call mpi_wait( mpi_rq(stp), mpi_st, mpi_er)
         l = 1
         do i = 1, ns
           do n = 1, s%D%pN
             do m = 0, var_Th%pM-1
               do j = var_Th%pThi_(src), var_Th%pThi_(src)+var_Th%pTh_(src)-1
                  a3c_Y3(j,m,n,i) = dcmplx(brecv(l,stp),brecv(l+1,stp))
                  l = l + 2
               end do
             end do
           end do
         end do
      end do

      do stp = 0, _Nth-1
         call mpi_wait( mpi_rq(_Nth+stp), mpi_st, mpi_er)
      end do
#endif
         
   end subroutine tra_p2s


!------------------------------------------------------------------------
!>  Convert spectral to real space
!------------------------------------------------------------------------
   subroutine tra_spec2phys(s, p)

      use variables, only : spec, phys, var_H
      use parameters, only : i_H1
      use legendre, only : leg_P

      type (spec), intent(in)  :: s !< spectral data to transform
      type (phys), intent(out) :: p !< transformed data in physical space

      integer, parameter :: Th2 = (i_Th+1)/2, Th_ = i_Th/2
      double complex :: c, Y(Th2, 0:2*_pHc-1, 0:1)
      integer :: n, i, im
      integer, save :: counter=0
      integer, save :: counter1=0
      _loop_sumth_vars      
      				! for each x_n ...
      ! legendre transform
#ifndef _SHTNS
         do n = 1, s%D%pN
            			! sum over l at each th_j -> Y_jm
            Y = 0d0
            _loop_sumth_begin
               c = dcmplx( s%Re(nh,n), s%Im(nh,n) )
               Y(:,m,s__) = Y(:,m,s__) + c * leg_P(:,nh)
            _loop_sumth_end
            a3c_Y3(  :Th2,:,n,1) = Y(:,:,0) + Y(:,:,1) 
            a3c_Y3(Th2+1:,:,n,1) = Y(Th_:1:-1,:,0) - Y(Th_:1:-1,:,1)
         end do
#else
         do n = 1, s%D%pN
            do i = 0, 2*_pHc-1
               m=m_pos_array(i)
               m_=m_array(i)
               im=m_/i_Mp ! shtns im defined by m=im*mres, so im = m/mres
               shtns_slm(0:length(i)) = cmplx( s%Re(m_start(m):m_start(m+1),n), s%Im(m_start(m):m_start(m+1),n) )
               call SH_to_spat_ml(shtns_c, im, shtns_slm, shtns_tra(1:i_Th), i_L-1)
               if (m_>0) then
                  a3c_Y3(:,m,n,1) = shtns_tra(1:i_Th)*sqrt(2.0_dp)
               else
                  a3c_Y3(:,m,n,1) = shtns_tra(1:i_Th)
               end if
            end do
         end do
#endif
            			! fft
      call tra_s2p(s,1)
      p%Re = a3d_X3(:,:,:,1)
   
   end subroutine tra_spec2phys


!------------------------------------------------------------------------
!>  Convert real to spectral space
!!  Ensure output s%L and s%M are set before call.
!------------------------------------------------------------------------
   subroutine tra_phys2spec(p, s)
      use variables, only : phys, spec

      type (phys), intent (in  )  :: p !< physcal scalar field to transform
      type (spec), intent (inout) :: s !< transformed data in spectral space

      integer, parameter :: Th2 = (i_Th+1)/2, Th_ = i_Th/2+1
      double complex :: Y(Th2, 0:2*_pHc-1, 0:1)
      integer :: n, i, im
      integer, save :: counter=0
      _loop_sumth_vars
        ! fft
      a3d_X3(:,:,:,1) = p%Re
      call tra_p2s(s,1)

      ! Choose between native legendre transform and shnts
#ifndef _SHTNS
         do n = 1, s%D%pN
            do m = 0, 2*_pHc-1
               Y(:,m,0) = a3c_Y3(:Th2,m,n,1) + a3c_Y3(i_Th:Th_:-1,m,n,1)
               Y(:,m,1) = a3c_Y3(:Th2,m,n,1) - a3c_Y3(i_Th:Th_:-1,m,n,1)
               if(Th2==Th_)  Y(Th2,m,0) = Y(Th2,m,0) * 0.5d0
            end do
            ! integrate over theta, Gauss weights
            call tra_leg_trans_p2s_native(Y, s%re(:,n), s%im(:,n))
         end do
#else
         do n=1, s%D%pN
            ! want to loop over all m's here, and set the correct 
            ! indices for s%Re and s%Im for outputs to go in. 
            ! m is the local m index.
            ! m_ is the actual order of spherical harmonic
            do i = 0, 2*_pHc-1
               shtns_slm=0.0_dpi
               m=m_pos_array(i)
               m_=m_array(i)
               im=m_/i_Mp ! shtns im defined by m=im*mres, so im = m/mres
               !shtns_tra(:) = a3c_Y3(i_th:1:-1,m,n,1)
               !shtns_tra(:) = a3c_Y3(:,m,n,1)
               call spat_to_SH_ml(shtns_c, im, a3c_Y3(:,m,n,1), shtns_slm, i_L-1)
               ! normalise all coefficients by 1/sqrt(2) where m != 0. 
               ! `m` in the code is a local index for m, rather than the acutal
               ! spherical harmonic order, so can't check explicitly for m==0 on
               ! all processors. So if process is not in a group that contains
               ! m=0, just normalise all coefficients!
               !if ( (m>0) .or. ( mpi_rnk/_Nr/=0/_Nr) ) then
               if (m_>0) then
                  s%Re(m_start(m):m_start(m+1),n)=dble(shtns_slm(0:length(i))) * one_div_sqrt_two!(1.0_dpi/sqrt(2.0_dpi))
                  s%Im(m_start(m):m_start(m+1),n)=dimag(shtns_slm(0:length(i))) * one_div_sqrt_two!(1.0_dpi/sqrt(2.0_dpi))
               else
                  s%Re(m_start(m):m_start(m+1),n)=dble(shtns_slm(0:length(i)))
                  s%Im(m_start(m):m_start(m+1),n)=dimag(shtns_slm(0:length(i)))
               end if
            end do
         end do
#endif
   end subroutine tra_phys2spec

!> Legendre transform, physical to spectral
   subroutine tra_leg_trans_p2s_native(in, out_re, out_im) ! using symmetry
      use parameters, only : i_H1
      use variables, only : var_H
      use legendre, only : leg_gwP

      double complex :: in(:,0:,0:) !< fft'd data ready for legendre transform
      double precision :: out_re(0:) !< real part of legendre transformed data
      double precision :: out_im(0:) !< imaginary part of legendre transformed data
      
      _loop_sumth_vars
      _loop_sumth_begin   
         out_re(nh) = sum( dble(in(:,m,s__)) * leg_gwP(:,nh))    
         out_im(nh) = sum(dimag(in(:,m,s__)) * leg_gwP(:,nh))
         if (m_==0) out_im(nh) = 0d0
      _loop_sumth_end

   end subroutine tra_leg_trans_p2s_native



!------------------------------------------------------------------------
!>  Spectral to real space, qst -> r,th,ph
!!
!! Used for anelastic version of code.
!------------------------------------------------------------------------
   subroutine tra_qst2rtp_vel(mes_oc, q,s,t, qr, sr, tr, r,th,ph,dur,qv1,dutr,dupr,qv2)
      use meshs, only : rdom
      use variables, only : var_H, spec, phys
      use parameters, only : i_H1
      use legendre, only : leg_P, leg_llP_, leg_llsin1mP, leg_qva, leg_qvb

      type(rdom), intent(in) :: mes_oc !< outer core domain
      type (spec), intent(in)  :: q !< q component of qst spherical harmonic decomposition
      type (spec), intent(in)  :: s !< s component of qst spherical harmonic decomposition
      type (spec), intent(in)  :: t !< t component of qst spherical harmonic decomposition
      type (spec), intent(in)  :: qr,sr,tr
      type (phys), intent(out) :: r !< radial component of physical field
      type (phys), intent(out) :: th !< theta component of physical field
      type (phys), intent(out) :: ph !< phi component of physical field
      type (phys), intent(out) :: dur,qv1,dutr,dupr,qv2
      
      integer, parameter  :: Th2 = (i_Th+1)/2, Th_ = i_Th/2 
      double complex :: c1,c2,c3,c4,c5,c6,c7,c8
      double complex ::  Y1(Th2,0:2*_pHc-1, 0:1), Y2(Th2,0:2*_pHc-1, 0:1) 
      double complex ::  Y3(Th2,0:2*_pHc-1, 0:1), Y4(Th2,0:2*_pHc-1, 0:1)
      double precision :: rinv_
      integer :: n
      _loop_sumth_vars      
         			! Ar = q
      do n = 1, s%D%pN
         Y1 = 0d0
         Y2 = 0d0
         _loop_sumth_begin
            c1 = dcmplx( q%Re(nh,n), q%Im(nh,n) )
            c2 = dcmplx( qr%Re(nh,n), qr%Im(nh,n) )
            Y1(:,m,s__) = Y1(:,m,s__) + c1 * leg_P(:,nh)
            Y2(:,m,s__) = Y2(:,m,s__) + c2 * leg_P(:,nh)
         _loop_sumth_end  
            a3c_Y8(  :Th2,:,n,1) = Y1(:,:,0) + Y1(:,:,1)
            a3c_Y8(Th2+1:,:,n,1) = Y1(Th_:1:-1,:,0) - Y1(Th_:1:-1,:,1)
            a3c_Y8(  :Th2,:,n,4) = Y2(:,:,0) + Y2(:,:,1)
            a3c_Y8(Th2+1:,:,n,4) = Y2(Th_:1:-1,:,0) - Y2(Th_:1:-1,:,1)
      end do     
                                ! At = s dth - t/sin dph
      do n = 1, s%D%pN
         Y1 = 0d0
         Y2 = 0d0
         _loop_sumth_begin
            c1 = dcmplx( s%Re(nh,n),  s%Im(nh,n) )
            c2 = dcmplx( t%Im(nh,n), -t%Re(nh,n) )
            c3 = dcmplx( sr%Re(nh,n),  sr%Im(nh,n) )
            c4 = dcmplx( tr%Im(nh,n), -tr%Re(nh,n) ) 
            Y1(:,m,1-s__) = Y1(:,m,1-s__) + c1 * leg_llP_(:,nh)  
            Y1(:,m,s__) = Y1(:,m,s__)     + c2 * leg_llsin1mP(:,nh)
            Y2(:,m,1-s__) = Y2(:,m,1-s__) + c3 * leg_llP_(:,nh) 
            Y2(:,m,s__) = Y2(:,m,s__)     + c4 * leg_llsin1mP(:,nh)
          _loop_sumth_end 

           a3c_Y8(  :Th2,:,n,2) = Y1(:,:,0) + Y1(:,:,1)
           a3c_Y8(Th2+1:,:,n,2) = Y1(Th_:1:-1,:,0) - Y1(Th_:1:-1,:,1)
           a3c_Y8(  :Th2,:,n,6) = Y2(:,:,0) + Y2(:,:,1)
           a3c_Y8(Th2+1:,:,n,6) = Y2(Th_:1:-1,:,0) - Y2(Th_:1:-1,:,1)           
      end do
                                ! Aph = s/sin dph + t dth
      do n = 1, s%D%pN
         Y1 = 0d0
         Y2 = 0d0
         Y3 = 0d0
         Y4 = 0d0
         rinv_=mes_oc%r(n+mes_oc%pNi-1,-1)   
         _loop_sumth_begin
            c1 = dcmplx( -s%Im(nh,n), s%Re(nh,n) )
            c2 = dcmplx(  t%Re(nh,n), t%Im(nh,n) )
            c3 = dcmplx( -sr%Im(nh,n), sr%Re(nh,n) )
            c4 = dcmplx(  tr%Re(nh,n), tr%Im(nh,n) )
            c5 = dcmplx( -s%Re(nh,n)*rinv_, -s%Im(nh,n)*rinv_ )
            c6 = dcmplx(  t%Im(nh,n)*rinv_, -t%Re(nh,n)*rinv_ )
            c7 = dcmplx(  t%Re(nh,n)*rinv_, t%Im(nh,n)*rinv_ )
            c8 = dcmplx(  s%Im(nh,n)*rinv_, -s%Re(nh,n)*rinv_ )
            Y1(:,m,s__)   = Y1(:,m,s__)   + c1 * leg_llsin1mP(:,nh)
            Y1(:,m,1-s__) = Y1(:,m,1-s__) + c2 * leg_llP_(:,nh)
            Y2(:,m,s__)   = Y2(:,m,s__)   + c3 * leg_llsin1mP(:,nh)
            Y2(:,m,1-s__) = Y2(:,m,1-s__) + c4 * leg_llP_(:,nh)
            Y3(:,m,s__)   = Y3(:,m,s__)   + c5 * leg_qva(:,nh)
            Y3(:,m,1-s__) = Y3(:,m,1-s__) + c6 * leg_qvb(:,nh)
            Y4(:,m,s__)   = Y4(:,m,s__)   + c7 * leg_qva(:,nh)
            Y4(:,m,1-s__) = Y4(:,m,1-s__) + c8 * leg_qvb(:,nh)
         _loop_sumth_end      
         a3c_Y8(  :Th2,:,n,3) = Y1(:,:,0) + Y1(:,:,1)
         a3c_Y8(Th2+1:,:,n,3) = Y1(Th_:1:-1,:,0) - Y1(Th_:1:-1,:,1)
         a3c_Y8(  :Th2,:,n,7) = Y2(:,:,0) + Y2(:,:,1)
         a3c_Y8(Th2+1:,:,n,7) = Y2(Th_:1:-1,:,0) - Y2(Th_:1:-1,:,1)
         a3c_Y8(  :Th2,:,n,5) = Y3(:,:,0) + Y3(:,:,1)
         a3c_Y8(Th2+1:,:,n,5) = Y3(Th_:1:-1,:,0) - Y3(Th_:1:-1,:,1)
         a3c_Y8(  :Th2,:,n,8) = Y4(:,:,0) + Y4(:,:,1)
         a3c_Y8(Th2+1:,:,n,8) = Y4(Th_:1:-1,:,0) - Y4(Th_:1:-1,:,1) 
      end do

      call tra_s2p_vel(s,8)
      r%Re  = a3d_X8(:,:,:,1)
      th%Re = a3d_X8(:,:,:,2)
      ph%Re = a3d_X8(:,:,:,3)
      dur%Re  = a3d_X8(:,:,:,4)
      qv1%Re  = a3d_X8(:,:,:,5)
      dutr%Re = a3d_X8(:,:,:,6)
      dupr%Re = a3d_X8(:,:,:,7)
      qv2%Re  = a3d_X8(:,:,:,8)

   end subroutine tra_qst2rtp_vel

!------------------------------------------------------------------------
!>  Spectral to real space, qst -> r,th,ph
!------------------------------------------------------------------------
   subroutine tra_qst2rtp(q,s,t, r,th,ph)
      use variables, only : spec, phys, var_H
      use parameters, only : i_H1, i_L
      use legendre, only : leg_P, leg_llP_, leg_llsin1mP

      type (spec), intent(in)  :: q !< q component of qst spherical harmonic decomposition
      type (spec), intent(in)  :: s !< s component of qst spherical harmonic decomposition
      type (spec), intent(in)  :: t !< t component of qst spherical harmonic decomposition
      type (phys), intent(out) :: r !< radial component of physical field
      type (phys), intent(out) :: th !< theta component of physical field
      type (phys), intent(out) :: ph !< phi component of physical field
      
      integer, parameter  :: Th2 = (i_Th+1)/2, Th_ = i_Th/2
      double complex ::  Y1(Th2,0:2*_pHc-1, 0:1), Y2(Th2,0:2*_pHc-1, 0:1)
      double complex ::  Y3(Th2,0:2*_pHc-1, 0:1)
      double precision :: l_array(0:i_L-1)
      double complex :: c1,c2,c3,c4
      integer :: n, i, im
      double precision :: q_fac
      integer, save :: counter = 0
      _loop_sumth_vars      
         			

      ! native legendre
      !---------------------------------------------------------------
#ifndef _SHTNS
                                 ! Ar = q
         do n = 1, q%D%pN
            Y1 = 0d0
            _loop_sumth_begin
               c1 = dcmplx( q%Re(nh,n), q%Im(nh,n) )
               Y1(:,m,s__) = Y1(:,m,s__) + c1 * leg_P(:,nh)
            _loop_sumth_end
            a3c_Y3(  :Th2,:,n,1) = Y1(:,:,0) +  Y1(:,:,1)
            a3c_Y3(Th2+1:,:,n,1) = Y1(Th_:1:-1,:,0) - Y1(Th_:1:-1,:,1)
         end do
                                ! Ath = s dth - t/sin dph
                                ! Aph = s/sin dph + t dth
         do n = 1, s%D%pN
            Y2 = 0d0
            Y3 = 0d0
            _loop_sumth_begin
               c1 = dcmplx( s%Re(nh,n),  s%Im(nh,n) )
               c2 = dcmplx( t%Im(nh,n), -t%Re(nh,n) )
               c3 = dcmplx( -s%Im(nh,n), s%Re(nh,n) )
               c4 = dcmplx(  t%Re(nh,n), t%Im(nh,n) )
               Y2(:,m,1-s__) = Y2(:,m,1-s__) + c1 * leg_llP_(:,nh)  
               Y2(:,m,s__)   = Y2(:,m,s__)   + c2 * leg_llsin1mP(:,nh)
               Y3(:,m,s__)   = Y3(:,m,s__)   + c3 * leg_llsin1mP(:,nh)
               Y3(:,m,1-s__) = Y3(:,m,1-s__) + c4 * leg_llP_(:,nh)
            _loop_sumth_end
            a3c_Y3(  :Th2,:,n,2) = Y2(:,:,0) +  Y2(:,:,1)
            a3c_Y3(Th2+1:,:,n,2) = Y2(Th_:1:-1,:,0) - Y2(Th_:1:-1,:,1)
            a3c_Y3(  :Th2,:,n,3) = Y3(:,:,0) +  Y3(:,:,1)
            a3c_Y3(Th2+1:,:,n,3) = Y3(Th_:1:-1,:,0) - Y3(Th_:1:-1,:,1)
         end do
#else
      ! Use shtns
         do n = 1, q%D%pN
            do i = 0, 2*_pHc-1
               m=m_pos_array(i) ! m is local m index.
               m_=m_array(i) ! m_ is the actual value of m
               im=m_/i_Mp ! shtns im defined by m=im*mres, so im = m/mres

               ! Normalisation for leeds qst definition vs shtns qst definition
               ! factor of sqrt(2) for r (m>0), sqrt(2/(l(l+1))) for s (m>0), -sqrt(2/(l(l+1))) for t (m>0), sqrt(1/(l(l+1))) for s (m=0), and -sqrt(1/(l(l+1))) for t(m=0).
               ! Need to multiply by these factors before calling Shqst_to_spat-ml.
               q_fac=1.0_dp
               if (m_>0) then
                  q_fac=sqrt(2.0_dp)
               end if
               
               ! normalise for shtns and put in array ready for shtns legendre transform
               if (m_>0) then 
                  shtns_qst_q(0:length(i)) = cmplx( q_fac * q%Re(m_start(m):m_start(m+1),n), q_fac * q%Im(m_start(m):m_start(m+1),n) )
                  shtns_qst_s(0:length(i)) = cmplx( s%Re(m_start(m):m_start(m+1),n) * shtns_sqrt_2_ll1(m_:) , s%Im(m_start(m):m_start(m+1),n) * shtns_sqrt_2_ll1(m_:) )
                  shtns_qst_t(0:length(i)) = cmplx( -1*t%Re(m_start(m):m_start(m+1),n) *shtns_sqrt_2_ll1(m_:) , -1*t%Im(m_start(m):m_start(m+1),n) * shtns_sqrt_2_ll1(m_:) )
               else
                  shtns_qst_q(0:length(i)) = cmplx( q_fac * q%Re(m_start(m):m_start(m+1),n), q_fac * q%Im(m_start(m):m_start(m+1),n) )
                  shtns_qst_s(0:length(i)) = cmplx( s%Re(m_start(m):m_start(m+1),n) * shtns_sqrt_ll1(m_:) , s%Im(m_start(m):m_start(m+1),n) * shtns_sqrt_ll1(m_:) )
                  shtns_qst_t(0:length(i)) = cmplx( -1*t%Re(m_start(m):m_start(m+1),n) * shtns_sqrt_ll1(m_:) , -1*t%Im(m_start(m):m_start(m+1),n) *  shtns_sqrt_ll1(m_:) )
               end if
               ! carry out shtns legendre transform
               call SHqst_to_spat_ml(shtns_c, im, shtns_qst_q(0:i_L), shtns_qst_s(0:i_L), shtns_qst_t(0:i_L), shtns_tra_r(1:i_Th), shtns_tra_th(1:i_Th), shtns_tra_ph(1:i_Th), i_L-1)
               ! copy output to fftw arrays ready for fft and parallel transpose
               a3c_Y3(:,m,n,1)  = shtns_tra_r(1:i_Th)
               a3c_Y3(:,m,n,2)  = shtns_tra_th(1:i_Th)
               a3c_Y3(:,m,n,3)  = shtns_tra_ph(1:i_Th)
            end do
         end do
#endif

      ! fftw transform and parallel transpose
      call tra_s2p(s,3)
      ! copy fftw arrays to physical variables
      r%Re  = a3d_X3(:,:,:,1)
      th%Re = a3d_X3(:,:,:,2)
      ph%Re = a3d_X3(:,:,:,3)

   end subroutine tra_qst2rtp


!------------------------------------------------------------------------
!>  Real to spectral space, r,th,ph -> q,s,t
!!  Ensure output q,s,t%L and q,s,t%M are set before call.
!------------------------------------------------------------------------
   subroutine tra_rtp2qst(r,th,ph, q,s,t)
      use variables, only : phys, spec, var_H
      use parameters, only : i_H1
      use legendre, only : leg_gwllP_, leg_gwllsin1mP, leg_gwP

      type (phys), intent(in) :: r !< radial component of physical field
      type (phys), intent(in) :: th !< theta component of physical field
      type (phys), intent(in) :: ph !< phi component of physical field
      type (spec), intent(out)  :: q !< q component of qst spherical harmonic decomposition
      type (spec), intent(out)  :: s !< s component of qst spherical harmonic decomposition
      type (spec), intent(out)  :: t !< t component of qst spherical harmonic decomposition
    
      integer, parameter :: Th2 = (i_Th+1)/2, Th_ = i_Th/2+1
      double complex :: Y1(Th2, 0:2*_pHc-1, 0:1), Y2(Th2, 0:2*_pHc-1, 0:1)
      double complex :: Y3(Th2, 0:2*_pHc-1, 0:1)
      integer :: n, i, im
      integer, save :: counter=0
      _loop_sumth_vars

      a3d_X3(:,:,:,1) = r%Re
      a3d_X3(:,:,:,2) = th%Re
      a3d_X3(:,:,:,3) = ph%Re
      call tra_p2s(s,3)      
      
      ! q  integrate q.r = q
#ifndef _SHTNS
      !-----------------------------------------------------------------
      ! native transform

      ! leg_P(j,nh) = legendre polynomial for l(nh),m(nh) where j is index for all theta   pa( (l*l+1)/2+m_+1 )
      ! leg_gw = legendre gauss weights
      ! leg_th = theta points
      ! leg_gwP = leg_P * leg_gw * (2l+1)/4  
      ! so , leg_gwP is Plm * w * (2l+1)/4
      !-----------------------------------------------------------------
         do n = 1, s%D%pN
            do m = 0, 2*_pHc-1
               Y1(:,m,0) = a3c_Y3(:Th2,m,n,1) + a3c_Y3(i_Th:Th_:-1,m,n,1)
               Y1(:,m,1) = a3c_Y3(:Th2,m,n,1) - a3c_Y3(i_Th:Th_:-1,m,n,1)
               if(Th2==Th_)  Y1(Th2,m,0) = Y1(Th2,m,0) * 0.5d0
            end do
            _loop_sumth_begin
               q%Re(nh,n) = sum( dble(Y1(:,m,s__))*leg_gwP(:,nh))
               q%Im(nh,n) = sum(dimag(Y1(:,m,s__))*leg_gwP(:,nh))
               if(m_==0) s%Im(nh,n) = 0d0
            _loop_sumth_end
         end do

      ! Rest of transforms
      do n = 1, s%D%pN
         do m = 0, 2*_pHc-1
            Y2(:,m,0) = a3c_Y3(:Th2,m,n,2) + a3c_Y3(i_Th:Th_:-1,m,n,2)
            Y2(:,m,1) = a3c_Y3(:Th2,m,n,2) - a3c_Y3(i_Th:Th_:-1,m,n,2)
            if(Th2==Th_)  Y2(Th2,m,0) = Y2(Th2,m,0) * 0.5d0
         end do


         do m = 0, 2*_pHc-1
            Y3(:,m,0) = a3c_Y3(:Th2,m,n,3) + a3c_Y3(i_Th:Th_:-1,m,n,3)
            Y3(:,m,1) = a3c_Y3(:Th2,m,n,3) - a3c_Y3(i_Th:Th_:-1,m,n,3)
            if(Th2==Th_)  Y3(Th2,m,0) = Y3(Th2,m,0) * 0.5d0
         end do
         _loop_sumth_begin
         			! s.th = dth P

            ! leg_gwllp_ : leg_gw(j) * (2l+1)/4  / sqrt(l(l+1)) * leg_P_(j,nh)

            s%Re(nh,n) =  sum( dble(Y2(:,m,1-s__))*leg_gwllP_(:,nh))
            s%Im(nh,n) =  sum(dimag(Y2(:,m,1-s__))*leg_gwllP_(:,nh))
            			! t.th = - 1/sin(th) dph P
            t%Re(nh,n) = -sum(dimag(Y2(:,m,s__))*leg_gwllsin1mP(:,nh))
            t%Im(nh,n) =  sum( dble(Y2(:,m,s__))*leg_gwllsin1mP(:,nh))
         			! s.ph = 1/sin(th) dph P
            s%Re(nh,n) =  &
               s%Re(nh,n) + sum(dimag(Y3(:,m,s__))*leg_gwllsin1mP(:,nh))
            s%Im(nh,n) =  &
               s%Im(nh,n) - sum( dble(Y3(:,m,s__))*leg_gwllsin1mP(:,nh))
            if(m_==0) s%Im(nh,n) = 0d0
            			! t.ph = dth P
            t%Re(nh,n) =  &
               t%Re(nh,n) + sum( dble(Y3(:,m,1-s__))*leg_gwllP_(:,nh))
            t%Im(nh,n) =  &
               t%Im(nh,n) + sum(dimag(Y3(:,m,1-s__))*leg_gwllP_(:,nh))
            if(m_==0) t%Im(nh,n) = 0d0
         _loop_sumth_end
      end do
      !-----------------------------------------------------------------
      ! Shtns
      !-----------------------------------------------------------------
#else
         do n=1, s%D%pN
            do i = 0,  2*_pHc-1
               shtns_qst_q=0.0_dpi
               shtns_qst_s=0.0_dpi
               shtns_qst_t=0.0_dpi
               m=m_pos_array(i)
               m_=m_array(i)
               im=m_/i_Mp ! shtns im defined by m=im*mres, so im = m/mres, mres = i_Mp

               call spat_to_SHqst_ml(shtns_c, im, &
                                 & a3c_Y3(:,m,n,1), a3c_Y3(:,m,n,2), a3c_Y3(:,m,n,3),& !r, th, ph, in
                                 & shtns_qst_q, shtns_qst_s, shtns_qst_t,&  ! q,s,t, out
                                 & i_L-1 )
               if (m_>0) then 
                  q%Re(m_start(m):m_start(m+1),n)=dble(shtns_qst_q(0:length(i))) * one_div_sqrt_two!
                  q%Im(m_start(m):m_start(m+1),n)=dimag(shtns_qst_q(0:length(i))) * one_div_sqrt_two
                  s%Re(m_start(m):m_start(m+1),n)=dble(shtns_qst_s(0:length(i)))  / shtns_sqrt_2_ll1(m_:m_+length(i))
                  s%Im(m_start(m):m_start(m+1),n)=dimag(shtns_qst_s(0:length(i))) / shtns_sqrt_2_ll1(m_:m_+length(i))
                  t%Re(m_start(m):m_start(m+1),n)=dble(shtns_qst_t(0:length(i))) *(-1.0d0) / shtns_sqrt_2_ll1(m_:m_+length(i))
                  t%Im(m_start(m):m_start(m+1),n)=dimag(shtns_qst_t(0:length(i))) *(-1.0d0) / shtns_sqrt_2_ll1(m_:m_+length(i))
               else 
                  q%Re(m_start(m):m_start(m+1),n)=dble(shtns_qst_q(0:length(i))) 
                  q%Im(m_start(m):m_start(m+1),n)=dimag(shtns_qst_q(0:length(i))) 
                  s%Re(m_start(m):m_start(m+1),n)=dble(shtns_qst_s(0:length(i))) / shtns_sqrt_ll1(m_:m_+length(i))
                  s%Im(m_start(m):m_start(m+1),n)=dimag(shtns_qst_s(0:length(i))) / shtns_sqrt_ll1(m_:m_+length(i))
                  t%Re(m_start(m):m_start(m+1),n)=dble(shtns_qst_t(0:length(i))) *(-1.0d0)/ shtns_sqrt_ll1(m_:m_+length(i))
                  t%Im(m_start(m):m_start(m+1),n)=dimag(shtns_qst_t(0:length(i))) *(-1.0d0)/ shtns_sqrt_ll1(m_:m_+length(i))
               end if
            end do
         end do
#endif

      

   end subroutine tra_rtp2qst


!------------------------------------------------------------------------
!>  Spectral to real space, A -> (grad A)t,p,  th,ph components only
!------------------------------------------------------------------------
   subroutine tra_grad(s, t,p)
      use variables, only : spec, phys, var_H
      use parameters, only : i_H1
      use legendre, only : leg_P_, leg_sin1mp

      type (spec), intent(in)  :: s !< s component of qst spherical harmonic decomposition
      type (phys), intent(out) :: t !< theta component of physical field
      type (phys), intent(out) :: p !< phi component of physical field
      
      integer, parameter :: Th2 = (i_Th+1)/2, Th_ = i_Th/2
      double complex :: c, Y(Th2, 0:2*_pHc-1, 0:1)
      double precision :: r1(s%D%pN)
      integer :: n
      _loop_sumth_vars      
      
      do n = 1, s%D%pN
         r1(n) = s%D%r(n+s%D%pNi-1,-1)
      end do
                                ! Ath = 1/r dth a
      do n = 1, s%D%pN
         Y = 0d0
         _loop_sumth_begin
            c = dcmplx( s%Re(nh,n), s%Im(nh,n) )
            c = c * r1(n)
            Y(:,m,1-s__) = Y(:,m,1-s__) + c * leg_P_(:,nh)
         _loop_sumth_end
         a3c_Y3(  :Th2,:,n,1) = Y(:,:,0) + Y(:,:,1) 
         a3c_Y3(Th2+1:,:,n,1) = Y(Th_:1:-1,:,0) - Y(Th_:1:-1,:,1)
      end do
                                ! Aph = 1/rsin dph a
      do n = 1, s%D%pN
         Y = 0d0
         _loop_sumth_begin
            c = dcmplx( -s%Im(nh,n), s%Re(nh,n) )
            c = c * r1(n)
            Y(:,m,s__) = Y(:,m,s__)  +  c * leg_sin1mP(:,nh)
         _loop_sumth_end
         a3c_Y3(  :Th2,:,n,2) = Y(:,:,0) + Y(:,:,1) 
         a3c_Y3(Th2+1:,:,n,2) = Y(Th_:1:-1,:,0) - Y(Th_:1:-1,:,1)
      end do

      call tra_s2p(s,2)
      t%Re = a3d_X3(:,:,:,1)
      p%Re = a3d_X3(:,:,:,2)     

   end subroutine tra_grad


!------------------------------------------------------------------------
!>  theta derivative
!------------------------------------------------------------------------
   subroutine tra_dth(s, p)
      use variables, only : spec, phys, var_H
      use parameters, only : i_H1
      use legendre, only : leg_P_

      type (spec), intent(in)  :: s !< spectral field to differentiate
      type (phys), intent(out) :: p !< ds/dtheta, in physical space

      integer, parameter :: Th2 = (i_Th+1)/2, Th_ = i_Th/2
      double complex :: c, Y(Th2, 0:2*_pHc-1, 0:1)
      integer :: n
      _loop_sumth_vars      
      
      do n = 1, s%D%pN
         Y = 0d0
         _loop_sumth_begin
            c = dcmplx( s%Re(nh,n), s%Im(nh,n) )
            Y(:,m,1-s__) = Y(:,m,1-s__) + c * leg_P_(:,nh)
         _loop_sumth_end
         a3c_Y3(  :Th2,:,n,1) = Y(:,:,0) + Y(:,:,1) 
         a3c_Y3(Th2+1:,:,n,1) = Y(Th_:1:-1,:,0) - Y(Th_:1:-1,:,1)
      end do
      call tra_s2p(s,1)
      p%Re = a3d_X3(:,:,:,1)     
   
   end subroutine tra_dth


!------------------------------------------------------------------------
!> Surface integral:
!!
!! \f$ \int_0^{2\pi} \int_0^\pi A \sin(\theta) d\theta d\phi = 4 \pi A00 \f$
!------------------------------------------------------------------------
   subroutine tra_surfintegral(p, S)
      use variables, only : var_th
      use parameters, only : d_Pi
      use dyn_mpi
      use legendre, only : leg_gwP0
      double precision, intent(in)  :: p(i_pTh, 0:i_Ph-1) !< field over a spherical surface (parallelised over theta)
      double precision, intent(out) :: S !< integral of field over spherical surface
      double precision :: Y(i_pTh), S_
      integer :: m, j, r

      j = var_Th%pTh
      Y = 0d0
      do m = 0, i_Ph-1
         Y(:j) = Y(:j) + p(:j,m)
      end do
      S_ = dot_product( Y(:j) , leg_gwP0(var_Th%pThi:var_Th%pThi+j-1) )
      S_ = S_ * 4d0*d_PI / dble(i_Ph)
#ifdef _MPI
      call mpi_reduce(S_, S, 1, mpi_double_precision, mpi_sum, &
          0, th_comm, mpi_er)
#else
      S=S_
#endif
   end subroutine tra_surfintegral


!*************************************************************************
 end module transform
!*************************************************************************
