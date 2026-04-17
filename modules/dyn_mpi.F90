#include "../parallel.h"
!*************************************************************************
!> MPI fortran initialisation.
!!
!! Most of MPI code is contained in `parallel.h`
!!*************************************************************************
 module dyn_mpi
!*************************************************************************
#ifdef _MPI
   use mpi
#endif
   implicit none
   save
   
#ifdef _MPI
   integer :: mpi_er, mpi_tg, mpi_rq(0:2*_Np), mpi_st(mpi_status_size)
   integer :: th_comm !< mpi communicator for theta
   integer :: r_comm !< mpi communicator for radius
#endif
   integer :: mpi_rnk !< global mpi rank
   integer :: mpi_sze !< global mpi size
   integer :: th_rnk !< theta mpi rank
   integer :: r_rnk !< radial mpi rank
   integer :: th_size !< theta mpi size
   integer :: r_size !< radial mpi size
   
 contains

!------------------------------------------------------------------------
!>  initialise the mpi variables.
!------------------------------------------------------------------------
   subroutine mpi_precompute()
      mpi_rnk = 0
      mpi_sze = 1
      th_rnk = 0
      th_size = 1
      r_rnk = 0
      r_size = 1
#ifdef _MPI
      call mpi_init(mpi_er)
      call mpi_comm_rank(mpi_comm_world, mpi_rnk, mpi_er)
      call mpi_comm_size(mpi_comm_world, mpi_sze, mpi_er)
      if(mpi_rnk == 0) then
         print *, '_Np =', _Np
         print *, '_Nth =', _Nth
         print *, '_Nr =', _Nr
         print *, 'mpi_sze =', mpi_sze  
      end if
      
      if(mpi_sze /= _Np) stop 'mpi_precompute: incorrect num procs'
      call split_mpi()
#endif
   end subroutine mpi_precompute

   !> Split mpi into r and th groups
   !!
   !! 
   subroutine split_mpi()
#ifdef _MPI
      call mpi_comm_split(mpi_comm_world, modulo(mpi_rnk,_Nr), mpi_rnk, th_comm, mpi_er)
      call mpi_comm_split(mpi_comm_world, mpi_rnk/_Nr, mpi_rnk, r_comm, mpi_er)
      call mpi_comm_rank(th_comm,th_rnk, mpi_er)
      call mpi_comm_rank(r_comm,r_rnk, mpi_er)
      call mpi_comm_size(th_comm, th_size, mpi_er)
      call mpi_comm_size(r_comm,r_size, mpi_er)
      write(*,*) "mpi_rnk:",mpi_rnk, ", r_rnk:", r_rnk,", modulo(mpi_rnk,_Nr):",modulo(mpi_rnk,_Nr), ", th_rnk:", th_rnk,", mpi_rnk/_Nr",mpi_rnk/_Nr, ", r_size:", r_size,", th_size:",th_size
#endif
   end subroutine split_mpi



!-------------------------------------------------------------------------
!> Get system time
!!
!! Used to get system time for d_start and d_stop
!! for accurate calculation of simulation duration/performance.
!! Use MPI_WTIME for mpi runs, fall back to system_clock where
!! in serial.
!!-------------------------------------------------------------------------
#ifndef _MPI
   subroutine get_time(d_current_time)
      double precision, intent(out) :: d_current_time !< system time
      integer :: count, count_max
      real :: count_rate
      call system_clock(count, count_rate, count_max)
      d_current_time = dble(count/count_rate)
   end subroutine get_time
#else
   subroutine get_time(d_current_time)
      double precision, intent(out) :: d_current_time !< system time
      d_current_time =  MPI_WTIME()
   end subroutine get_time
#endif
 
   
!*************************************************************************
 end module dyn_mpi
!*************************************************************************
