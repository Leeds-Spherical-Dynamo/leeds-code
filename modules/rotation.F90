!*************************************************************************
!> This module deals with the rotation of the inner core.
!! It computes rotation rate and viscous and magnetic torques at each 
!! time step.
!!*************************************************************************
 module rotation
!*************************************************************************

   implicit none
   save

   double precision, private :: velTorq
   double precision, private :: magTorq

 contains

!------------------------------------------------------------------------
!>  Set initial values for rotation rate rot_omega 
!------------------------------------------------------------------------
   subroutine rot_precompute(mes_oc,rot_omega, rot_inertia, rot_velTorq, rot_magTorq)

      use meshs, only : rdom
      use parameters, only : d_Pi

      type(rdom), intent(in) :: mes_oc
      double precision, intent(out) :: rot_omega !< Rotation rate
      double precision, intent(out) :: rot_inertia !< Moment of inertia set to \f$8\pi/15 r_i^5\f$
      double precision, intent(out) :: rot_velTorq !< Viscous torque
      double precision, intent(out) :: rot_magTorq !< Magnetic torque
      rot_omega   = 0d0
      rot_inertia = (8d0*d_PI/15d0) * mes_oc%r(1,1)**5
      rot_velTorq = 0d0
      rot_magTorq = 0d0
   end subroutine rot_precompute



!------------------------------------------------------------------------
!> Predictor step.
!------------------------------------------------------------------------
   subroutine rot_predictor(rot_inertia, rot_velTorq, rot_magTorq, rot_omega)

      use parameters, only : b_rot_force_0
      use timestep, only : tim_dt

      double precision, intent(in) :: rot_inertia !< Moment of inertia set to \f$8\pi/15 r_i^5\f$
      double precision, intent(in) :: rot_velTorq !< Viscous torque
      double precision, intent(in) :: rot_magTorq !< Magnetic torque
      double precision, intent(inout) :: rot_omega !< Rotation rate
      
      if (b_rot_force_0) then 
         rot_omega = 0d0
         return
      else
         rot_omega = rot_omega  &
            +  (tim_dt/rot_inertia) * (rot_velTorq + rot_magTorq) 

         velTorq = rot_velTorq
         magTorq = rot_magTorq
      end if
   end subroutine rot_predictor
   
   
!------------------------------------------------------------------------
!> Corrector step.
!------------------------------------------------------------------------
   subroutine rot_corrector(rot_inertia, rot_velTorq, rot_magTorq, rot_omega)

      use parameters, only : b_rot_force_0, d_implicit
      use timestep, only : tim_dt

      double precision, intent(in) :: rot_inertia !< Moment of inertia set to \f$8\pi/15 r_i^5\f$
      double precision, intent(in) :: rot_velTorq !< Viscous torque
      double precision, intent(in) :: rot_magTorq !< Magnetic torque
      double precision, intent(inout) :: rot_omega !< Rotation rate

      if (b_rot_force_0) then 
         return
      else
         rot_omega = rot_omega  &
            + (tim_dt*d_implicit/rot_inertia)  &
            * ((rot_velTorq + rot_magTorq) - ( velTorq + magTorq ) )

         velTorq = rot_velTorq
         magTorq = rot_magTorq
      end if
   end subroutine rot_corrector
   

   
!*************************************************************************
 end module rotation
!*************************************************************************
