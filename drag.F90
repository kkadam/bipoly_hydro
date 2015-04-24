!********************************************************************
!*
!*  DRAG
!*
!********************************************************************
subroutine drag(drag_factor, timestep)
implicit none
include 'runhydro.h'
include 'mpif.h'
!********************************************************************
!*
!   vel converts the conserved kinematic variables into "velocities"
!   for determining the ammount of material advected
!   for the radial direction:
!        u  =  s / <rho>
!
!   for the vertical direction:
!        w  =  t / <rho>
!
!   for the azimuthal direction:
!        jn  =  a / rho, the specific angular momentum
!
!   vel also limits the velocities by the following algorithm:
!
!   --> if the density in the cell is greater then den_cutoff and the
!        velocity is greater than vmax, reset the velocity to vamx
!   --> if the density falls below den_cutoff and the velocity 
!       exceeds vlim_factor times the maximum sound speed then reset
!       the velocity to the velocity limit.  vlim_factor is given
!       a value in setup
!*
!********************************************************************
!*
!*   Subroutine arguments
!*

real :: drag_factor, timestep

!*
!********************************************************************
!*
!*   Global variables
!*

real,dimension(numr_dd,numz_dd,numphi) :: pot, rho 
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) ::  u, w, jn
common /velocities/ u, w, jn

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

!*
!********************************************************************
!*
!*   Local Variables
!*

real, dimension(numr_dd,numz_dd,numphi) :: rhjkl

integer :: J, K, L

!*
!*
!******************************************************************   
!  initialize the local variables 

! center rho at azimuthal face centers
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl(J,K,L) = 0.5 * (rho(J,K,L) + rho(J,K,L-1))
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      rhjkl(J,K,philwb) = 0.5 * (rho(J,K,philwb) + rho(J,K,phiupb))
   enddo
enddo         

do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         jn(J,K,L) = jn(J,K,L) - omega_frame*rhf(J)*rhf(J)*   &
	             drag_factor*timestep/cirp
         a(J,K,L) = jn(J,K,L)*rhjkl(J,K,L)
      enddo 
   enddo
enddo

! enforce boundary conditions on jn and a
call drag_bc

!   communicate the shared values of jn, a that may have been reset
call comm(jn)
call comm(a)

return
end subroutine drag
