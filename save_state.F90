!*************************************************************
!*
!*  SAVE_STATE
!*
!*************************************************************
subroutine save_state
implicit none
include 'runhydro.h'
!*************************************************************
!*
!  save copies s, t, a, and rho into temp arrays so they
!  are not changed by the first call to advect in a 
!  timestep cycle.  The first call to advect in a cycle
!  is only used to get time centered velocities to do a 
!  full timestep of advection.
!* 
!*************************************************************
!*
!*   Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: s1, t1, a1, rho1
common /save/ s1, t1, a1, rho1

!*
!*************************************************************
s1 = s
t1 = t
a1 = a
rho1 = rho
     
return
end subroutine save_state
