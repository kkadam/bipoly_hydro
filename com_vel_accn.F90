!************************************************************************
!*  
!*  COM_VEL_ACCN
!*
!************************************************************************
subroutine com_vel_accn
implicit none
include 'runhydro.h'
!************************************************************************
!*
! 06/09/2003 : Modified the radial and azimuthal components of the 
!              acceleration of the com. These components are now 
!              calculated in terms of the grid angle phi(L) (defined
!              in setup_frac.F).
!              The x and y components of the acceleration of the com
!              are the same for every cell at a given timestep.
!*
!************************************************************************
!*
!*   Global variables

real,dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd,numz_dd,numphi) :: temps, tempa
common /source_temp/ temps, tempa

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

!*
!************************************************************************
!*
!*   Local variables

integer :: I, J, K, L

#ifdef SHORT 
integer*8 :: ierror
#else
integer :: ierror
#endif

!*
!************************************************************************

! Calculate the cartesian components of the velocity of the com in the
! rotating frame

! v_com is indeterminate if the delt(2) and delt(3) are zero. In this
! case v_com is forced to be zero. This case is only encountered in the
! initial couple of timesteps (tstep = 1 and 2).

do J = 1, 3
!   v_com(J) = v_com(J) + a_com(J)*delt(2)

   if (delt(3) .eq. 0.0) then
       v_com(J) = 0.0
   else
       v_com(J) = com(1,J)*(2*delt(2)*delt(3)+delt(3)**2)           &
                  - com(2,J)*(delt(2)+delt(3))**2                   &
                  + com(3,J)*(delt(2))**2
       v_com(J) = v_com(J)/(delt(2)*delt(3))
       v_com(J) = v_com(J)/(delt(2)+delt(3))
   endif
enddo

! Calculate the corresponding cylindrical components

! 02/20/2002 : Modified calculation for cylindrical components of velocity
!cylin_v_com(1) = v_com(1)*cos(phi_com) + v_com(2)*sin(phi_com)
!cylin_v_com(2) = v_com(2)*cos(phi_com) - v_com(1)*sin(phi_com)
!cylin_v_com(3) = v_com(3)

! 06/17/2003 : Modified
!cylin_v_com(1) = (com(1,1)*v_com(1) + com(1,2)*v_com(2))*R_com_inv
!cylin_v_com(2) = (-com(1,2)*v_com(1) + com(1,1)*v_com(2))*R_com_inv
!cylin_v_com(3) = v_com(3)

do L = 1, numphi
   cylin_v_com(1,L) = v_com(1)*cos_cc(L) + v_com(2)*sin_cc(L)
   cylin_v_com(2,L) = v_com(2)*cos_cc(L) - v_com(1)*sin_cc(L)
   cylin_v_com(3,L) = v_com(3)
enddo


! Calculate the cartesian components of the acceleration of the com in 
! the rotating frame
!   a_com(J) = (com(1,J)-com(1,J))/(delt(2)*delt(2)) +          &
!                ( (com(1,J)-com(3,J))/(delt(2)+delt(3)) )/delt(2)
!   a_com(J) = a_com(J)/1.5

! a_com is indeterminate if the delt(2) and delt(3) are zero. In this
! case a_com is forced to be zero. This case is only encountered in the
! initial couple of timesteps (tstep = 1 and 2).

do J = 1, 3
   if (delt(3) .eq. 0.0) then
       a_com(J) = 0.0
   else
       a_com(J) = 2*( com(1,J)*delt(3) - com(2,J)*(delt(2)+delt(3))  &
                  + com(3,J)*delt(2) )
       a_com(J) = a_com(J)/(delt(2)*delt(3))
       a_com(J) = a_com(J)/(delt(2)+delt(3))
   endif
enddo

! Calculate the corresponding cylindrical components 

! 02/20/2002 : Modified calculation for cylindrical components of acceleration
!cylin_a_com(1) = a_com(1)*cos(phi_com) + a_com(2)*sin(phi_com)
!cylin_a_com(2) = a_com(2)*cos(phi_com) - a_com(1)*sin(phi_com)
!cylin_a_com(3) = a_com(3)

! 06/09/2003 : Modified
!cylin_a_com(1) = (com(1,1)*a_com(1) + com(1,2)*a_com(2))*R_com_inv
!cylin_a_com(2) = (-com(1,2)*a_com(1) + com(1,1)*a_com(2))*R_com_inv
!cylin_a_com(3) = a_com(3)

do L = 1, numphi
   cylin_a_com(1,L) = a_com(1)*cos_cc(L) + a_com(2)*sin_cc(L)
   cylin_a_com(2,L) = a_com(2)*cos_cc(L) - a_com(1)*sin_cc(L)
   cylin_a_com(3,L) = a_com(3)
enddo

return
end subroutine com_vel_accn
