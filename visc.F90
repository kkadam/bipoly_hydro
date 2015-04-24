!************************************************************************
!*  
!*  VISC
!*
!************************************************************************
subroutine visc
implicit none
include 'runhydro.h'
!************************************************************************
!*
!
!*
!************************************************************************
!*
!*  Subroutine Arguments
!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep


!*
!************************************************************************
!*
!*   Local Variables

real, dimension(numr_dd,numz_dd,numphi) :: qrr, qzz, qpp

integer :: J, K, L

!*
!************************************************************************
!  initialize the local variables

! fill in the viscous stress tensor
call stress(qrr,qzz,qpp)

! update radial momentum density, S
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         s(J,K,L) = s(J,K,L) - dt*drinv*(qrr(J,K,L) - qrr(J-1,K,L))
      enddo
   enddo
enddo

! update the vertical momentum density, T
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         t(J,K,L) = t(J,K,L) - dt*dzinv*(qzz(J,K,L) - qzz(J,K-1,L))
      enddo
   enddo
enddo

! update the angular momentum density, A
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         a(J,K,L) = a(J,K,L) - dt*dphiinv*(qpp(J,K,L) - qpp(J,K,L-1))
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      a(J,K,philwb) = a(J,K,philwb) - dt*dphiinv*                      &
                      (qpp(J,K,philwb) - qpp(J,K,phiupb))
   enddo
enddo

call source_bc

!  share these updated values at processor edges with neighbors
call comm(s)
call comm(t)
call comm(a)

return
end subroutine visc
