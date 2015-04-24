!************************************************************************
!*  
!*  STRESS
!*
!************************************************************************
subroutine stress(qrr,qzz,qpp)
implicit none
include 'runhydro.h'
!************************************************************************
!*
!
!  Stress fills in the diagonal terms of the viscous stress tensor
!
!*
!************************************************************************
!*
!*  Subroutine Arguments

real, dimension(numr_dd,numz_dd,numphi) :: qrr, qzz, qpp

!*
!************************************************************************
!*
!*   Global variables

real,dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn  

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi 

real :: pi, grav, cfl_factor, den_cutoff, vlim_factor, viscosity
common /constants/ pi, grav, cfl_factor, den_cutoff, vlim_factor,    &
                   viscosity

!*
!************************************************************************
!*
!*   Local Variables

real :: dv

integer :: J, K, L

!*
!************************************************************************
!  initialize the local variables
qrr = 0.0
qzz = 0.0
qpp = 0.0

! construct the three diagonal components of the stress
! tensor... ignore the shear terms
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         dv = u(J+1,K,L) - u(J,K,L)
	 if (dv < 0.0 ) qrr(J,K,L) = viscosity*rho(J,K,L)*dv*dv

	 dv = w(J,K+1,L) - w(J,K,L)
	 if( dv < 0.0 ) qzz(J,K,L) = viscosity*rho(J,K,L)*dv*dv

	 dv = rhfinv(J)*(jn(J,K,L+1) - jn(J,K,L))
	 if( dv < 0.0 ) qpp(J,K,L) = viscosity*rho(J,K,L)*dv*dv
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      dv = u(J+1,K,phiupb) - u(J,K,phiupb)
      if( dv < 0.0 ) qrr(J,K,phiupb) = viscosity*rho(J,K,phiupb)*dv*dv

      dv = w(J,K+1,phiupb) - w(J,K,phiupb)
      if( dv < 0.0 ) qzz(J,K,phiupb) = viscosity*rho(J,K,phiupb)*dv*dv

      dv = rhfinv(J)*(jn(J,K,philwb) - jn(J,K,phiupb))
      if( dv < 0.0 ) qpp(J,K,phiupb) = viscosity*rho(J,K,phiupb)*dv*dv
   enddo
enddo


!  share these updated values at processor edges with neighbors
call comm_dir(qrr,1)
call comm_dir(qzz,2)

return
end subroutine stress
