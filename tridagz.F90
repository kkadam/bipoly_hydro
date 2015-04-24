!***************************************************************************
!*
!*  TRIDAGZ
!*
!***************************************************************************
subroutine tridagz(az,bz,cz,knownz,potz)
implicit none
include 'runhydro.h'
!***************************************************************************
!*
!  tridagz solves for potz from the linear tridiagonal system of equations
!  with bz being the diagonal elements of the matrix, az and cz are the
!  off-diagonal elements and knownz is the right hand side.  The code
!  comes from section 2.4, page 43 of Numerical Recipes in Fortran, 2nd ed.
!*
!***************************************************************************
!*
!*  Subroutine Arguments

real :: az, cz

real, dimension(numz) :: bz

real, dimension(numr_dd_z,numz,numphi_dd) :: knownz, potz

!*
!***************************************************************************
!*
!*  Local Variables

real, dimension(numz) :: bet, gam  

integer :: J, K, L, rad_lwr_bnd, rad_upr_bnd

!*
!***************************************************************************
!  initialize the local vraiables
bet = 0.0
gam = 0.0
rad_lwr_bnd = 2
rad_upr_bnd = numr_dd_z - 1

! setup
bet(2) = bz(2)
do L = 1, numphi_dd
   do J = rad_lwr_bnd, rad_upr_bnd
      potz(J,2,L) = knownz(J,2,L) / bet(2)
   enddo
enddo
 
!  decomposition and forward substitution
do K = 3, numz-1
   gam(K) = cz / bet(K-1)
   bet(K) = bz(K) - az*gam(K)
   do L = 1, numphi_dd
      do J = rad_lwr_bnd, rad_upr_bnd
         potz(J,K,L) = (knownz(J,K,L) - az*potz(J,K-1,L))/ bet(K)
      enddo
   enddo
enddo

! back subsitution
do L = 1, numphi_dd
   do K = numz-2, 2, -1
      do J = rad_lwr_bnd, rad_upr_bnd
         potz(J,K,L) = potz(J,K,L) - gam(K+1)*potz(J,K+1,L)
      enddo
   enddo
enddo

return
end subroutine tridagz
