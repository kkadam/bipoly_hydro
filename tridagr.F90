!***************************************************************************
!*
!*  TRIDAGR
!*
!***************************************************************************
subroutine tridagr(ar,br,cr,knownr,potr)
implicit none
include 'runhydro.h'
!***************************************************************************
!*
!  tridagr solves for potr from the linear tridiagonal system of equations
!  with br being the diagonal elements of the matrix, ar and cr are the
!  off-diagonal elements and knownr is the right hand side.  The code
!  comes from section 2.4, page 43 of Numerical Recipes in Fortran, 2nd ed.
!*
!***************************************************************************
!*
!*  Subroutine Arguments

real, dimension(numr) :: ar, cr

real, dimension(numr,numphi_dd) :: br

real, dimension(numr,numz_dd,numphi_dd) :: knownr, potr

!*
!***************************************************************************
!*
!*  Local Variables

real, dimension(numr,numphi) :: bet, gam  

integer :: J, K, L

!*
!***************************************************************************
!  initialize the local variables
gam = 0.0
bet = 0.0

! setup
do L = 1, numphi_dd
   bet(2,L) = br(2,L)
enddo
do L = 1, numphi_dd
   do K = zlwb, zupb
      potr(2,K,L) = knownr(2,K,L) / bet(2,L)
   enddo
enddo

!  decomposition and forward substitution
do L = 1, numphi_dd
   do J = 3, numr-1
      gam(J,L) = cr(J-1) / bet(J-1,L)
      bet(J,L) = br(J,L) - ar(J)*gam(J,L)
      do K = zlwb, zupb
         potr(J,K,L) = (knownr(J,K,L)-ar(J)*potr(J-1,K,L)) / bet(J,L)
      enddo
   enddo
enddo

! back subsitution
do L = 1, numphi_dd
   do K = zlwb, zupb
      do J = numr-2, 2, -1
         potr(J,K,L) = potr(J,K,L) - gam(J+1,L)*potr(J+1,K,L)
      enddo
   enddo
enddo

return
end subroutine tridagr
