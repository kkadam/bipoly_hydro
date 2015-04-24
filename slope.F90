!************************************************************************
!*
!*  SLOPE
!*
!************************************************************************
subroutine slope( q, qj, qk, ql )
implicit none
include 'runhydro.h'
!************************************************************************
!*
!  slope calculates the linear slope for each conserved density
!  that is used in the advection subroutine to interpolate the
!  conserved densities to the cell interface
!
!  for each conserved density calculate its slope in each coordinate
!  direction and reset the slope to zero if the cell represents a local
!  extreema or if the slope is within the parameter limit of being
!  zero
!
!*
!************************************************************************
!*
!*  Subroutine Arguments

real, dimension(numr_dd,numz_dd,numphi) :: q, qj, qk, ql

!*
!************************************************************************
!*
!*  Global Variables

!*
!***************************************************************************
!* 
!*  Local Variables

real :: temp

real, parameter :: limit = 1.0e-20

integer :: J, K, L

!*
!****************************************************************************
!  initialize the local variables
temp = 0.0
 
!  slopes in radial direction
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         qj(J,K,L) = (q(J+1,K,L)-q(J,K,L))*(q(J,K,L)-q(J-1,K,L))
         temp = q(J+1,K,L) - q(J-1,K,L)
         if( qj(J,K,L) < 0.0 .or. abs(temp) < limit ) then
            qj(J,K,L) = 0.0
         else
            qj(J,K,L) = 2.0 * qj(J,K,L) / temp
         endif
      enddo
   enddo
enddo

!  slopes in vertical direction
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         qk(J,K,L) = (q(J,K+1,L) - q(J,K,L))*(q(J,K,L) - q(J,K-1,L))      
         temp = q(J,K+1,L) - q(J,K-1,L)
         if( qk(J,K,L) < 0.0 .or. abs(temp) < limit ) then
            qk(J,K,L) = 0.0
         else
            qk(J,K,L) = 2.0 * qk(J,K,L) / temp
         endif
      enddo
   enddo
enddo

!  slopes in azimuthal direction
do L = philwb+1,phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         ql(J,K,L) = (q(J,K,L+1) - q(J,K,L))*(q(J,K,L) - q(J,K,L-1))
         temp = q(J,K,L+1) - q(J,K,L-1)
         if( ql(J,K,L) < 0.0 .or. abs(temp) < limit ) then
            ql(J,K,L) = 0.0
         else
            ql(J,K,L) = 2.0 * ql(J,K,L) / temp
        endif
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      ql(J,K,1) = (q(J,K,2) - q(J,K,1))*(q(J,K,1) - q(J,K,numphi))
      temp = q(J,K,2) - q(J,K,numphi)
      if( ql(J,K,1) < 0.0 .or. abs(temp) < limit ) then
         ql(J,K,1) = 0.0
      else
         ql(J,K,1) = 2.0 * ql(J,K,1) / temp
      endif
      ql(J,K,numphi) = (q(J,K,1) - q(J,K,numphi))*(q(J,K,numphi) - q(J,K,numphi-1))
      temp = q(J,K,1) - q(J,K,numphi-1)
      if( ql(J,K,numphi) < 0.0 .or. abs(temp) < limit ) then
         ql(J,K,numphi) = 0.0
      else
         ql(J,K,numphi) = 2.0 * ql(J,K,numphi) / temp
      endif
   enddo
enddo

! need to communicate here too...but only need to pass slopes in 
! J direction in that direction and slopes in K direction in
! the K direction
call comm_dir(qj,1)
call comm_dir(qk,2)

return
end subroutine slope
