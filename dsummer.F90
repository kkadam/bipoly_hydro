!*******************************************************************************
!*
!*  DSUMMER
!*
!*******************************************************************************
subroutine dsummer(temp,x,l1flag,xl1,sum1,sum2)
implicit none
include 'runhydro.h'
!*******************************************************************************
!*
!  dsummer sums up the array temp and assigns a given element in the sum to
!  star1 or star2 (sum1 or sum2) depending on where that element lies relative
!  to the L1 point and does not sum over elements that lie beyoned each 
!  star's Roche lobe
!*
!*******************************************************************************
!*
!*  Subroutine Arguments

real, dimension(numr_dd,numz_dd,numphi) :: temp

real, dimension(numr_dd,numphi) :: x

real :: xl1, sum1, sum2

integer :: l1flag

!*
!*******************************************************************************
!*
!*  Global Variables

!*
!*******************************************************************************
!*
!*  Local Variables

integer :: J, K, L

!*
!*******************************************************************************
sum1 = 0.0
sum2 = 0.0

if(l1flag <= 1) then
! L1 on *1 side
   do L = 1, numphi_by_four
      do K = zlwb, zupb
         do J = rlwb, rupb
            if( x(j,l) >= xl1 ) then
               sum1 = sum1 + temp(J,K,L)
            else
               sum2 = sum2 + temp(J,K,L)
            endif
         enddo
      enddo
   enddo 
   do L = 3*numphi_by_four+1, numphi
      do K = zlwb, zupb
         do J = rlwb, rupb
            if( x(J,L) >= xl1 ) then
               sum1 = sum1 + temp(J,K,L)
            else
               sum2 = sum2 + temp(J,K,L)
            endif 
         enddo
      enddo
   enddo 
   do L = numphi_by_four+1, 3*numphi_by_four
      do K = zlwb, zupb
         do J = rlwb, rupb
            sum2 = sum2 + temp(J,K,L)
         enddo
      enddo
   enddo 
else
! l1flag = 2 case, L1 point is on the star 2 side
   do L = 1, numphi_by_four
      do K = zlwb, zupb
         do J = rlwb, rupb
            sum1 = sum1 + temp(J,K,L)
         enddo
      enddo
   enddo 
   do L = 3*numphi_by_four+1, numphi
      do K = zlwb, zupb
         do J = rlwb, rupb
            sum1 = sum1 + temp(J,K,L)
         enddo
      enddo 
   enddo 
   do L = numphi_by_four+1, 3*numphi_by_four
      do K = zlwb, zupb
         do J = rlwb, rupb
            if( x(J,L) <= xl1 ) then
               sum2 = sum2 + temp(J,K,L)
            else
               sum1 = sum1 + temp(J,K,L)
            endif
         enddo 
      enddo 
   enddo
endif

return 
end subroutine dsummer
