!**************************************************************************
!*
!*  FLUX
!*
!**************************************************************************
subroutine flux(quant, slope, varea, deltat, volinv, dir, retval)
implicit none
include 'runhydro.h'
!**************************************************************************
!*
!*  Subroutine Arguments
      
real quant(numr_dd,numz_dd,numphi)
real slope(numr_dd,numz_dd,numphi)
real varea(numr_dd,numz_dd,numphi)
real retval(numr_dd,numz_dd,numphi)
real volinv(numr_dd)
real deltat
integer dir

!*
!**************************************************************************
!*
!*  Local Variables

real :: qdot
integer ::  J, K, L

!*
!**************************************************************************
!  initialize the local variables
qdot = 0.0
retval = 0.0

!  change in cell's content of quant due to advection in radial direction 
if( dir == 1 ) then
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb-1, rupb
            if( varea(J,K,L) > 0.0 ) then
               qdot = (quant(J,K,L) + 0.5*slope(J,K,L))*             &
                      varea(J,K,L)
            else
               qdot = (quant(J+1,K,L) - 0.5*slope(J+1,K,L))*         &
                      varea(J,K,L)
            endif
            retval(J,K,L) = retval(J,K,L) - deltat*volinv(J)*qdot
            retval(J+1,K,L) = retval(J+1,K,L) + deltat*volinv(J+1)*  &
                              qdot
         enddo
      enddo
   enddo
   retval(rlwb-1,zlwb:zupb,:) = 0.0
   retval(rupb+1,zlwb:zupb,:) = 0.0
endif 

!  same for vertical direction
if( dir == 2 ) then
   do L = philwb, phiupb
      do K = zlwb-1, zupb
         do J = rlwb, rupb
            if( varea(J,K,L) > 0.0 ) then
               qdot = (quant(J,K,L) + 0.5*slope(J,K,L))*             &
                      varea(J,K,L)
            else
               qdot = (quant(J,K+1,L) - 0.5*slope(J,K+1,L))*         &
                      varea(J,K,L)
            endif
            retval(J,K,L) = retval(J,K,L) - deltat*volinv(J)*qdot
            retval(J,K+1,L) = retval(J,K+1,L) + deltat*volinv(J)*    &
                              qdot
         enddo
      enddo
   enddo
   retval(rlwb:rupb,zlwb-1,:) = 0.0
   retval(rlwb:rupb,zupb+1,:) = 0.0
endif
     
!  same for azimuthal direction
if( dir == 3 ) then
   do L = philwb, phiupb-1
      do K = zlwb, zupb
         do J = rlwb, rupb
            if( varea(J,K,L) > 0.0 ) then
               qdot = (quant(J,K,L) + 0.5*slope(J,K,L))*            &
                      varea(J,K,L)
            else
               qdot = (quant(J,K,L+1) - 0.5*slope(J,K,L+1))*        &
                              varea(J,K,L)
            endif
            retval(J,K,L) = retval(J,K,L) - deltat*volinv(J)*qdot
            retval(J,K,L+1) = retval(J,K,L+1) + deltat*volinv(J)*   &
                              qdot
         enddo
      enddo
   enddo
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( varea(J,K,numphi) > 0.0 ) then
            qdot = (quant(J,K,numphi) + 0.5*slope(J,K,numphi))*     &
                   varea(J,K,numphi)
         else
            qdot = (quant(J,K,1) - 0.5*slope(J,K,1))*               &
                   varea(J,K,numphi)
         endif
         retval(J,K,numphi) = retval(J,K,numphi) - deltat*volinv(J)*   &
                              qdot
         retval(J,K,1) = retval(J,K,1) + deltat*volinv(J)*qdot
      enddo
   enddo
endif

return
end subroutine flux
