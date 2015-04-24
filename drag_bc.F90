!************************************************************************
!*
!*  DRAG_BC
!*
!************************************************************************
subroutine drag_bc
implicit none
include 'runhydro.h'
!************************************************************************
!*
!
!*
!**********************************************************************
!*
!*   Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

logical :: iam_on_top, iam_on_bottom, iam_on_axis,           &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                &
             in_neighbor, out_neighbor, root,                &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                  &
           in_neighbor, out_neighbor, root,                  &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,           &
                        iam_on_bottom, iam_on_axis,          &
                        iam_on_edge, down_neighbor,          &
                        up_neighbor, in_neighbor,            &
                        out_neighbor, root, column_num,      &
                        row_num, pe_grid, iam_root,          &
                        REAL_SIZE, INT_SIZE


!*
!************************************************************************
!*
!*   Local Variables

real :: sum

integer :: J, K, L 

!*
!************************************************************************
!  initialize the local variables

if( iam_on_axis ) then
   do K = zlwb, zupb
      sum = 0.0
      do L = philwb, phiupb
         sum = sum + jn(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         jn(rlwb,K,L) = sum
      enddo
      sum = 0.0
      do L = philwb, phiupb
         sum = sum + a(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         a(rlwb,K,L) = sum
      enddo
   enddo
   if( isym == 3 ) then
      do L = philwb, phiupb
         do K = zlwb, zupb
            jn(rlwb-1,K,L) =   jn(rlwb,K,L)
            a(rlwb-1,K,L)  =    a(rlwb,K,L)
	 enddo
      enddo
   else
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            jn(rlwb-1,K,L)               =  jn(rlwb,K,L+numphi_by_two)
            jn(rlwb-1,K,L+numphi_by_two) =  jn(rlwb,K,L)
	 enddo
      enddo
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            a(rlwb-1,K,L)               = a(rlwb,K,L+numphi_by_two)
            a(rlwb-1,K,L+numphi_by_two) = a(rlwb,K,L)
	 enddo
      enddo
   endif
endif

if( iam_on_bottom ) then
   if( isym == 1 ) then
      if( boundary_condition(1) == 1 ) then
         ! wall boundary condition
	 do L = philwb, phiupb
	    do J = rlwb, rupb
               jn(J,zlwb-1,L) = jn(J,zlwb,L)
               a(J,zlwb-1,L)  =  a(J,zlwb,L)
	    enddo
	 enddo
      else if( boundary_condition(1) == 2 ) then
         ! free boundary condition
	 do L = philwb, phiupb
	    do J = rlwb, rupb
               jn(J,zlwb-1,L) = jn(J,zlwb,L)
               a(J,zlwb-1,L)  =  a(J,zlwb,L)
	    enddo
	 enddo
      else  
         ! dirichlet bc case
	 do L = philwb, phiupb
	    do J = rlwb, rupb
	       jn(J,zlwb-1,L) = 0.0
	       a(J,zlwb-1,L)  = 0.0
	    enddo
	 enddo
      endif
   else
      ! if isym /= 1 then force wall boundary
      ! condition at the bottom of the grid
      do L = philwb, phiupb
         do J = rlwb, rupb
            jn(J,zlwb-1,L) = jn(J,zlwb,L)
            a(J,zlwb-1,L)  =  a(J,zlwb,L)
         enddo
      enddo
   endif
endif

if( iam_on_edge ) then
   if( boundary_condition(2) == 1 ) then
      ! wall boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            jn(rupb+1,K,L) = jn(rupb,K,L)
            a(rupb+1,K,L)  =  a(rupb,K,L)
         enddo
      enddo
   else if( boundary_condition(2) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            jn(rupb+1,K,L) = jn(rupb,K,L)
            a(rupb+1,K,L)  =  a(rupb,K,L)
         enddo
      enddo
   else	! Dirichlet bc
      do L = philwb, phiupb
         do K = zlwb, zupb
            jn(rupb+1,K,L) = 0.0
	    a(rupb+1,K,L)  = 0.0
         enddo
      enddo
   endif
endif

if( iam_on_top ) then
   if( boundary_condition(3) == 1 ) then
      ! wall boundary condition 
      do L = philwb, phiupb
         do J = rlwb, rupb
            jn(J,zupb+1,L) = jn(J,zupb,L)
            a(J,zupb+1,L)  = a(J,zupb,L)
         enddo
      enddo
   else if( boundary_condition(3) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            jn(J,zupb+1,L) = jn(J,zupb,L)
            a(J,zupb+1,L)  = a(J,zupb,L)
         enddo
      enddo
   else	! Dirichlet bc
      do L = philwb, phiupb
         do J = rlwb, rupb
            jn(J,zupb+1,L) = 0.0
            a(J,zupb+1,L)  = 0.0
	 enddo
      enddo
   endif
endif

return
end subroutine drag_bc
