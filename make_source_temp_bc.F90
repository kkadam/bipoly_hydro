!************************************************************************
!*
!*  MAKE_SOURCE_TEMP_BC
!*
!************************************************************************
subroutine make_source_temp_bc
implicit none
include 'runhydro.h'
!************************************************************************
!*
!
!*
!**********************************************************************
!*
!*   Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: temps, tempa
common /source_temp/ temps, tempa

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
   do L = philwb, phiupb
      do K = zlwb, zupb
         temps(rlwb,K,L) = 0.0
      enddo
   enddo
   do K = zlwb, zupb
      sum = 0.0
      do L = philwb, phiupb
         sum = sum + tempa(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         tempa(rlwb,K,L) = sum
      enddo
   enddo
   if( isym == 3 ) then
      do L = philwb, phiupb
         do K = zlwb, zupb
            temps(rlwb-1,K,L) = - temps(rlwb+1,K,L)
            tempa(rlwb-1,K,L) =   tempa(rlwb,K,L)
	 enddo
      enddo
   else
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            temps(rlwb-1,K,L)               = - temps(rlwb+1,K,L+numphi_by_two)
            temps(rlwb-1,K,L+numphi_by_two) = - temps(rlwb+1,K,L)
	 enddo
      enddo
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            tempa(rlwb-1,K,L)               = tempa(rlwb,K,L+numphi_by_two)
            tempa(rlwb-1,K,L+numphi_by_two) = tempa(rlwb,K,L)
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
               temps(J,zlwb-1,L) = temps(J,zlwb,L)
               tempa(J,zlwb-1,L) = tempa(J,zlwb,L)
	    enddo
	 enddo
      else if( boundary_condition(1) == 2 ) then
         ! free boundary condition
	 do L = philwb, phiupb
	    do J = rlwb, rupb
               temps(J,zlwb-1,L) = temps(J,zlwb,L)
               tempa(J,zlwb-1,L) = tempa(J,zlwb,L)
	    enddo
	 enddo
      else  
         ! dirichlet bc case
	 do L = philwb, phiupb
	    do J = rlwb, rupb
	       temps(J,zlwb-1,L) = 0.0
	       tempa(J,zlwb-1,L) = 0.0
	    enddo
	 enddo
      endif
   else
      ! if isym /= 1 then force wall boundary
      ! condition at the bottom of the grid
      do L = philwb, phiupb
         do J = rlwb, rupb
            temps(J,zlwb-1,L) = temps(J,zlwb,L)
            tempa(J,zlwb-1,L) = tempa(J,zlwb,L)
         enddo
      enddo
   endif
endif

if( iam_on_edge ) then
   if( boundary_condition(2) == 1 ) then
      ! wall boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            temps(rupb+1,K,L) = - temps(rupb-1,K,L)
            temps(rupb,K,L)   =   0.0
            tempa(rupb+1,K,L) =   tempa(rupb,K,L)
         enddo
      enddo
   else if( boundary_condition(2) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            temps(rupb+1,K,L) = temps(rupb,K,L)
            tempa(rupb+1,K,L) = tempa(rupb,K,L)
         enddo
      enddo
   else	! Dirichlet bc
      do L = philwb, phiupb
         do K = zlwb, zupb
            temps(rupb+1,K,L) = 0.0
            temps(rupb,K,L)   = 0.5*(abs(temps(rupb,K,L)) +      &
                                         temps(rupb,K,L))
	    tempa(rupb+1,K,L) = 0.0
         enddo
      enddo
   endif
endif

if( iam_on_top ) then
   if( boundary_condition(3) == 1 ) then
      ! wall boundary condition 
      do L = philwb, phiupb
         do J = rlwb, rupb
            temps(J,zupb+1,L) = temps(J,zupb,L)
            tempa(J,zupb+1,L) = tempa(J,zupb,L)
         enddo
      enddo
   else if( boundary_condition(3) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            temps(J,zupb+1,L) = temps(J,zupb,L)
            tempa(J,zupb+1,L) = tempa(J,zupb,L)
         enddo
      enddo
   else	! Dirichlet bc
      do L = philwb, phiupb
         do J = rlwb, rupb
            temps(J,zupb+1,L) = 0.0
            tempa(J,zupb+1,L) = 0.0
	 enddo
      enddo
   endif
endif

return
end subroutine make_source_temp_bc
