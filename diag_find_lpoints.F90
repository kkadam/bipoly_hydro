!************************************************************************
!*
!*  DIAG_FIND_LPOINTS
!*
!************************************************************************
subroutine diag_find_lpoints(factor,rholoc1,irholoc1,rholoc2,irholoc2,  &
                             lpoints,roche_vol1,roche_vol2)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!   diag_find_lpoints finds the locations of the Lagrange points
!   along the line of centers, the outer extent of each *s Roche
!   lobe, and the potential at each of these points.  Also 
!   calculates the total volume occupied by each lobe. 
!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: rho_diag, potr, vr, vz, vphi
real, dimension(numr_dd,numphi) :: x, y, xin, yin
common /diag/ rho_diag, potr, vr, vz, vphi, x, y, xin, yin

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                      &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                           &
             in_neighbor, out_neighbor, root,                           &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                             &
           in_neighbor, out_neighbor, root,                             &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                      &
                        iam_on_bottom, iam_on_axis,                     &
                        iam_on_edge, down_neighbor,                     &
                        up_neighbor, in_neighbor,                       &
                        out_neighbor, root, column_num,                 &
                        row_num, pe_grid, iam_root,                     &
                        REAL_SIZE, INT_SIZE

!*
!************************************************************************
!*
!*   Local Variables

real, dimension(numr_dd,numz_dd,numphi) :: sum_it

real, dimension(12,numr_procs*numz_procs) :: comm_val

real, dimension(2) :: pass, passed

real :: rchtest, curvature, error_sum

integer, dimension(3) :: location

integer, dimension(1) :: index

#ifdef SHORT
integer*8, dimension(MPI_STATUS_SIZE) :: istatus

integer*8 :: ierror, I
#else
integer, dimension(MPI_STATUS_SIZE) :: istatus

integer :: ierror, I
#endif

integer :: phiflag, l1flag, iflag

integer :: J, K, L

!*
!************************************************************************
!*
!*  Subroutine Arguments

real, dimension(4,5) :: lpoints

real, dimension(3) :: rholoc1, rholoc2

real :: roche_vol1, roche_vol2, factor

integer, dimension(3) :: irholoc1, irholoc2

!*
!************************************************************************
!  initialize the local variables
lpoints = 0.0
roche_vol1 = 0.0
roche_vol2 = 0.0
sum_it = 0.0
pass = 0.0
passed = 0.0
rchtest = 0.0
curvature = 0.0
error_sum = 0.0
location = 0
index = 0
istatus = 0
ierror = 0
iflag = 0
phiflag = 0
l1flag = 0

!------------------------------------------------------------------------
! find the L1 point 
comm_val = 0.0
if( (zhf(zlwb) <= rholoc1(2)) .and. (rholoc1(2) <= zhf(zupb)) ) then
   K = irholoc1(2)
   L = irholoc1(3)
   do J = rlwb, rupb
      if( rhf(J) <= rholoc1(1) ) then
         rchtest = (potr(J+1,K,L) - potr(J,K,L))*                       &
                   (potr(J,K,L) - potr(J-1,K,L))
         if( rchtest < 0.0 ) then
            curvature = potr(J+1,K,L) + potr(J-1,K,L)                   &
                        - 2.0 * potr(J,K,L)
            if( curvature < 0.0 ) then
               !found L1
               comm_val(1,iam+1) = potr(J,K,L)
               comm_val(2,iam+1) = x(J,L)
               comm_val(3,iam+1) = y(J,L)
               comm_val(4,iam+1) = zhf(K)
               comm_val(5,iam+1) = 1.0
               comm_val(6,iam+1) = 1.0
               exit
            endif
         endif
      endif
   enddo 
endif

if( comm_val(1,iam+1) >= 0.0 ) then
   if( (zhf(zlwb) <= rholoc2(2)) .and. (rholoc2(2) <= zhf(zupb)) ) then
      K = irholoc2(2)
      L = irholoc2(3)
      do J = rlwb, rupb
         if( rhf(J) <= rholoc2(1) ) then
            rchtest = (potr(J+1,K,L) - potr(J,K,L))*                      &
                      (potr(J,K,L) - potr(J-1,K,L))
            if( rchtest < 0.0 ) then
               curvature = potr(J+1,K,L) + potr(J-1,K,L)                  &
                           -2.0 * potr(J,K,L)
               if( curvature < 0.0 ) then
                  ! found L1
                  comm_val(1,iam+1) = potr(J,K,L)
                  comm_val(2,iam+1) = x(J,L)
                  comm_val(3,iam+1) = y(J,L)
                  comm_val(4,iam+1) = zhf(K)
                  comm_val(5,iam+1) = 2.0
                  comm_val(6,iam+1) = 1.0
                  exit
               endif
            endif
         endif
      enddo
   endif 
endif
 
if( iam_root ) then
   do I = 1, numprocs - 1
      call mpi_recv(comm_val(1:6,I+1),6,REAL_SIZE,I,1000+I,              &
                    MPI_COMM_WORLD,istatus,ierror)
   enddo

   error_sum = sum(comm_val(6,:))

   if( error_sum > 1.0 ) then
      ! found multiple L1 points
      write(6,*) 'Error: found multiple Lagrange points!'
      write(6,*) comm_val
      comm_val(1:5,1) = 0.0
   else
      ! everything worked
      do I = 1, numprocs
         if( comm_val(6,I) > 0.0 ) then
            comm_val(1,1) = comm_val(1,I)
            comm_val(2,1) = comm_val(2,I)
            comm_val(3,1) = comm_val(3,I)
            comm_val(4,1) = comm_val(4,I)
            comm_val(5,1) = comm_val(5,I)
            exit
         endif
      enddo

   endif
else
   call mpi_send(comm_val(1:6,iam+1),6,REAL_SIZE,root,1000+iam,         &
                 MPI_COMM_WORLD,ierror)
endif

call mpi_bcast(comm_val(1:5,1),5,REAL_SIZE,root,MPI_COMM_WORLD,ierror)
lpoints(1,1) = comm_val(1,1)
lpoints(2,1) = comm_val(2,1)
lpoints(3,1) = comm_val(3,1)
lpoints(4,1) = comm_val(4,1)
       
if( comm_val(5,1) == 0.0 ) then
   l1flag = 0
else if( comm_val(5,1) == 1.0 ) then
   l1flag = 1
else if( comm_val(5,1) == 2.0 ) then
   l1flag = 2
endif
       
!------------------------------------------------------------------------
! find other quantities for the Roche lobes: position of the
! L2, L3 points and the outer edges of the Roche lobes.  These 
! values may be useful for diagnosing problems with the
! Roche lobe specific integrations.

comm_val = 0.0

if( (zhf(zlwb) <= rholoc1(2)) .and. (rholoc1(2) <= zhf(zupb)) ) then
   K = irholoc1(2)
   L = irholoc1(3)
   do J = rlwb, rupb
      if( rhf(J) > rholoc1(1) ) then
         rchtest = (potr(J-1,K,L) - potr(J,K,L))*                       &
                   (potr(J,K,L) - potr(J+1,K,L))
         if( rchtest < 0.0 ) then
            curvature = potr(J-1,K,L) + potr(J+1,K,L) -                 &
                        2.0*potr(J,K,L)
            if( curvature < 0.0 ) then
               ! found 'L2'
               comm_val(1,iam+1) = x(J,L)
               comm_val(2,iam+1) = y(J,L)
               comm_val(3,iam+1) = potr(J,K,L)
            endif
         endif
         ! find outer boundary of Primary's Roche Lobe
         if( (potr(J,K,L) < lpoints(1,1)) .and.                         &
             (lpoints(1,1) < potr(J+1,K,L)) ) then
             comm_val(4,iam+1) = x(J,L)
             comm_val(5,iam+1) = y(J,L)
             comm_val(6,iam+1) = potr(J,K,L)
         endif 
      endif 
   enddo 
endif

!  look for quantities on secondary side
if( (zhf(zlwb) <= rholoc2(2)) .and. (rholoc2(2) <= zhf(zupb)) ) then
   K = irholoc2(2)
   L = irholoc2(3)
   do J = rlwb, rupb
      if( rhf(J) > rholoc2(1) ) then
         rchtest = (potr(J-1,K,L) - potr(J,K,L))*                         &
                   (potr(J,K,L) - potr(J+1,K,L))
         if( rchtest < 0.0 ) then
            curvature = potr(J-1,K,L) + potr(J+1,K,L) -                   &
                        2.0*potr(J,K,L)
            if( curvature < 0.0 ) then
               ! found 'L3'
               comm_val(7,iam+1) = x(J,L)
               comm_val(8,iam+1) = y(J,L)
               comm_val(9,iam+1) = potr(J,K,L)
            endif 
         endif
         ! find outer boundary of secondary's Roche lobe
         if( (potr(J,K,L) < lpoints(1,1)) .and.                            &
             (lpoints(1,1) < potr(J+1,K,L)) ) then
             comm_val(10,iam+1) = x(J,L) 
             comm_val(11,iam+1) = y(J,L)
             comm_val(12,iam+1) = potr(J,K,L)
         endif
      endif
   enddo 
endif

! communicate these values
if( iam_root ) then
   do I = 1, numprocs - 1
      call mpi_recv(comm_val(:,I),12,REAL_SIZE,I,2000+I,            &
                    MPI_COMM_WORLD,istatus,ierror)
      if( comm_val(1,I) /= 0.0 ) then
         lpoints(1,2) = comm_val(3,I)
         lpoints(2,2) = comm_val(1,I)
         lpoints(3,2) = comm_val(2,I)
         lpoints(4,2) = rholoc1(2)
      endif 
      if( comm_val(4,I) /= 0.0 ) then
         lpoints(1,3) = comm_val(6,I)
         lpoints(2,3) = comm_val(4,I)
         lpoints(3,3) = comm_val(5,I)
         lpoints(4,3) = rholoc1(2)
      endif 
      if( comm_val(7,I) /= 0.0 ) then
         lpoints(1,4) = comm_val(9,I)
         lpoints(2,4) = comm_val(7,I)
         lpoints(3,4) = comm_val(8,I)
         lpoints(4,4) = rholoc2(2)
      endif 
      if( comm_val(10,I) /= 0.0 ) then
         lpoints(1,5) = comm_val(12,I)
         lpoints(2,5) = comm_val(10,I)
         lpoints(3,5) = comm_val(11,I)
         lpoints(4,5) = rholoc2(2)
      endif       
   enddo 
else
   call mpi_send(comm_val(:,iam+1),12,REAL_SIZE,root,2000+iam,      &
                 MPI_COMM_WORLD,ierror)
endif

!------------------------------------------------------------------------
! sum up the total volume in each Roche lobe
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( potr(J,K,L) <= lpoints(1,1) ) then
            sum_it(J,K,L) = rhf(J)
         else
            sum_it(J,K,L) = 0.0
         endif	
      enddo
   enddo
enddo

call dsummer(sum_it,x,l1flag,lpoints(2,1),roche_vol1,roche_vol2)

pass(1) = roche_vol1
pass(2) = roche_vol2

call mpi_reduce(pass,passed,2,REAL_SIZE,MPI_SUM,root,                  &
                MPI_COMM_WORLD,ierror)

if( iam_root ) then
   roche_vol1 = factor*passed(1)
   roche_vol2 = factor*passed(2)
endif

!--------------------------------------------------------------------------
 
return
end subroutine diag_find_lpoints
