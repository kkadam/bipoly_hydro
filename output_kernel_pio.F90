!*************************************************************************
!*
!*  OUTPUT_KERNEL
!*
!*************************************************************************
subroutine output_kernel(q, record_number, unit_number)
implicit none
include 'runhydro.h'
include 'mpif.h'
!*************************************************************************
!*
!*  Global Variables

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
!************************************************************************* 
!*
!*   Local variables

real, dimension(numr_dd-2,numz_dd-2,numphi) :: output_array

real :: time1, time2

integer :: J, K, L

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

!*
!*************************************************************************
!*
!*   Subroutine Arguments

real, dimension(numr_dd,numz_dd,numphi) :: q

integer :: record_number, unit_number


!*
!*************************************************************************

ierror = 0

! copy the data into the output array
do L = 1, numphi
   do K = 1, numz_dd-2
      do J = 1, numr_dd-2
         output_array(J,K,L) = q(J+1,K+1,L)
      enddo
   enddo
enddo

! make sure everyone is on the same page before we start
call mpi_barrier(MPI_COMM_WORLD,ierror)

time1 = mpi_wtime()

! write out the array
write(unit_number,rec=record_number) output_array

time2 = mpi_wtime()

! let everyone catch up before we plow ahead
call mpi_barrier(MPI_COMM_WORLD,ierror)

time2 = mpi_wtime()

if ( iam_root ) then
   write(6,*) 'completed dump: ', record_number, unit_number,  ' in time: ', time2 - time1
endif

return
end subroutine output_kernel
