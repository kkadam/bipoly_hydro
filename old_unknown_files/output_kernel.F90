!*************************************************************************
!*
!*  OUTPUT_KERNEL
!*
!*************************************************************************
!
!  2/28/2005 added a timer for the I/O and also implemented a temporary
!  array that excludes the ghost zone values.  When using the Intel 
!  compiler v 8.0 noted that the output will die if there is an array
!  syntax statement in the write command (something that looks like
!  write(...) rho(rlwb:rupb,zlwb:zupb,:) )
! 
!*************************************************************************
subroutine output_kernel(q,template,frnum)
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

character(len=54) :: filename

#ifdef SHORT
integer*8 :: token_tag, ierror

integer*8, dimension(MPI_STATUS_SIZE) :: istatus
#else
integer :: token_tag, ierror

integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif

integer :: token, record_length

integer :: J, K, L


!*
!*************************************************************************
!*
!*   Subroutine Arguments

real, dimension(numr_dd,numz_dd,numphi) :: q

character(len=50) :: template

integer :: frnum

!*
!*************************************************************************
! Initialize local variables
record_length = 0
token = 1
token_tag = 100
ierror = 0
istatus = 0

do L = 1, numphi
   do K = 1, numz_dd-2
      do J = 1, numr_dd-2
         output_array(J,K,L) = q(J+1,K+1,L)
      enddo
   enddo
enddo

! let fortan determine the size of each record for
! these files.  The record length should be the
! size in bytes of each of the following data sections
inquire(iolength=record_length) output_array

! make sure everyone is on the same page before we start
call mpi_barrier(MPI_COMM_WORLD,ierror)

time1 = mpi_wtime()

! write out the density array
write(filename,'(a,i4)') trim(template),frnum
if( iam_root ) then

   open(unit=50,file=trim(filename),form='unformatted',access='direct',recl=record_length) 

   write(50,rec=iam+1) output_array

   close(50)

   call mpi_send(token,1,INT_SIZE,iam+1,token_tag,MPI_COMM_WORLD,ierror)

else
   call mpi_recv(token,1,INT_SIZE,iam-1,token_tag,MPI_COMM_WORLD,istatus,ierror)

   open(unit=50,file=trim(filename),form='unformatted',access='direct',recl=record_length)

   write(50,rec=iam+1) output_array

   close(50)

   if( iam /= numprocs-1 ) then
      call mpi_send(token,1,INT_SIZE,iam+1,token_tag,MPI_COMM_WORLD,ierror)
   endif
endif

! let everyone catch up before we plow ahead
call mpi_barrier(MPI_COMM_WORLD,ierror)

time2 = mpi_wtime()

if ( iam_root ) then
   write(6,*) 'wrote out dump: ', trim(filename), ' in time: ', time2 - time1
endif

return
end subroutine output_kernel
