!*************************************************************************
!*
!*  OUTPUT_DEBUG
!*
!*************************************************************************
subroutine output_debug(frnum,output_array)
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

character(len=54) :: output_file

character(len=50) :: output_template

integer :: output_recl, token

#ifdef SHORT
integer*8 :: token_tag, ierror

integer*8, dimension(MPI_STATUS_SIZE) :: istatus
#else
integer :: token_tag, ierror

integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif

!*
!*************************************************************************
!*
!*   Subroutine Arguments

real, dimension(numr_dd,numz_dd,numphi) :: output_array

!integer, dimension(numr_dd,numz_dd,numphi) :: tmp_array

integer :: frnum

!*
!*************************************************************************
! Initialize local variables
output_recl = 0
token = 1
token_tag = 100
ierror = 0
istatus = 0
output_template =   'file'

!coercion to float, output of int arrays is weird on the T3E...
!output_array = tmp_array

! let fortan determine the size of each record for
! these files.  The record length should be the
! size in bytes of each of the following data sections
inquire(iolength=output_recl) output_array(rlwb:rupb,zlwb:zupb,:)

! create the filenames for the files every pe is going
write(output_file,'(a,i4)') trim(output_template),frnum

! pass a token from processor to processor, letting
! each one write to the files in turn
if( iam_root ) then

   open(unit=60,file=trim(output_file),form='unformatted',               &
        access='direct',recl=output_recl) 

   write(60,rec=iam+1) output_array(rlwb:rupb,zlwb:zupb,:)

   close(60)

   call mpi_send(token,1,INT_SIZE,iam+1,token_tag,                       &
                 MPI_COMM_WORLD,ierror)
else
   call mpi_recv(token,1,INT_SIZE,iam-1,token_tag,                       &
                 MPI_COMM_WORLD,istatus,ierror)

   open(unit=60,file=trim(output_file),form='unformatted',               &
        access='direct',recl=output_recl)

   write(60,rec=iam+1) output_array(rlwb:rupb,zlwb:zupb,:)

   close(60)

   if( iam /= numprocs-1 ) then
      call mpi_send(token,1,INT_SIZE,iam+1,token_tag,                    &
                    MPI_COMM_WORLD,ierror)
   endif
endif

! let everyone catch up before we plow ahead
call mpi_barrier(MPI_COMM_WORLD,ierror)


return
end subroutine output_debug
