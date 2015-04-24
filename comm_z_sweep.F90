!********************************************************************
!*
!*  COMM_Z_SWEEP
!*
!********************************************************************
subroutine comm_z_sweep(quant) 
implicit none
include 'runhydro.h'
include 'mpif.h'
!********************************************************************
!*
!  comm_z_sweep is closely related to comm and is used in the
!  adi solver, helmadi.  comm_z_sweep updates data across 
!  horizontal boundaries in the pe grid when the input array
!  has all vertical data in local memory, the radial index block 
!  distributed across the vertical dimension of the pe grid and
!  the azimuthal index block distributed across the horizontal
!  dimension of the pe grid
!*
!********************************************************************
!*
!*  Subroutine Arguments
  
real, dimension(numr_dd_z,numz,numphi_dd) :: quant

!*
!********************************************************************
!*
!*   Global Variables

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                  &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                       &
             in_neighbor, out_neighbor, root,                       &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                         &
           in_neighbor, out_neighbor, root,                         &
           REAL_SIZE, INT_SIZE, numprocs 
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                  &
                        iam_on_bottom, iam_on_axis,                 &
                        iam_on_edge, down_neighbor,                 &
                        up_neighbor, in_neighbor,                   &
                        out_neighbor, root, column_num,             &
                        row_num, pe_grid, iam_root,                 &
                        REAL_SIZE, INT_SIZE

!*
!********************************************************************
!*
!*   Local Variables

real, dimension(numz, numphi_dd) :: r_send, r_recv

#ifdef SHORT
integer*8, dimension(MPI_STATUS_SIZE) :: istatus 
  
integer*8 :: message_length, ierror, irequest_send, irequest_recv
#else
integer, dimension(MPI_STATUS_SIZE) :: istatus

integer :: message_length, ierror, irequest_send, irequest_recv
#endif

integer :: K, L

!*
!********************************************************************
!  initialize the local variables
r_send = 0.0
r_recv = 0.0
message_length = 0
ierror = 0
irequest_send = 0
istatus = 0

! pass quant at top vertical boundary
message_length = numphi_dd * numz  
if( iam_on_bottom ) then	! iam can only send

   do L = 1, numphi_dd
      do K = 1, numz
         r_send(K,L) = quant(numr_dd_z-1,K,L)
      enddo
   enddo
   call mpi_isend(r_send,message_length,REAL_SIZE,up_neighbor,      &
                  100,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_send, istatus, ierror)

else if( iam_on_top ) then	! iam can only receive

   call mpi_irecv(r_recv,message_length,REAL_SIZE,down_neighbor,      &
                 100,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = 1, numphi_dd
      do K = 1, numz
         quant(1,K,L) = r_recv(K,L)
      enddo
   enddo

else				! iam sends then receives

   do L = 1, numphi_dd
      do K = 1, numz
         r_send(K,L) = quant(numr_dd_z-1,K,L)
      enddo
   enddo
   call mpi_irecv(r_recv,message_length,REAL_SIZE,down_neighbor,     &
                 100,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_isend(r_send,message_length,REAL_SIZE,up_neighbor,       &
                  100,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = 1, numphi_dd
      do K = 1, numz
         quant(1,K,L) = r_recv(K,L)
      enddo
   enddo
   call mpi_wait(irequest_send, istatus, ierror)

endif

! pass quant at bottom vertical boundary 
if( iam_on_bottom ) then	! iam can only receive

   call mpi_irecv(r_recv,message_length,REAL_SIZE,up_neighbor,        &
                 200,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = 1, numphi_dd
      do K = 1, numz
         quant(numr_dd_z,K,L) = r_recv(K,L)
      enddo
   enddo

else if( iam_on_top ) then	! iam can only send

   do L = 1, numphi_dd
      do K = 1, numz
         r_send(K,L) = quant(2,K,L)
      enddo
   enddo
   call mpi_isend(r_send,message_length,REAL_SIZE,down_neighbor,      &
                 200,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_send, istatus, ierror)

else				! iam sends then receives

   do L = 1, numphi_dd
      do K = 1, numz
         r_send(K,L) = quant(2,K,L)
      enddo
   enddo
   call mpi_irecv(r_recv,message_length,REAL_SIZE,up_neighbor,         &
                 200,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_isend(r_send,message_length,REAL_SIZE,down_neighbor,       &
                  200,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = 1, numphi_dd
      do K = 1, numz
         quant(numr_dd_z,K,L) = r_recv(K,L)
      enddo
   enddo
   call mpi_wait(irequest_send, istatus, ierror)

endif

return
end subroutine comm_z_sweep
