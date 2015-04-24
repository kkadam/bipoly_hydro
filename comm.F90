!********************************************************************
!*
!*  COMM
!*
!********************************************************************
subroutine comm(quant) 
implicit none
include 'runhydro.h'
include 'mpif.h'
!********************************************************************
!*
!  comm fills in guard cells between neighboring processors
!
!  the order of message passing is as follows:
!
!  -> send data at furthest physical radial slice to outward
!     neighbor to fill in their innermost guard cell slice
!
!  -> send data at innermost physical radial slice to inward
!     neighbor to fill in their outermost radial guard cell 
!     slice
!
!  -> send data at uppermost physical vertical slice to upward 
!     neighbor to fill in their lowermost vertical guard cell
!     slice
!
!  -> send data at lowermost physical vertical slice to downward
!     neighbor to fill in their uppermost vertical guard cell
!     slice
!
!*
!********************************************************************
!*
!*  Subroutine Arguments
  
real, dimension(numr_dd,numz_dd,numphi) :: quant

!*
!********************************************************************
!*
!*   Global Variables

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                  &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                      &
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

real, dimension(numz_dd,numphi) :: r_send, r_recv

real, dimension(numr_dd,numphi) :: z_send, z_recv

#ifdef SHORT
integer*8, dimension(MPI_STATUS_SIZE) :: istatus 
  
integer*8 :: message_length, ierror, irequest_send, irequest_recv
#else
integer, dimension(MPI_STATUS_SIZE) :: istatus

integer :: message_length, ierror, irequest_send, irequest_recv
#endif

integer :: J, K, L

!*
!********************************************************************
!  initialize the local variables
message_length = 0
ierror = 0
irequest_send = 0
istatus = 0

! pass quant at outer radial boundary
message_length = numphi * numz_dd  
if( iam_on_axis ) then	! iam can only send

   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         r_send(K,L) = quant(rupb,K,L)
      enddo
   enddo 
   call mpi_isend(r_send,message_length,REAL_SIZE,out_neighbor,       &
                  100,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_send, istatus, ierror)

else if( iam_on_edge ) then	! iam can only receive

   call mpi_irecv(r_recv,message_length,REAL_SIZE,in_neighbor,         &
                  100,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         quant(rlwb-1,K,L) = r_recv(K,L)
      enddo
   enddo

else				! iam sends then receives

   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         r_send(K,L) = quant(rupb,K,L)
      enddo
   enddo
   call mpi_irecv(r_recv,message_length,REAL_SIZE,in_neighbor,        &
                  100,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_isend(r_send,message_length,REAL_SIZE,out_neighbor,       &
                  100,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         quant(rlwb-1,K,L) = r_recv(K,L)
      enddo
   enddo
   call mpi_wait(irequest_send, istatus, ierror)

endif

! pass quant at inner radial boundary 
if( iam_on_axis ) then	! iam can only receive

   call mpi_irecv(r_recv,message_length,REAL_SIZE,out_neighbor,         &
                 200,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         quant(rupb+1,K,L) = r_recv(K,L)
      enddo
   enddo	

else if( iam_on_edge ) then	! iam can only send

   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         r_send(K,L) = quant(rlwb,K,L)
      enddo
   enddo	
   call mpi_isend(r_send,message_length,REAL_SIZE,in_neighbor,         &
                  200,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_send, istatus, ierror)

else				! iam sends then receives

   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         r_send(K,L) = quant(rlwb,K,L)
      enddo
   enddo
   call mpi_irecv(r_recv,message_length,REAL_SIZE,out_neighbor,        &
                  200,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_isend(r_send,message_length,REAL_SIZE,in_neighbor,         &
                  200,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do K = zlwb-1, zupb+1
         quant(rupb+1,K,L) = r_recv(K,L)
      enddo
   enddo 
   call mpi_wait(irequest_send, istatus, ierror)

endif

!  pass quant at top vertical boundary
message_length = numphi*numr_dd
if( iam_on_top ) then		! iam can only receive

   call mpi_irecv(z_recv,message_length,REAL_SIZE,down_neighbor,    &
                  300,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         quant(J,zlwb-1,L) = z_recv(J,L)
      enddo
   enddo 

else if( iam_on_bottom ) then	! iam can only send

   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         z_send(J,L) = quant(J,zupb,L)
      enddo
   enddo 
   call mpi_isend(z_send,message_length,REAL_SIZE,up_neighbor,      &
                  300,MPI_COMM_WORLD,irequest_send,ierror) 
   call mpi_wait(irequest_send, istatus, ierror)

else				! iam sends then receives

   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         z_send(J,L) = quant(J,zupb,L)
      enddo
   enddo
   call mpi_irecv(z_recv,message_length,REAL_SIZE,down_neighbor,     &
                  300,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_isend(z_send,message_length,REAL_SIZE,up_neighbor,       &
                  300,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         quant(J,zlwb-1,L) = z_recv(J,L)
      enddo
   enddo
   call mpi_wait(irequest_send, istatus, ierror)

endif

!  pass quant at bottom vertical boundary
if( iam_on_top ) then		! iam can only send

   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         z_send(J,L) = quant(J,zlwb,L)
      enddo
   enddo
   call mpi_isend(z_send,message_length,REAL_SIZE,down_neighbor,      &
                  400,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_send, istatus, ierror)

else if( iam_on_bottom ) then	! iam can only receive

   call mpi_irecv(z_recv,message_length,REAL_SIZE,up_neighbor,        &
                  400,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         quant(J,zupb+1,L) = z_recv(J,L)
      enddo
   enddo

else				! iam  sends then receives

   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         z_send(J,L) = quant(J,zlwb,L)
      enddo
   enddo
   call mpi_irecv(z_recv,message_length,REAL_SIZE,up_neighbor,         &
                  400,MPI_COMM_WORLD,irequest_recv,ierror)
   call mpi_isend(z_send,message_length,REAL_SIZE,down_neighbor,       &
                  400,MPI_COMM_WORLD,irequest_send,ierror)
   call mpi_wait(irequest_recv, istatus, ierror)
   do L = philwb, phiupb
      do J = rlwb-1, rupb+1
         quant(J,zupb+1,L) = z_recv(J,L)
      enddo
   enddo
   call mpi_wait(irequest_send, istatus, ierror)

endif

return
end subroutine comm
