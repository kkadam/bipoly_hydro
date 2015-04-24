!*****************************************************************************
!*
!*  SWAP_Z_R
!*
!*****************************************************************************
subroutine swap_z_r(xz, xr)
implicit none
include 'runhydro.h'
include 'mpif.h'
!*****************************************************************************
!*
!  swap_z_r takes an array that is distributed with:
!                            r in blocks across numz_procs
!                            phi in blocks across numr_procs
!                            z in local memory
!
!  and transforms it into an array distributed with:
!                           z in blocks across numz_procs
!                           phi in blocks across numr_procs
!                           r in local memory
!
!  communication proceeds independently for each column of 
!  processors
!
!  assume that numphi is divisible by numr_procs
!
!  assume that (numr-2) is evenly divisible by numz_procs
!
!*
!*****************************************************************************
!*
!*  Subroutine Arguments

real, dimension(numr_dd_z, numz, numphi_dd) :: xz

real, dimension(numr, numz_dd, numphi_dd) :: xr

!*
!*****************************************************************************
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
!*****************************************************************************
!*
!*  Local Variables
 
real, dimension(numr_dd_z,numz_dd,numphi_dd) :: send_buffer, receive_buffer

#ifdef SHORT 
integer*8, dimension(numz_procs) :: send_to, receive_from

integer*8, dimension(MPI_STATUS_SIZE) :: istatus

integer*8 :: I, ierror, message_length, irequest_send, irequest_recv
#else
integer, dimension(numz_procs) :: send_to, receive_from

integer, dimension(MPI_STATUS_SIZE) :: istatus

integer :: I, ierror, message_length, irequest_send, irequest_recv
#endif

integer :: iam_mod, send_to_mod, receive_from_mod, index , q, counter,     &
           temp, J, K, L

!*
!*****************************************************************************
!  initialize the local variables
send_to = 0
receive_from = 0
temp = 0
iam_mod = 0
send_to_mod = 0
receive_from_mod = 0
message_length = 0
ierror = 0 
index = 0
istatus = 0
irequest_send = 0
send_buffer = 0.0
receive_buffer = 0.0
xr = 0.0

iam_mod = row_num

!  copy iam's own section of xz to xr
do L = 1, numphi_dd
   do K = 1, numz_dd
      do J = 1, numr_dd_z
         xr(J+iam_mod*(numr_dd_z-2),K,L) = xz(J,K+iam_mod*(numz_dd-2),L)
      enddo
   enddo
enddo 

!  now do the message passing to exchange the rest of xz and xr
!  have to communicate numz_procs - 1 times to fill in the rest of
!  xr on all processors.  

message_length = numr_dd_z * numz_dd * numphi_dd

index = row_num+1

do I = 1, numz_procs - 1

   ! cshift by + I
   do q = 1, numz_procs - I
      send_to(q) = pe_grid(column_num+1,q+I)
   enddo
   counter = 1
   do q = numz_procs - I + 1, numz_procs
      send_to(q) = pe_grid(column_num+1,counter)
      counter = counter + 1
   enddo

   ! cshift by - I
   do q = numz_procs, I + 1, -1
      receive_from(q) = pe_grid(column_num+1,q-I)
   enddo
   counter = numz_procs
   do q = I, 1, -1
      receive_from(q) = pe_grid(column_num+1,counter)
      counter = counter - 1
   enddo

   temp = send_to(index)
   send_to_mod = 0
   do while( temp >= numr_procs )
      send_to_mod = send_to_mod + 1
      temp = temp - numr_procs
   enddo

   temp = receive_from(index)
   receive_from_mod = 0
   do while( temp >= numr_procs )
      receive_from_mod = receive_from_mod + 1
      temp = temp - numr_procs
   enddo 

   do L = 1, numphi_dd
      do K = 1, numz_dd
         do J = 1, numr_dd_z
            send_buffer(J,K,L) = xz(J,K+send_to_mod*(numz_dd-2),L)
         enddo
      enddo
   enddo

   call mpi_irecv(receive_buffer, message_length, REAL_SIZE, receive_from(index),  &
                  I*5000+iam, MPI_COMM_WORLD, irequest_recv, ierror)

   call mpi_isend(send_buffer, message_length, REAL_SIZE, send_to(index),          &
                  I*5000+send_to(index), MPI_COMM_WORLD, irequest_send, ierror)

  
   call mpi_wait(irequest_recv, istatus, ierror)
 
   do L = 1, numphi_dd
      do K = 1, numz_dd
         do J = 1, numr_dd_z
            xr(J+receive_from_mod*(numr_dd_z-2),K,L) = receive_buffer(J,K,L)
         enddo
      enddo
   enddo

   call mpi_wait(irequest_send, istatus, ierror)

enddo
        
return
end subroutine swap_z_r
