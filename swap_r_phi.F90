!*************************************************************************
!*
!*  SWAP_R_PHI
!*
!*************************************************************************
subroutine swap_r_phi(xr, xphi)
implicit none
include 'runhydro.h'
include 'mpif.h'
!*************************************************************************
!*
!  swap_phi_r takes an array that is dstributed with:
!                              z in blocks across numz_procs
!                              phi in blocks across numr_procs
!                              r in local memory
!
!  and trasforms it into an array distributed with:
!                              z in blocks across numz_procs
!                              r in blocks across numr_procs
!                              phi in local memory
!
!  communication proceeds independently for each row of processors
!
!  assume that numphi is evenly divisible by numr_procs
!*
!*************************************************************************
!*
!*  Subroutine Arguments

real, dimension(numr, numz_dd, numphi_dd) :: xr

real, dimension(numr_dd, numz_dd, numphi) :: xphi

!*
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
!*  Local Variables

real, dimension(numr_dd,numz_dd,numphi_dd) :: send_buffer, receive_buffer

#ifdef SHORT
integer*8, dimension(numr_procs) :: send_to, receive_from

integer*8, dimension(MPI_STATUS_SIZE) :: istatus 

integer*8 :: I, ierror, message_length, irequest_send, temp_int,         &
             irequest_recv
#else
integer, dimension(numr_procs) :: send_to, receive_from

integer, dimension(MPI_STATUS_SIZE) :: istatus

integer :: I, ierror, message_length, irequest_send, irequest_recv
#endif

integer :: iam_mod, send_to_mod, receive_from_mod, index, q, counter,   &
           J, K, L

!*
!*************************************************************************
!  initialize the local variables
send_to = 0
receive_from = 0
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
xphi = 0.0

!  first copy iam's own section of xphi into xr
#ifdef SHORT
temp_int = numr_procs
iam_mod = mod(iam,temp_int)
#else
iam_mod = mod(iam,numr_procs)
#endif

do L = 1, numphi_dd
   do K = 1, numz_dd
      do J = 1, numr_dd
         xphi(J,K,L+iam_mod*numphi_dd) = xr(J+iam_mod*(numr_dd-2),K,L)
      enddo
   enddo
enddo 

message_length = numr_dd * numz_dd * numphi_dd

index = column_num + 1

do I = 1, numr_procs - 1

   ! cshift by + I
   do q = 1, numr_procs - I
      send_to(q) = pe_grid(q+I,row_num+1)
   enddo
   counter = 1
   do q = numr_procs - I + 1, numr_procs
      send_to(q) = pe_grid(counter,row_num+1)
      counter = counter + 1
   enddo

   ! cshift by - I
   do q = numr_procs, I + 1, -1
      receive_from(q) = pe_grid(q-I,row_num+1)
   enddo
   counter = numr_procs
   do q = I, 1, -1
      receive_from(q) = pe_grid(counter,row_num+1)
      counter = counter - 1
   enddo

#ifdef SHORT
   send_to_mod = mod(send_to(index),temp_int)

   receive_from_mod = mod(receive_from(index),temp_int)
#else
   send_to_mod = mod(send_to(index),numr_procs)

   receive_from_mod = mod(receive_from(index),numr_procs)
#endif

   do L = 1, numphi_dd
      do K = 1, numz_dd
         do J = 1, numr_dd
            send_buffer(J,K,L) = xr(J+send_to_mod*(numr_dd-2),K,L)
         enddo
      enddo
   enddo

   call mpi_irecv(receive_buffer, message_length, REAL_SIZE, receive_from(index), &
                 I*2000+iam, MPI_COMM_WORLD, irequest_recv, ierror)

   call mpi_isend(send_buffer, message_length, REAL_SIZE, send_to(index),         &
                  I*2000+send_to(index), MPI_COMM_WORLD, irequest_send, ierror)

   call mpi_wait(irequest_recv, istatus, ierror)

   do L = 1, numphi_dd
      do K = 1, numz_dd
         do J = 1, numr_dd
            xphi(J,K,L+receive_from_mod*numphi_dd) = receive_buffer(J,K,L)
         enddo
      enddo
   enddo

   call mpi_wait(irequest_send, istatus, ierror)

enddo

return
end subroutine swap_r_phi
