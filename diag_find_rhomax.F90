!************************************************************************
!*
!*  DIAG_FIND_RHOMAX
!*
!************************************************************************
subroutine diag_find_rhomax(rhomax1,rhomax2,rholoc1,rholoc2,          &
                            irholoc1,irholoc2)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!
!   diag_find_rhomax finds the location of the maximum density on
!   the grid and then finds the maximum density location for the
!   companion star.  The labelling is that the star that holds the
!   maximum density value is star1, the companion is then star2.
!
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

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                     &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                          &
             in_neighbor, out_neighbor, root,                          &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                            &
           in_neighbor, out_neighbor, root,                            &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                     &
                        iam_on_bottom, iam_on_axis,                    &
                        iam_on_edge, down_neighbor,                    &
                        up_neighbor, in_neighbor,                      &
                        out_neighbor, root, column_num,                &
                        row_num, pe_grid, iam_root,                    &
                        REAL_SIZE, INT_SIZE

!*
!************************************************************************
!*
!*   Local Variables

real, dimension(numr_dd,numz_dd,numphi) :: temp

real, dimension(4,numr_procs*numz_procs) :: comm_val

integer, dimension(3,numr_procs*numz_procs) :: icomm_val

integer, dimension(3) :: location

integer, dimension(1) :: index

#ifdef SHORT
integer*8, dimension(MPI_STATUS_SIZE) :: istatus

integer*8 :: ierror, I
#else
integer, dimension(MPI_STATUS_SIZE) :: istatus

integer :: ierror, I
#endif

integer :: itemploc1, itemploc2

integer :: J, K, L

!*
!************************************************************************
!*
!*  Subroutine Arguments

real, dimension(3) :: rholoc1, rholoc2

real :: rhomax1, rhomax2

integer, dimension(3) :: irholoc1, irholoc2

!*
!************************************************************************
!  initialize the local variables
temp = 0.0
comm_val = 0.0
icomm_val = 0
location = 0
index = 0
istatus = 0
ierror = 0
itemploc1 = 0
itemploc2 = 0

!------------------------------------------------------------------------

! find the global maximum density and its location
! this will define what is labelled star1.  star2 
! is expected to be pi around from this point in azimuth
location = maxloc(rho_diag(rlwb:rupb,zlwb:zupb,:))
icomm_val(1,iam+1) = location(1) + 1
icomm_val(2,iam+1) = location(2) + 1
icomm_val(3,iam+1) = location(3)
comm_val(1,iam+1) = rho_diag(location(1)+1,location(2)+1,location(3))
comm_val(2,iam+1) = rhf(location(1)+1)
comm_val(3,iam+1) = zhf(location(2)+1)
comm_val(4,iam+1) = phi(location(3))

if( iam_root ) then
   do I = 1, numprocs-1
      call mpi_recv(comm_val(1:4,I+1),4,REAL_SIZE,I,1000+I,             &
                    MPI_COMM_WORLD,istatus,ierror)
   enddo
else
   call mpi_send(comm_val(1:4,iam+1),4,REAL_SIZE,root,1000+iam,         &
                 MPI_COMM_WORLD,ierror)
endif

if( iam_root ) then
   do I = 1, numprocs - 1
      call mpi_recv(icomm_val(:,I+1),3,INT_SIZE,I,2000+I,               &
                    MPI_COMM_WORLD,istatus,ierror)
   enddo
   index = maxloc(comm_val(1,:))
   rhomax1 = comm_val(1,index(1))
   rholoc1(1) = comm_val(2,index(1))
   rholoc1(2) = comm_val(3,index(1))
   rholoc1(3) = comm_val(4,index(1))
   irholoc1(1) = icomm_val(1,index(1))
   irholoc1(2) = icomm_val(2,index(1))
   irholoc1(3) = icomm_val(3,index(1))
else
   call mpi_send(icomm_val(:,iam+1),3,INT_SIZE,root,                    &
                 2000+iam,MPI_COMM_WORLD,ierror)
endif

call mpi_bcast(irholoc1,3,INT_SIZE,root,MPI_COMM_WORLD,ierror)

call mpi_bcast(rholoc1,3,REAL_SIZE,root,MPI_COMM_WORLD,ierror)

!------------------------------------------------------------------------

! form the density and potential arrays that are cshifted so
! the line of centers is along the x axis
temp = cshift(rho_diag,dim=3,shift=irholoc1(3)-1)

!------------------------------------------------------------------------

! find the maximum density for the other star
comm_val = 0.0
icomm_val = 0

location = maxloc(temp(rlwb:rupb,zlwb:zupb,                             &
                       numphi_by_four:3*numphi_by_four))
!itemploc1 is the azimuthal zone number relative to the
! cshifted density array
itemploc1 = location(3) + numphi_by_four - 1
! itemploc2 is the azimuthal zone number relative to the
! original density array
itemploc2 = mod(itemploc1 + irholoc1(3) - 1,numphi)
icomm_val(1,iam+1) = location(1) + 1
icomm_val(2,iam+1) = location(2) + 1
icomm_val(3,iam+1) = itemploc2
comm_val(1,iam+1) = temp(location(1)+1,location(2)+1,itemploc1)
comm_val(2,iam+1) = rhf(location(1)+1)
comm_val(3,iam+1) = zhf(location(2)+1)
comm_val(4,iam+1) = phi(itemploc2)

if( iam_root ) then
   do I = 1, numprocs-1
      call mpi_recv(comm_val(1:4,I+1),4,REAL_SIZE,I,3000+I,             &
                    MPI_COMM_WORLD,istatus,ierror)
   enddo
else
   call mpi_send(comm_val(1:4,iam+1),4,REAL_SIZE,root,3000+iam,         &
                 MPI_COMM_WORLD,ierror)
endif

if( iam_root ) then
   do I = 1, numprocs - 1
      call mpi_recv(icomm_val(:,I+1),3,INT_SIZE,I,4000+I,               &
                    MPI_COMM_WORLD,istatus,ierror)
   enddo
   index = maxloc(comm_val(1,:))
   rhomax2 = comm_val(1,index(1))
   rholoc2(1) = comm_val(2,index(1))
   rholoc2(2) = comm_val(3,index(1))
   rholoc2(3) = comm_val(4,index(1))
   irholoc2(1) = icomm_val(1,index(1))
   irholoc2(2) = icomm_val(2,index(1)) 
   irholoc2(3) = icomm_val(3,index(1))
else
   call mpi_send(icomm_val(:,iam+1),3,INT_SIZE,root,                    &
                 4000+iam,MPI_COMM_WORLD,ierror)
endif

call mpi_bcast(irholoc2,3,INT_SIZE,root,MPI_COMM_WORLD,ierror)

call mpi_bcast(rholoc2,3,REAL_SIZE,root,MPI_COMM_WORLD,ierror)

return
end subroutine diag_find_rhomax
