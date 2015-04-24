!********************************************************************
!*
!*  VISC_TIME
!*
!********************************************************************
subroutine visc_time
implicit none
include 'runhydro.h'
include 'mpif.h'
!********************************************************************
!*
! calculates the limiting time step from the diffusive viscosity
! terms introduced to stabilize the scheme
!
!  deltat = 0.25  Min[ (deltaXi rho / Qii)**0.5 ]
!
!  get at this by finding the max value of Qii / ( dXi rho)
!
!*
!********************************************************************
!*
!*   Global variables
!*

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi,                          &
        drinv, dzinv, dphiinv

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                   &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                        &
             in_neighbor, out_neighbor, root,                        &
             REAL_SIZE, INT_SIZE, numprocs 
#else
integer :: iam, down_neighbor, up_neighbor,                          &
           in_neighbor, out_neighbor, root,                          &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                   &
                        iam_on_bottom, iam_on_axis,                  &
                        iam_on_edge, down_neighbor,                  &
                        up_neighbor, in_neighbor,                    &
                        out_neighbor, root, column_num,              &
                        row_num, pe_grid, iam_root,                  &
                        REAL_SIZE, INT_SIZE

!*
!********************************************************************
!*
!*   Local Variables
!*

real, dimension(numr_dd,numz_dd,numphi) :: qrr, qzz, qpp

real, dimension(numr_dd,numz_dd,numphi) :: deltj,  deltk, deltl

real, dimension(numprocs,7) :: global_diagnostics 

real :: time_max, temp_max, rhoinv

#ifdef SHORT
integer*8, dimension(MPI_STATUS_SIZE) :: istatus

integer*8 :: ierror
#else
integer, dimension(MPI_STATUS_SIZE) :: istatus

integer :: ierror
#endif

integer :: I, J, K, L

integer, dimension(1) :: index

integer, dimension(3) :: loc

!*
!*
!********************************************************************
!  initialize the local variables
global_diagnostics = 0.0
time_max = 0.0
temp_max = 0.0
loc = 0
istatus = 0
ierror = 0
index = 0

  100  format(10e16.8)

! populate  the viscous stress tensor
call stress(qrr,qzz,qpp)

!    deltx = inverse of square diffusive time across cell in a given
!            coordinate direction
!
!    deltj = qrr / dR rho
!
!
!    deltk = qzz / dZ rho
!
!
!    deltl = qpp / R dPhi rho
!
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhoinv = 1.0 / rho(J,K,L)
         deltj(J,K,L) = drinv * rhoinv * qrr(J,K,L)
         deltK(J,K,L) = dzinv * rhoinv * qzz(J,K,L)
         deltl(J,K,L) = rhfinv(J) * dphiinv * rhoinv * qpp(J,K,L)
      enddo
   enddo
enddo 

!    now each processor finds their own max value of 1 / dt
time_max = maxval(deltj(rlwb:rupb,zlwb:zupb,:))
loc = maxloc(deltj(rlwb:rupb,zlwb:zupb,:))

temp_max = maxval(deltk(rlwb:rupb,zlwb:zupb,:))
if( temp_max > time_max ) then
    time_max = temp_max
    loc = maxloc(deltk(rlwb:rupb,zlwb:zupb,:))
endif

temp_max = maxval(deltl(rlwb:rupb,zlwb:zupb,:))
if( temp_max > time_max ) then
   time_max = temp_max
   loc = maxloc(deltl(rlwb:rupb,zlwb:zupb,:))
endif

!    use the information in loc to prepare a diagnostic message
!    to go to root along with our local max value of 1 / dt

global_diagnostics(iam+1,1) = time_max
global_diagnostics(iam+1,2) = rhf(loc(1)+1)
global_diagnostics(iam+1,3) = dz*(float(loc(2)) - 0.5)
global_diagnostics(iam+1,4) = dphi*(float(loc(3))-1.0)
global_diagnostics(iam+1,5) = deltj(loc(1)+1,loc(2)+1,loc(3))
global_diagnostics(iam+1,6) = deltk(loc(1)+1,loc(2)+1,loc(3))
global_diagnostics(iam+1,7) = deltl(loc(1)+1,loc(2)+1,loc(3))

!   if iam the root pe collect the global diagnotic array, find the max 
!   value of 1 / dt and which pe sent that value.  Then calculate the time
!   step and write out the diagnostic information from the pe that
!   limited the time step.
if( iam_root ) then
   do I = 1, numprocs-1
      call mpi_recv(global_diagnostics(I+1,:),7,REAL_SIZE,I,100+I,          &
                    MPI_COMM_WORLD, istatus, ierror)
   enddo
   index = maxloc(global_diagnostics(:,1))
   if( global_diagnostics(index(1),1) > 0.0 ) then
      dt_visc = 0.25 * sqrt(1.0/global_diagnostics(index(1),1))
   else
      dt_visc = 1.0e8
   endif
!DEBUG   write(30,100) dt_visc, global_diagnostics(index(1),:)
else
   call mpi_send(global_diagnostics(iam+1,:),7,REAL_SIZE,root,100+iam,      &
                 MPI_COMM_WORLD,ierror) 
endif

return
end subroutine visc_time
