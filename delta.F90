!********************************************************************
!*
!*  DELTA
!*
!********************************************************************
subroutine delta
implicit none
include 'runhydro.h'
include 'mpif.h'
!********************************************************************
!*
!  This subroutine calculates the global timestep to evolve the 
!  the fluid with.  Get this number by having each processor find 
!  ts maximum value of 
!  
!  dt^-1_i =   (Sound Speed + |V_i|) / dx_i  i = 1, 3
!
!  This value is sent by each processor to root which finds the
!  max of the individual processor values and broadcasts the
!  global timestep back to all processors. root also prints out
!  a diagnostic message every timestep 
!
!  Final timestep is:
!
!  dt  =  cfl_factor  MIN ( dx_i / ( C + |V_i| ) )
!
!  where cfl_factor is set by setup and is typically 0.5
!
!  In calculating the signal velocities, use cell centered radial and
!  vertical velocity approximations
!*
!********************************************************************
!*
!*   Global variables
!*

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi,                          &
        drinv, dzinv, dphiinv

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

integer :: isoadi, call_pot, zero_out
common /flags/ isoadi, call_pot, zero_out

real :: pi, grav, cfl_factor, den_cutoff, vlim_factor, viscosity
common /constants/ pi, grav, cfl_factor, den_cutoff, vlim_factor,  &
                   viscosity

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

real, dimension(numr_dd,numz_dd,numphi) :: c, deltj,  deltk, deltl

real :: adiso, vel, time_max, temp_max, dt_old

real, dimension(numprocs,8) :: global_diagnostics

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
adiso = 0.0
time_max = 0.0
temp_max = 0.0
global_diagnostics = 0.0
loc = 0
istatus = 0
ierror = 0
index = 0

  100  format(10e16.8)

! get the limiting timestep from the diffusive viscous terms
call visc_time

!   form the sound speed
!
!   c^2 = (dp / drho)|s = kappa gamma rho^(gamma-1) = gamma p / rho 
!
!       = gamma (gamma-1) eps rho / rho = gamma (gamma-1) eps
!
adiso = 1.0
if( isoadi == 2 .or. isoadi == 3 ) adiso = gamma
   ! isoadi is a relic from the hpf code that allowed isothermal
   ! or adiabatic equations of state
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         c(J,K,L) = sqrt(adiso*p(J,K,L)/rho(J,K,L))
      enddo
   enddo
enddo

!    deltj = inverse of signal time across cell in radial direction
!
!    deltj = ( c + |<u>| ) / dR
!
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         vel = 0.5*( u(J,K,L) + u(J+1,K,L) )
         deltj(J,K,L) = drinv*(c(J,K,L) + abs(vel))
      enddo
   enddo
enddo 

!    deltk = ( c + |<w>| ) / dZ
!
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         vel = 0.5*( w(J,K,L) + w(J,K+1,L) )
         deltk(J,K,L) = dzinv*(c(J,K,L) + abs(vel))
      enddo
   enddo
enddo 

!    deltl = ( c + |jn / R| ) / (R dphi)
!
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         vel = 0.5 * rhfinv(J) * ( jn(J,K,L) + jn(J,K,L+1) )
         deltl(J,K,L) = rhfinv(J)*dphiinv*( c(J,K,L) + abs(vel) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      vel = 0.5 * rhfinv(J) * ( jn(J,K,phiupb) + jn(J,K,philwb) )
      deltl(J,K,phiupb) = rhfinv(J)*dphiinv*( c(J,K,phiupb) + abs(vel) )
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
global_diagnostics(iam+1,3) = zhf(loc(2)+1)
global_diagnostics(iam+1,4) = dphi*(float(loc(3))-1.0)
global_diagnostics(iam+1,5) = c(loc(1)+1,loc(2)+1,loc(3))
global_diagnostics(iam+1,6) = deltj(loc(1)+1,loc(2)+1,loc(3))*dr
global_diagnostics(iam+1,7) = deltk(loc(1)+1,loc(2)+1,loc(3))*dz
global_diagnostics(iam+1,8) = deltl(loc(1)+1,loc(2)+1,loc(3))*dphi*rhf(loc(1)+1)

!   if iam the root pe collect the global diagnotic array, find the max 
!   value of 1 / dt and which pe sent that value.  Then calculate the time
!   step and write out the diagnostic information from the pe that
!   limited the time step.
if( iam_root ) then
   do I = 1, numprocs-1
      call mpi_recv(global_diagnostics(I+1,:),8,REAL_SIZE,I,100+I,          &
                    MPI_COMM_WORLD, istatus, ierror)
   enddo
   index = maxloc(global_diagnostics(:,1))
   ! save the old timestep
   dt_old = dt
   ! have the Courant limited timestep in dt
   dt = cfl_factor / global_diagnostics(index(1),1)
   ! if the viscous time constraint is more stringent use it
   if( dt > dt_visc ) dt = dt_visc
   ! if the timestep grew by more than 20% limit it to 20% growth
   if( dt > 1.2 * dt_old .and. dt_old > 0.0 ) dt = 1.2 * dt_old
   write(30,100) dt, global_diagnostics(index(1),:)
else
   call mpi_send(global_diagnostics(iam+1,:),8,REAL_SIZE,root,100+iam,      &
                 MPI_COMM_WORLD,ierror) 
endif

!   finally, broadcast the timestep increment
call mpi_bcast(dt,1,REAL_SIZE,root,MPI_COMM_WORLD,ierror)  

return
end subroutine delta
