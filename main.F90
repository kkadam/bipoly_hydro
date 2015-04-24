!************************************************************
!*
!*  MAIN
!*
!************************************************************
program main
implicit none
include 'runhydro.h'
include 'mpif.h'
include 'main_formats.h'
!************************************************************
!*
!*  Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real :: dt, time, dt_visc
integer :: tstep 
common /timestep/ dt, time, dt_visc, tstep

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi,                   &    
                             drinv, dzinv, dphiinv

integer :: isoadi, call_pot, zero_out
common /flags/ isoadi, call_pot, zero_out

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

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
!************************************************************      
!*
!*   Local variables

real, dimension(numr_dd-2,numz_dd-2,numphi) :: dummy

integer, dimension(numr_procs*numz_procs,2) :: pe_coord

integer :: tstart, tstop, do_diag

integer :: one, two, conts_tag, conts_tstep

integer :: I, J, N, frnum, frnum_init

real :: intrvl, rho_boundary, q

real :: time1, time2, timef

real :: drag_factor

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

integer :: record_length

character(len=50) :: template, filename

!*
!************************************************************      
!  Initialize MPI
call mpi_init(ierror)

! Initialize all data in common blocks
call initialize

! set the default data sizes for MPI communications
#ifdef SHORT
REAL_SIZE = MPI_REAL4
INT_SIZE = MPI_INTEGER4
#else
REAL_SIZE = MPI_REAL8
INT_SIZE = MPI_INTEGER
#endif

! TIMING
time1 = mpi_wtime()

! Initialize local variables
pe_coord = 0
tstart = 0
tstop = 0
do_diag = 1
one = 1
two = 2
conts_tag = 1
conts_tstep = 5000
ierror = 0
frnum = 0
frnum_init = 0
intrvl = 0.0
ierror = 0

open(unit=6, file='output/drive.out', &
     form='formatted',status='unknown', position='append')

!  Set up logical and integer variables that describe the
!  processor grid and the message passing pattern

iam_on_top = .false.
iam_on_bottom = .false.
iam_on_axis = .false.
iam_on_edge = .false.
iam_root = .false.
root = 0

call mpi_comm_rank(MPI_COMM_WORLD, iam, ierror)
numprocs = numr_procs * numz_procs
if( iam == 0 ) iam_root = .true.

!  Make sure the input data in runhdro.h doesn't violate
!  assumptions made in the code implementaion.  If there
!  is a problem might as well stop now...
if( mod(numr_procs,2) /= 0 ) then
   ! not en even number of processors in radial
   ! direction to avoid deadlock in message passing
   if( iam_root ) write(6,100)
   stop
endif

if( mod(numz_procs,2) /= 0 ) then
   ! not an even number of processors in vertical
   ! direction to avoid deadlock in message passing
   if( iam_root ) write(6,110)
   stop
endif

if( mod(numphi,numr_procs) /= 0 ) then
   ! to do data swapping in helmadi the azimuthal
   ! dimension of the data has to distribute evenly
   ! across numr_procs
   if( iam_root ) write(6,120)
   stop
endif

if( mod(numr-2,numz_procs) /= 0 ) then
   ! helmadi's communications also require that
   ! radial data be divisible
   ! evenly by numz_procs.  Checking if I picked
   ! a valid pair of numr and numz_procs
   if( iam_root ) write(6,130)
   stop
endif 

!  Setup the logical pe grid
N = 1
do J = 1, numz_procs
   do I = 1, numr_procs
      pe_coord(N,1) = I-1
      pe_coord(N,2) = J-1
      pe_grid(I,J) = N-1
      N = N + 1
   enddo
enddo

!  Determine which processors have external boundary zones
!  in place of a guard cell layer that will require special
!  treatment
if( pe_coord(iam+1,2) == 0 ) iam_on_bottom = .true.

if( pe_coord(iam+1,2) == numz_procs-1 ) iam_on_top = .true.

if( pe_coord(iam+1,1) == 0 ) iam_on_axis = .true.

if( pe_coord(iam+1,1) == numr_procs-1 ) iam_on_edge = .true. 

!  Determine passing partners for comm and related routines
if( iam_on_top ) then
   up_neighbor = -1
else
   up_neighbor = iam + numr_procs
endif

if( iam_on_bottom ) then
   down_neighbor = -1
else
   down_neighbor = iam - numr_procs
endif

if( iam_on_axis ) then
   in_neighbor = -1
else
   in_neighbor = iam - 1
endif

if( iam_on_edge ) then
   out_neighbor = -1
else
   out_neighbor = iam + 1
endif

!  Determine whether pe is in an even row or column to order
!  message passing
row_num = pe_coord(iam+1,2)
 
column_num = pe_coord(iam+1,1)

call setup(tstart, tstop, do_diag, intrvl, frnum, rho_boundary, q)

! the timestep file holds the diagnostic output from delta
open(unit=30, file='output/timestep', &
     form='formatted',status='unknown', position='append')

frnum_init = frnum
inquire(iolength=record_length) dummy

do I = 1, 8

   if ( I == 1 ) template = 'output/data/frame'
   if ( I == 2 ) template = 'output/data/frac1.'
   if ( I == 3 ) template = 'output/data/frac2.'
   if ( I == 4 ) template = 'output/data/tau'
   if ( I == 5 ) template = 'output/data/pot'
   if ( I == 6 ) template = 'output/data/velr'
   if ( I == 7 ) template = 'output/data/velz'
   if ( I == 8 ) template = 'output/data/velphi'

   if ( iam < 10 ) then
      write(filename,'(a,i4,a,i1)') trim(template), frnum, '_', iam
   else if ( iam >= 10 .and. iam < 100 ) then
      write(filename,'(a,i4,a,i2)') trim(template), frnum, '_', iam
   else if ( iam >= 100 .and. iam < 1000 ) then
      write(filename,'(a,i4,a,i3)') trim(template), frnum, '_', iam
   else if ( iam >= 1000 .and. iam < 10000 ) then
      write(filename,'(a,i4,a,i4)') trim(template), frnum, '_', iam
   endif

   if ( I == 1 ) then
      open(unit=50,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   else if ( I == 2 ) then
      open(unit=51,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   else if ( I == 3 ) then
      open(unit=52,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   else if ( I == 4 ) then
      open(unit=53,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   else if ( I == 5 ) then
      open(unit=54,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   else if ( I == 6 ) then
      open(unit=55,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   else if ( I == 7 ) then
      open(unit=56,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   else if ( I == 8 ) then
      open(unit=57,file=trim(filename),form='unformatted',status='new', &
           access='direct',recl=record_length)
   endif

enddo

if( tstart == 1 ) then
   ! Initialize delt which holds the timestep size from the previous
   ! couple of timesteps. Note that this has to be done only for the
   ! initial run i.e. tstart = 1
   delt = 0.0

   ! if this is the very first timestep in a simulation
   ! we have to initialize the temps and tempa arrays
   ! for sourcing s and a
   call delta
   dt = 0.5 * dt
   call make_source_temp

! 05/24/2003 :
! Calculate the com coordinates (first call)
   call calc_com_coords

! 05/24/2003 :
! Calculate the velocity and acceleration of the com (first call)
   call com_vel_accn

   ! and also, calculate diagnostics for the initial model
   call dodiag(rho_boundary,q,frnum)
   ! and dump the initial data arrays for imaging
   if( iam_root ) then
      write(6,*) 'frnum = ',frnum,' time = ',time
   endif
   call output_files(frnum,frnum_init)
   frnum = frnum + 1
endif

!TIMING
time2 = mpi_wtime()
if( iam_root ) then
   write(6,*) 'Startup Time ',(time2-time1)*1.0e-3
endif

! Set the ammount of drag of angular momentum per timestep
! with drag_factor = 0.01 the system will lose 1% of its 
! angular momentum per orbit
!drag_factor = 1.0e-2
drag_factor = reallyadrag


!TIMING 
time1 = mpi_wtime()

do tstep = tstart, tstop          ! Start timestep loop

   call delta

   dt = 0.5 * dt

   call source(one)

   call vel

   if( iam_on_axis .and. (zero_out > 0) ) then
      s(1:zero_out,:,:) = 0.0
      u(1:zero_out,:,:) = 0.0
      t(1:zero_out,:,:) = 0.0
      w(1:zero_out,:,:) = 0.0
      a(1:zero_out,:,:) = 0.0
      jn(1:zero_out,:,:) = 0.0
   endif

   call save_state

   call advect(one)

   call vel

   dt = 2.0*dt

   call advect(two)

   call visc

   if( iam_on_axis .and. (zero_out > 0) ) then
      s(1:zero_out,:,:) = 0.0
      u(1:zero_out,:,:) = 0.0
      t(1:zero_out,:,:) = 0.0
      w(1:zero_out,:,:) = 0.0
      a(1:zero_out,:,:) = 0.0
      jn(1:zero_out,:,:) = 0.0
      rho(1:zero_out,:,:) = densmin
      tau(1:zero_out,:,:) = taumin
   endif

   call state

   if( call_pot == 1 ) then
      call potential_solver
   endif

   dt = 0.5*dt

   call source(two)

   call vel

   if( iam_on_axis .and. (zero_out > 0) ) then
      s(1:zero_out,:,:) = 0.0
      u(1:zero_out,:,:) = 0.0
      t(1:zero_out,:,:) = 0.0
      w(1:zero_out,:,:) = 0.0
      a(1:zero_out,:,:) = 0.0
      jn(1:zero_out,:,:) = 0.0
   endif

   dt = 2.0*dt

   if ( time / cirp >=  dragtime ) then
      drag_factor = 0.0
   endif

   call drag(drag_factor,dt)

! 05/24/2003 :
! Calculate the com
   call calc_com_coords

! 05/24/2003 :
! Calculate the velocity and acceleration of the com
   call com_vel_accn

   if( mod(tstep,do_diag) == 0 ) call dodiag(rho_boundary,q,frnum)

   time = time + dt

   if( intrvl*time/cirp > (frnum-1000) ) then  !ready to write out rho
       if( iam_root ) then
           write(6,*) 'frnum = ',frnum,' time = ',time
       endif
       call output_files(frnum,frnum_init)
       frnum = frnum + 1
    endif

    if( mod(tstep,conts_tstep) == 0 ) then
       call ritecont_frac_new(frnum)
       conts_tag = conts_tag + 1
    endif

enddo                             !end of timestep loop

!TIMING
time2 = mpi_wtime()
if( iam_root ) then
   write(6,*) 'Hydrocode loop time ',(time2-time1)
endif

if( tstop > 100 ) then
   call ritecont(frnum)
endif

! close the data output files
close(50)
close(51)
close(52)
close(53)
close(54)
close(55)
close(56)
close(57)

!  Clean up after MPI
call mpi_finalize(ierror)

stop
end program main
