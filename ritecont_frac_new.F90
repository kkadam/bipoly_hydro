!************************************************************************
!*
!*  RITECONT_FRAC_NEW
!*
!************************************************************************
subroutine ritecont_frac_new(frnum)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!  ritecont writes out the state of the hydrocode to be used as
!  initial conditions for the next run.  All pes dump their entire
!  portion of the array (guard and live cells).  Arrays go into
!  file fort.12, the scalars for the simulation (time and current
!  frame number) go into fort.13
!
!  modified 5/17/2000 to write out the mass fraction array in
!  addition to the other 8 fluid arrays.
!
!  10/02/2002 : Modified the write statement so that each processor now
!               writes a separate continuation file. Previously each
!               processor wrote to a common continuation file (fort.12).
!               The write statements for fort.13 however is unchanged
!               since we only want root to write this file.
!               NOTE : Supermike needs files in little_endian format.
!*  
!************************************************************************
!*
!*  Subroutine Arguments

integer :: frnum

!*
!************************************************************************
!*
!*   Global Variables

real, dimension(numr_dd, numz_dd, numphi) :: pot, rho
common /poisson/    pot, rho

real, dimension(numr_dd, numz_dd, numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

real, dimension(numr_dd, numz_dd, numphi) :: p, tau
common /thermo/  p, tau

real, dimension(numr_dd, numz_dd, numphi) :: s, t, a
common /kinematic/     s, t, a 

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd, numz_dd, numphi) :: temps, tempa
common /source_temp/ temps, tempa

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep 

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
!*  Local variables

integer :: record_length, token

#ifdef SHORT
integer*8 :: token_tag, ierror

integer*8, dimension(MPI_STATUS_SIZE) :: istatus
#else
integer :: token_tag, ierror

integer, dimension(MPI_STATUS_SIZE) :: istatus
#endif

character :: dot
character(len=150) :: template12, template13, conts12, conts13
character(len=150) :: conts12_file, conts13_file

!*
!*************************************************************************
!  initialize the local variable
dot = '.'
record_length = 0
token = 1
token_tag = 100
ierror = 0
istatus = 0

template12 = 'output/conts/fort.12.'
template13 = 'output/conts/fort.13.'

if (tstep .lt. 100000) then
    write(conts12, '(a,i5)') trim(template12), tstep
    write(conts13_file, '(a,i5)') trim(template13), tstep
elseif ( (tstep .ge. 100000) .and. (tstep .lt. 1000000) ) then
    write(conts12, '(a,i6)') trim(template12), tstep
    write(conts13_file, '(a,i6)') trim(template13), tstep
else
    write(conts12, '(a,i7)') trim(template12), tstep
    write(conts13_file, '(a,i7)') trim(template13), tstep
endif

write(conts12, '(a,a1)') trim(conts12), dot

  110 format(1h1,///,' sssssssssssssssssssssssssssssssssssssssssss',/,  &
      ' s',41x,'s',/,' s  model from time step number',i6,' has   s',/, &
      ' s  been stored on disk.  ',/,                                   &
      ' sssssssssssssssssssssssssssssssssssssssssss',/////)             

inquire(iolength=record_length) rho

! 10/02/2002 : Modified the write statements for fort.12

! create the filenames for the files every processor is going to write out
if (iam .lt. 10) then
    write(conts12_file,'(a,i1)') trim(conts12),iam
elseif ( (iam .ge. 10) .and. (iam .le. 99) ) then
    write(conts12_file,'(a,i2)') trim(conts12),iam
else
    write(conts12_file,'(a,i3)') trim(conts12),iam
endif

open(unit=12,file=trim(conts12_file),form='unformatted', status='unknown')
!open(unit=12,file=trim(conts12_file),form='unformatted')

write(12) s, t, a, rho, tau, pot, tempa, frac1, frac2

if( iam_root ) then
   open(unit=13, file=trim(conts13_file), form='unformatted', status='unknown')

   write(13) time, frnum, phi_com, R_com, com, v_com, a_com, cylin_a_com,  &
             cylin_v_com, delt
   close(13)
   write(6,110) tstep
endif

call mpi_barrier(MPI_COMM_WORLD,ierror)

return
end subroutine ritecont_frac_new
