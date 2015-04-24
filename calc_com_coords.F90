!************************************************************************
!*
!* CALC_COM_COORDS
!*
!************************************************************************
subroutine calc_com_coords
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
! Calculates the center of mass of the system in the rotating frame 
! (cartesian coordinates x, y, z and cylindrical radius R). Also stores 
! the com and timestep sizes for the previous couple of timesteps.
!*
!************************************************************************
!*
!*   Global variables :

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd,numz_dd,numphi) :: rho_diag, potr, vr, vz, vphi
real, dimension(numr_dd,numphi) :: x, y, xin, yin
common /diag/ rho_diag, potr, vr, vz, vphi, x, y, xin, yin

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

real :: pi, grav, cfl_factor, den_cutoff, vlim_factor, viscosity
common /constants/ pi, grav, cfl_factor, den_cutoff, vlim_factor,  &
                   viscosity

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

logical :: iam_on_top, iam_on_bottom, iam_on_axis,               &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                    &
             in_neighbor, out_neighbor, root,                    &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                      &
           in_neighbor, out_neighbor, root,                      &
           REAL_SIZE, INT_SIZE, numprocs
#endif

integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,               &
                        iam_on_bottom, iam_on_axis,              &
                        iam_on_edge, down_neighbor,              &
                        up_neighbor, in_neighbor,                &
                        out_neighbor, root, column_num,          &
                        row_num, pe_grid, iam_root,              &
                        REAL_SIZE, INT_SIZE

!*
!************************************************************************
!*
!*   Local Variables

real, dimension(4) :: com_sum, syscom_sum

real :: factor, mtot

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

integer :: J, K, L

!*
!************************************************************************
!  initialize the local variables
com_sum = 0.0
syscom_sum = 0.0
mtot = 0.0
ierror = 0

!------------------------------------------------------------------------
! factor is the common volume multiplier for global sums
! modify the value to include effects due to symmetries
factor = dr * dz * dphi
if( isym == 2 ) then
   factor = 2.0 * factor
else if( isym == 3 ) then
   factor = 4.0 * factor
endif

!------------------------------------------------------------------------
! x and y are the cartesian coordinates of a cell center in the
! rotating frame
do L = philwb, phiupb
   do J = rlwb-1, rupb+1
      x(J,L) = rhf(J) * cos_cc(L)
      y(J,L) = rhf(J) * sin_cc(L)
   enddo
enddo

!------------------------------------------------------------------------
! form rho_diag, the density array with boundary zones 
! zeroed out so we don't get effects from material that
! has been advected off the grid
rho_diag = rho
if( iam_on_bottom ) rho_diag(:,zlwb-1:zlwb,:) = 0.0
if( iam_on_top )    rho_diag(:,zupb:zupb+1,:) = 0.0
if( iam_on_edge )   rho_diag(rupb:rupb+1,:,:) = 0.0

!-----------------------------------------------------------------------
! Store the center of mass of the total system for the previous couple
! of timesteps
do K = 1, 3
   do J = 3, 2, -1
      com(J,K) = com(J-1,K)
   enddo
   com(1,K) = 0.0
enddo

do J = 3, 2, -1
   delt(J) = delt(J-1)
enddo

!------------------------------------------------------------------------
! find the center of mass for the total system
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         com(1,1) = com(1,1) + rhf(J)*x(J,L)*rho_diag(J,K,L)
         com(1,2) = com(1,2) + rhf(J)*y(J,L)*rho_diag(J,K,L)
         com(1,3) = com(1,3) + rhf(J)*zhf(K)*rho_diag(J,K,L)
         mtot = mtot + rhf(J)*rho_diag(J,K,L)
      enddo
   enddo
enddo
com_sum(1) = factor*com(1,1)
com_sum(2) = factor*com(1,2)
com_sum(3) = factor*com(1,3)
com_sum(4) = factor*mtot

call mpi_reduce(com_sum(1:4),syscom_sum(1:4),4,REAL_SIZE,          &
                MPI_SUM,root,MPI_COMM_WORLD,ierror)

if( iam_root ) then
   mtot = syscom_sum(4)
   com(1,1) = syscom_sum(1)/mtot
   com(1,2) = syscom_sum(2)/mtot
   com(1,3) = syscom_sum(3)/mtot
endif

call mpi_bcast(com,9,REAL_SIZE,root,MPI_COMM_WORLD,ierror)

! Calculates the cylindrical radius of the com and its inverse
R_com = sqrt(com(1,1)**2 + com(1,2)**2)
R_com_inv = 1.0/R_com

! Calculate the angle the com makes with the x axis.
! Done only for the current timestep.
phi_com = atan(com(1,2)/com(1,1))

! Store the timestep size in delt
delt(1) = dt

return
end subroutine calc_com_coords
