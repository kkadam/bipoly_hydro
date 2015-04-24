!***********************************************************************
!*
!*  DEBUGIN
!*
!***********************************************************************
subroutine debugin
implicit none
include 'runhydro.h'
!***********************************************************************
!*
!  debugin sets up simple initial conditions for testing the
!  hydrocode package
!*
!***********************************************************************
!*
!*  Global Variable

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi,                             &
                             drinv, dzinv, dphiinv

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

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
!***********************************************************************
!*
!*  Local Variable

real :: x, y

integer :: J, K, L

!*
!***********************************************************************
!  initialize the local variables
x = 0.0
y = 0.0

!  Sod's problem initial conditions
do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
        ! x = rhf(J) * cos_cc(L)
        ! if( x >= 0.4 ) then
         if( zhf(K) <= 0.5 ) then
            rho(J,K,L) = 1.0
            tau(J,K,L) = 1.92416745
         else
            rho(J,K,L) = 0.125
            tau(J,K,L) = 0.371498572 
         endif              
      enddo
   enddo
enddo
      
s = 0.0
t = 0.0
a = 0.0
pot = 0.0

return
end subroutine debugin
