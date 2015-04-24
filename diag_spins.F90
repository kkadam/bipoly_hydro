!************************************************************************
!* 
!*  DIAG_SPINS
!*
!************************************************************************
subroutine diag_spins(factor,star1com,star2com,spin1,spin2,mask)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!
!   diag_spins calculates the z component of spin angular momentum 
!   of each star about their center of mass relative to an inertial
!   frame of reference
!
!   the formula is:
!
!   spin_i = rho ( (x - x_i) vy - (y - y_i) vx) )
!
!          = rho ( (x - x_i) (sin(phi)vr + cos(phi)vphi) -
!
!                  (y - y_i) (cos(phi)vr - sin(phi)vphi) )
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

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                      &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                           &
             in_neighbor, out_neighbor, root,                           &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                             &
           in_neighbor, out_neighbor, root,                             &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                      &
                        iam_on_bottom, iam_on_axis,                     &
                        iam_on_edge, down_neighbor,                     &
                        up_neighbor, in_neighbor,                       &
                        out_neighbor, root, column_num,                 &
                        row_num, pe_grid, iam_root,                     &
                        REAL_SIZE, INT_SIZE

!*
!************************************************************************
!*
!*   Local Variables

real, dimension(numr_dd,numz_dd,numphi) :: temp

real, dimension(2) :: pass, passed

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

integer :: J, K, L

!*
!************************************************************************
!*
!*  Subroutine Arguments

real, dimension(3) :: star1com, star2com

real :: factor, spin1, spin2

integer, dimension(numr_dd,numz_dd,numphi) :: mask

!*
!************************************************************************
!  initialize the local variables
pass = 0.0
passed = 0.0
ierror = 0
 
! add in the frame component of rotation to vphi
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         temp(J,K,L) = vphi(J,K,L) + rhf(J)*omega_frame
      enddo
   enddo
enddo

!------------------------------------------------------------------------------
! sum up the z component of spin angular momentum

spin1 = 0.0
spin2 = 0.0
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( mask(J,K,L) < 0 ) then
            spin2 = spin2 + rhf(J)*rho_diag(J,K,L)*(       &
                            (x(J,L) - star2com(1))*        &
                            (sin_cc(L)*vr(J,K,L) +         &
                             cos_cc(L)*temp(J,K,L)) -      &
                            (y(J,L) - star2com(2))*        &
                            (cos_cc(L)*vr(J,K,L) -         &
                             sin_cc(L)*temp(J,K,L)))
         else if( mask(J,K,L) > 0 ) then
            spin1 = spin1 + rhf(J)*rho_diag(J,K,L)*(       &
                            (x(J,L) - star1com(1))*        &
                            (sin_cc(L)*vr(J,K,L) +         &
                             cos_cc(L)*temp(J,K,L)) -      &
                            (y(J,L) - star1com(2))*        &
                            (cos_cc(L)*vr(J,K,L) -         &
                             sin_cc(L)*temp(J,K,L)))
         endif
      enddo
   enddo
enddo

!------------------------------------------------------------------------------
pass(1) = factor * spin1
pass(2) = factor * spin2

call mpi_reduce(pass,passed,2,REAL_SIZE,MPI_SUM,root,        &
                MPI_COMM_WORLD,ierror)

if( iam_root ) then
   spin1 = passed(1)
   spin2 = passed(2)
endif
       
return
end subroutine diag_spins
