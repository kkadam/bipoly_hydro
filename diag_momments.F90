!************************************************************************
!*
!*  DIAG_MOMMENTS
!*
!************************************************************************
subroutine diag_momments(factor,star1comin,star2comin,mom1,mom2,mask)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!   diag_momments calculates momments of the current density
!   distribution allowing one to construct the projections of
!   the density onto the functions r**l Ylm for 2 <= l <= 4,
!
!   take momments about the center of mass of the system
!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: rho_diag, potr, vr, vz, vphi
real, dimension(numr_dd,numphi) :: x, y, xin, yin
common /diag/ rho_diag,potr,vr,vz,vphi,x,y,xin,yin

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

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

real, dimension(62) :: pass, passed

real :: xl, yl, zl

real :: temp, xx, yy, zz, xxx, yyy, zzz

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

integer :: J, K, L

!*
!************************************************************************
!*
!*   Subroutine Arguments

real, dimension(31) :: mom1, mom2

real, dimension(3) :: star1comin, star2comin

real :: factor

integer, dimension(numr_dd,numz_dd,numphi) :: mask

!*
!************************************************************************
!  initialize the local variables
pass = 0.0
passed = 0.0
xl = 0.0
yl = 0.0
zl = 0.0
temp = 0.0
xx = 0.0
yy = 0.0
zz = 0.0
xxx = 0.0
yyy = 0.0
zzz = 0.0
ierror = 0

do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( mask(J,K,L) < 0 ) then
            xl = xin(J,L) - star2comin(1)
            yl = yin(J,L) - star2comin(2)
            zl = zhf(K) - star2comin(3)
            xx = xl*xl
            yy = yl*yl
            zz = zl*zl
            xxx = xx*xl
            yyy = yy*yl
            zzz = zz*zl
            temp = rhf(J)*rho_diag(J,K,L)
            ! quadratcis
            mom2(1) = mom2(1) + xx*temp
            mom2(2) = mom2(2) + yy*temp
            mom2(3) = mom2(3) + zz*temp
            mom2(4) = mom2(4) + xl*yl*temp
            mom2(5) = mom2(5) + xl*zl*temp
            mom2(6) = mom2(6) + yl*zl*temp
            ! cubics
            mom2(7) = mom2(7) + xxx*temp
            mom2(8) = mom2(8) + yyy*temp
            mom2(9) = mom2(9) + zzz*temp
            mom2(10) = mom2(10) + xx*yl*temp
            mom2(11) = mom2(11) + xx*zl*temp
            mom2(12) = mom2(12) + xl*yy*temp
            mom2(13) = mom2(13) + xl*zz*temp
            mom2(14) = mom2(14) + yy*zl*temp
            mom2(15) = mom2(15) + zz*yl*temp
            mom2(16) = mom2(16) + xl*yl*zl*temp
            ! quartics
            mom2(17) = mom2(17) + xx*xx*temp
            mom2(18) = mom2(18) + yy*yy*temp
            mom2(19) = mom2(19) + zz*zz*temp
            mom2(20) = mom2(20) + xxx*yl*temp
            mom2(21) = mom2(21) + xxx*zl*temp
            mom2(22) = mom2(22) + xx*yy*temp
            mom2(23) = mom2(23) + xx*zz*temp
            mom2(24) = mom2(24) + xx*yl*zl*temp
            mom2(25) = mom2(25) + yyy*xl*temp
            mom2(26) = mom2(26) + yyy*zl*temp
            mom2(27) = mom2(27) + yy*zz*temp
            mom2(28) = mom2(28) + yy*xl*zl*temp
            mom2(29) = mom2(29) + zzz*xl*temp
            mom2(30) = mom2(30) + zzz*yl*temp
            mom2(31) = mom2(31) + zz*xl*yl*temp
         else if( mask(J,K,L) > 0 ) then
            xl = xin(J,L) - star1comin(1)
            yl = yin(J,L) - star1comin(2)
            zl = zhf(K) - star1comin(3)
            xx = xl*xl
            yy = yl*yl
            zz = zl*zl
            xxx = xx*xl
            yyy = yy*yl
            zzz = zz*zl
            temp = rhf(J)*rho_diag(J,K,L)
            ! quadratics
            mom1(1) = mom1(1) + xx*temp
            mom1(2) = mom1(2) + yy*temp
            mom1(3) = mom1(3) + zz*temp
            mom1(4) = mom1(4) + xl*yl*temp
            mom1(5) = mom1(5) + xl*zl*temp
            mom1(6) = mom1(6) + yl*zl*temp
            ! cubics
            mom1(7) = mom1(7) + xxx*temp
            mom1(8) = mom1(8) + yyy*temp
            mom1(9) = mom1(9) + zzz*temp
            mom1(10) = mom1(10) + xx*yl*temp
            mom1(11) = mom1(11) + xx*zl*temp
            mom1(12) = mom1(12) + xl*yy*temp
            mom1(13) = mom1(13) + xl*zz*temp
            mom1(14) = mom1(14) + yy*zl*temp
            mom1(15) = mom1(15) + zz*yl*temp
            mom1(16) = mom1(16) + xl*yl*zl*temp
            ! quartics
            mom1(17) = mom1(17) + xx*xx*temp
            mom1(18) = mom1(18) + yy*yy*temp
            mom1(19) = mom1(19) + zz*zz*temp
            mom1(20) = mom1(20) + xxx*yl*temp
            mom1(21) = mom1(21) + xxx*zl*temp
            mom1(22) = mom1(22) + xx*yy*temp
            mom1(23) = mom1(23) + xx*zz*temp
            mom1(24) = mom1(24) + xx*yl*zl*temp
            mom1(25) = mom1(25) + yyy*xl*temp
            mom1(26) = mom1(26) + yyy*zl*temp
            mom1(27) = mom1(27) + yy*zz*temp
            mom1(28) = mom1(28) + yy*xl*zl*temp
            mom1(29) = mom1(29) + zzz*xl*temp
            mom1(30) = mom1(30) + zzz*yl*temp
            mom1(31) = mom1(31) + zz*xl*yl*temp
         endif	
      enddo
   enddo
enddo

pass(1:31) = factor * mom1
pass(32:62) = factor * mom2

call mpi_reduce(pass,passed,62,REAL_SIZE,MPI_SUM,root,             &
                MPI_COMM_WORLD,ierror)

if( iam_root ) then       
   mom1 = passed(1:31)
   mom2 = passed(32:62)
endif
   
return
end subroutine diag_momments
