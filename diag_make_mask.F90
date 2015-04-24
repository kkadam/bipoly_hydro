!************************************************************************
!* 
!*  DIAG_MAKE_MASK
!*
!************************************************************************
subroutine diag_make_mask(factor,rho_boundary,q,rholoc1,syscom,       &
                          etot1,etot2,ntot1,ntot2,mask,frnum)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!   diag_make_mask creates a mask array that tells us whether a given
!   cell belongs to either star1 or star2 or neither star.  We define
!   "belong" as follows:
!   > first find the surface that marks the boundary between density
!     values above and below rho_boundary
!   > given a paramteric relation for a clear division point between
!     the two stars divide the boundary zones into those that belong
!     to star1 and star2
!   > calculate the average total energy of all cells on the boundary
!     of either star
!   > for each cell on the grid, decide whether it is on star1's or
!     star2's side of the grid and if it is more tightly bound (more
!     negative total energy) than the fiducial boundary zone energy
!     it then belongs to the star in question and is marked
!
!     star2 cells have mask value of -1
!     star1 cells have mask vaule of 1
!     all other cells have vaule 0
!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: rho_diag, potr, vr, vz, vphi
real, dimension(numr_dd,numphi) :: x, y, xin, yin
common /diag/ rho_diag,potr,vr,vz,vphi,x,y,xin,yin

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

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

real, dimension(numr_dd,numz_dd,numphi) :: e

real, dimension(4) :: diag_sum, diag_summed

real, dimension(3) :: divider, side

real :: gammam1inv, test, eavg1, eavg2

integer, dimension(numr_dd,numz_dd,numphi) :: boundary

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
       
real, dimension(3) :: rholoc1, syscom

real :: factor, rho_boundary, q, etot1, etot2, ntot1, ntot2

integer, dimension(numr_dd,numz_dd,numphi) :: mask

integer :: frnum

!*
!************************************************************************
!  initialize the local variables
e = 0.0
diag_sum = 0.0
diag_summed = 0.0
divider(1) = syscom(1)*(1.0-q) + rholoc1(1)*cos(rholoc1(3))*q
divider(2) = syscom(2)*(1.0-q) + rholoc1(1)*sin(rholoc1(3))*q
divider(3) = syscom(3)*(1.0-q) + rholoc1(2)*q
side = 0.0
gammam1inv = 1.0 / (gamma - 1.0)
test = 0.0
eavg1 = 0.0
eavg2 = 0.0
boundary = 0
ierror = 0
      
!-----------------------------------------------------------------------
! mark cells that belong to the surface where rho makes the
! transistion across the rho_boundary value
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb-1
         if( rho_diag(J,K,L) > rho_boundary ) then
            ! J,K,L is inside a *
            if( rho_diag(J+1,K,L) <= rho_boundary ) then
               boundary(J,K,L) = 1
            endif	
         else
            ! J,K,L is outside the two *'s
            if( rho_diag(J+1,K,L) > rho_boundary ) then
               boundary(J+1,K,L) = 1
            endif 
         endif
      enddo
   enddo
enddo

!------------------------------------------------------------------------
! calculate total energy of every cell
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
               e(J,K,L) = factor*rhf(J)*( p(J,K,L)*gammam1inv +    &
                          rho_diag(J,K,L)*(potr(J,K,L) + 0.5*(     &
                          vr(J,K,L)*vr(J,K,L) +                    &
                          vphi(J,K,L)*vphi(J,K,L) +                &
                          vz(J,K,L)*vz(J,K,L) )))
        !DEBUGe(J,K,L) = factor*rhf(J)*rho_diag(J,K,L)*potr(J,K,L)
      enddo
   enddo
enddo 

!------------------------------------------------------------------------
! add up total energy of each boundary cell
! and form the average energy at the bounding
! surface of each star
!
! test is the following dot product:
!
!      (X(system com) - X(division point)) .
!      (X(J,K,L) - X(division point))
!  if test > 0 then on *2 side of the grid
!  if test <= 0 then on *1 side of the grid
!
side = syscom - divider
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( boundary(J,K,L) > 0 ) then
            test = side(1)*(x(J,L) - divider(1)) +                      &
                   side(2)*(y(J,L) - divider(2)) +                      &
                   side(3)*(zhf(K) - divider(3))
            if( test > 0.0 ) then
               ! point lies on *2 side
               etot2 = etot2 + e(J,K,L)
               ntot2 = ntot2 + 1.0
            else
               ! point lies on *1 side
               etot1 = etot1 + e(J,K,L)
               ntot1 = ntot1 + 1.0
            endif
         endif
      enddo
   enddo
enddo
diag_sum(1) = etot1
diag_sum(2) = ntot1
diag_sum(3) = etot2
diag_sum(4) = ntot2  

call mpi_reduce(diag_sum,diag_summed,4,REAL_SIZE,MPI_SUM,root,      &
                MPI_COMM_WORLD,ierror)

if( iam_root ) then
   etot1 = diag_summed(1)
   ntot1 = diag_summed(2)
   etot2 = diag_summed(3)
   ntot2 = diag_summed(4)
   diag_sum(1) = etot1 / ntot1
   diag_sum(2) = etot2 / ntot2
endif
 
call mpi_bcast(diag_sum(1:2),2,REAL_SIZE,root,MPI_COMM_WORLD,ierror)
eavg1 = diag_sum(1)
eavg2 = diag_sum(2)

!------------------------------------------------------------------------
! form the mask array
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
                test = side(1)*(x(J,L) - divider(1)) +               &
                side(2)*(y(J,L) - divider(2)) +                      &
                side(3)*(zhf(K) - divider(3))
         if( test > 0.0 ) then
            ! point lies on *2 side
            if( e(J,K,L) <= eavg2 ) mask(J,K,L) = -1
         else
            ! point lies on *1 side
            if( e(J,K,L) <= eavg1 ) mask(J,K,L) = 1
         endif	
      enddo
   enddo
enddo
 
!------------------------------------------------------------------------      
return
end subroutine diag_make_mask
