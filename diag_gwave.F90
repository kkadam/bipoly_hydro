!************************************************************************
!*
!*  DIAG_GWAVE
!*
!************************************************************************
subroutine diag_gwave(factor,ilm)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!   diag_gwave calculates the 6 nontrivial components of the second
!   time derivative of the moment of inertia for the binary system.
!   Use the trick of introducing the equation of continuity and
!   Euler's equations to eliminate the need to take numerical time
!   derivatives.  The following god-awful formula is from Rasio &
!   Shapiro, 1992, ApJ, vol 401, page 229
!
!   (dt)^2 I_lm = Integral ( 2 rho v_l v_m + 2 P delta_lm 
!
!               - 2 ( x_m (d_l) Phi + x_l (d_m) Phi ) ) dV
!
!    Use standard coordinate transformations to express this
!    in cyclindrical coordinates.  Once an observation point
!    is chosen the strain components (in transverse traceless 
!    guage) h+ and hx can be related to the calculated (dt)^2I_lm.
!
!    This calculation is done in "inertial" space, that is to
!    say the azimuthal coordinate is advanced by 2 Pi times
!    the number of orbits the frame has moved thus far in the
!    calculation and the azimuthal veloicty has the frame's
!    motion added back into it.
!
!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: rho_diag, potr, vr, vz, vphi
real, dimension(numr_dd,numphi) :: x, y, xin, yin
common /diag/ rho_diag, potr, vr, vz, vphi, x, y, xin, yin

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

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

real, dimension(numr_dd,numz_dd,numphi) :: vtemp

real, dimension(numr_dd,numphi) :: costemp, sintemp

real, dimension(6) :: ilm_summed

real :: ixx, iyy, izz, ixy, ixz, iyz

real :: dpotdr, dpotdz, dpotdphi

real :: vr2, vphi2, vrvphi, cos2phi

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

real, dimension(6) :: ilm

real :: factor

!*
!************************************************************************
!  initialize the local variables
vtemp = 0.0
costemp = 0.0
sintemp = 0.0
ilm_summed = 0.0
ixx = 0.0
iyy = 0.0
izz = 0.0
ixy = 0.0
ixz = 0.0
iyz = 0.0
dpotdr = 0.0
dpotdz = 0.0
dpotdphi = 0.0
vr2 = 0.0
vphi2 = 0.0
vrvphi = 0.0
cos2phi = 0.0
ierror = 0

!----------------------------------------------------------------------- 
! set up the cosine and sine arrays for the inertial space
do L = philwb, phiupb
   do J = rlwb, rupb
      costemp(J,L) = xin(J,L) * rhfinv(J)
      sintemp(J,L) = yin(J,L) * rhfinv(J)
   enddo
enddo

!----------------------------------------------------------------------- 
! add the frame motion into the azimuthal velocity
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         vtemp(J,K,L) = vphi(J,K,L) + rhf(J) * omega_frame
      enddo
   enddo
enddo

!----------------------------------------------------------------------- 
! do the six integrals
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         cos2phi = costemp(J,L)*costemp(J,L) -                         &
                   sintemp(J,L)*sintemp(J,L)
         dpotdr = 0.5*drinv*(pot(J+1,K,L) - pot(J-1,K,L))
         dpotdz = 0.5*dzinv*(pot(J,K+1,L) - pot(J,K-1,L))
         if( L == philwb ) then
            dpotdphi = 0.5*dphiinv*(pot(J,K,philwb+1) -                &
                                    pot(J,K,phiupb))
         else if( L == phiupb ) then
            dpotdphi = 0.5*dphiinv*(pot(J,K,philwb) -                  &
                                    pot(J,K,phiupb-1))
         else
            dpotdphi = 0.5*dphiinv*(pot(J,K,L+1) - pot(J,K,L-1))
         endif
         vr2 = vr(J,K,L)*vr(J,K,L)
         vphi2 = vtemp(J,K,L)*vtemp(J,K,L)
         vrvphi = vr(J,K,L)*vtemp(J,K,L)

         ixx = ixx + rhf(J) * ( p(J,K,L) + rho_diag(J,K,L)*            &
         (costemp(J,L)*(costemp(J,L)*vr2 -                             &
                        2.0*sintemp(J,L)*vrvphi -                      &
                        rhf(J)*costemp(J,L)*dpotdr +                   &
                        sintemp(J,L)*dpotdphi) +                       &
          sintemp(J,L)*sintemp(J,L)*vphi2) )

         iyy = iyy + rhf(J) * ( p(J,K,L) + rho_diag(J,K,L)*            &
         (sintemp(J,L)*(sintemp(J,L)*vr2 +                             &
                        2.0*costemp(J,L)*vrvphi -                      &
                        rhf(J)*sintemp(J,L)*dpotdr -                   &
                        costemp(J,L)*dpotdphi) +                       &
          costemp(J,L)*costemp(J,L)*vphi2) )

         izz = izz + rhf(J) * ( p(J,K,L) + rho_diag(J,K,L)*            &
         ( vz(J,K,L)*vz(J,K,L) - zhf(K)*dpotdz ) )

         ixy = ixy + rhf(J) * rho_diag(J,K,L) *                        &
         (2.0*costemp(J,L)*sintemp(J,L)*                               &
          ( vr2 - vphi2 - rhf(J)*dpotdr ) +                            &
          cos2phi*(2.0*vrvphi - dpotdphi) )

         ixz = ixz + rhf(J) * rho_diag(J,K,L) *                        &
         ( 2.0*vz(J,K,L)*(costemp(J,L)*vr(J,K,L) -                     &
                          sintemp(J,L)*vtemp(J,K,L)) -                 &
           costemp(J,L)*( rhf(J)*dpotdz + zhf(K)*dpotdr ) +            &
           zhf(K)*rhfinv(J)*sintemp(J,L)*dpotdphi )

         iyz = iyz + rhf(J) * rho_diag(J,K,L) *                        &
         ( 2.0*vz(J,K,L)*(sintemp(J,L)*vr(J,K,L) +                     &
                          costemp(J,L)*vtemp(J,K,L)) -                 &
          sintemp(J,L)*( rhf(J)*dpotdz + zhf(K)*dpotdr ) -             &
          zhf(K)*rhfinv(J)*costemp(J,L)*dpotdphi )

      enddo
   enddo
enddo

ilm(1) = 2.0 * factor * ixx
ilm(2) = 2.0 * factor * iyy
ilm(3) = 2.0 * factor * izz
ilm(4) = factor * ixy
ilm(5) = factor * ixz
ilm(6) = factor * iyz

!----------------------------------------------------------------------- 
! do the gloabl summation of I_lm
call mpi_reduce(ilm,ilm_summed,6,REAL_SIZE,MPI_SUM,root,   &
                MPI_COMM_WORLD,ierror)

if( iam_root ) then
   ilm = ilm_summed
endif

!----------------------------------------------------------------------- 
return
end subroutine diag_gwave
