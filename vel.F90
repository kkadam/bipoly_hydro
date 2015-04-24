!********************************************************************
!*
!*  VEL
!*
!********************************************************************
subroutine vel
implicit none
include 'runhydro.h'
include 'mpif.h'
!********************************************************************
!*
!   vel converts the conserved kinematic variables into "velocities"
!   for determining the ammount of material advected
!   for the radial direction:
!        u  =  s / <rho>
!
!   for the vertical direction:
!        w  =  t / <rho>
!
!   for the azimuthal direction:
!        jn  =  a / rho, the specific angular momentum
!
!   vel also limits the velocities by the following algorithm:
!
!   --> if the density in the cell is greater then den_cutoff and the
!        velocity is greater than vmax, reset the velocity to vamx
!   --> if the density falls below den_cutoff and the velocity 
!       exceeds vlim_factor times the maximum sound speed then reset
!       the velocity to the velocity limit.  vlim_factor is given
!       a value in setup
!*
!********************************************************************
!*
!*   Global variables
!*

real,dimension(numr_dd,numz_dd,numphi) :: pot, rho 
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) ::  u, w, jn
common /velocities/ u, w, jn

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: pi, grav, cfl_factor, den_cutoff, vlim_factor, viscosity
common /constants/ pi, grav, cfl_factor, den_cutoff, vlim_factor, &
                   viscosity

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
!********************************************************************
!*
!*   Local Variables
!*

real :: rhjkl

real :: vlim, emax, global_emax

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

integer :: J, K, L

!*
!*
!******************************************************************   
!  initialize the local variables 
vlim = 0.0
emax = 0.0
global_emax = 0.0
ierror = 0

!   setup for the velocity limiting:
!   if this is a polytropic evloution use max sound speed to limit
!   the velocity of low density material
emax = maxval(p(rlwb:rupb,zlwb:zupb,:)/rho(rlwb:rupb,zlwb:zupb,:))
call mpi_reduce(emax,global_emax,1,REAL_SIZE,MPI_MAX,root,            &
                MPI_COMM_WORLD,ierror) 
if( iam_root ) then
   vlim = vlim_factor * sqrt(gamma*global_emax)
endif
call mpi_bcast(vlim,1,REAL_SIZE,root,MPI_COMM_WORLD,ierror)


! convert s to u and reset if necessary
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl = 0.5 * rinv(J) *                             &
                 (rho(J,K,L)*(r(J) + 0.25*dr) +              &
                  rho(J-1,K,L)*(r(J) - 0.25*dr) )
         if( rhjkl > 0.0 ) then
	    u(J,K,L) = s(J,K,L) / rhjkl
            if( rhjkl >= den_cutoff ) then
	       if( abs(u(J,K,L)) >= vmax ) then
	          u(J,K,L) = sign(vmax,u(J,K,L))
                  s(J,K,L) = u(J,K,L) * rhjkl
	       endif
	    else
	       if( abs(u(J,K,L)) >= vlim ) then
	          u(J,K,L) = sign(vlim,u(J,K,L))
	          s(J,K,L) = u(J,K,L) * rhjkl
	       endif
            endif
	 else
	    u(J,K,L) = 0.0
	    s(J,K,L) = 0.0
	 endif
      enddo
   enddo
enddo

! convert t to w and reset if necessary
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl = 0.5 * (rho(J,K,L) + rho(J,K-1,L))
         if( rhjkl > 0.0 ) then
            w(J,K,L) = t(J,K,L) / rhjkl
	    if( rhjkl > den_cutoff ) then
	       if( abs(w(J,K,L)) >= 1.0 ) then
	          w(J,K,L) = sign(1.0,w(J,K,L))
		  t(J,K,L) = w(J,K,L) * rhjkl
	       endif
	    else
	       if( abs(w(J,K,L)) >= vlim ) then
	          w(J,K,L) = sign(vlim,w(J,K,L))
	          t(J,K,L) = w(J,K,L) * rhjkl
	       endif
	    endif
	 else
	    w(J,K,L) = 0.0
	    t(J,K,L) = 0.0
	 endif
      enddo
   enddo
enddo

! convert a to jn and reset if necessary
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl = 0.5 * (rho(J,K,L) + rho(J,K,L-1))
         if( rhjkl > 0.0 ) then
	    jn(J,K,L) = a(J,K,L) / rhjkl
	    if( rhjkl > den_cutoff ) then
	       if( abs(jn(J,K,L))*rhfinv(J) >= vmax ) then
	          jn(J,K,L) = rhf(J) * sign(vmax,jn(J,K,L))
		  a(J,K,L) = jn(J,K,L) * rhjkl
	       endif
	    else
	       if( abs(jn(J,K,L))*rhfinv(J) >= vlim ) then
	          jn(J,K,L) = rhf(J) * sign(vlim,jn(J,K,L))
	       a(J,K,L) = jn(J,K,L) * rhjkl
	       endif
	    endif
	 else
	    jn(J,K,L) = 0.0
	    a(J,K,L) = 0.0
	 endif
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      rhjkl = 0.5 * (rho(J,K,philwb) + rho(J,K,phiupb))
      if( rhjkl > 0.0 ) then
         jn(J,K,philwb) = a(J,K,philwb) / rhjkl
         if( rhjkl > den_cutoff ) then
            if( abs(jn(J,K,philwb))*rhfinv(J) >= vmax ) then
               jn(J,K,philwb) = rhf(J) * sign(vmax,jn(J,K,philwb))
               a(J,K,philwb) = jn(J,K,philwb) * rhjkl
            endif
         else
            if( abs(jn(J,K,philwb))*rhfinv(J) >= vlim ) then
               jn(J,K,philwb) = rhf(J) * sign(vlim,jn(J,K,philwb))
               a(J,K,philwb) = jn(J,K,philwb) * rhjkl
            endif
         endif
      else
         jn(J,K,philwb) = 0.0
         a(J,K,philwb) = 0.0
      endif
   enddo
enddo

! apply the boundary conditions
call vel_bc

!   communicate the shared values of s, t, a, u, w, jn that may have been reset
call comm(u)
call comm(w)
call comm(jn)
call comm(a)
call comm(s)
call comm(t)

return
end subroutine vel
