!************************************************************************
!*
!*  DODIAG.1
!*
!************************************************************************
subroutine dodiag(rho_boundary,q,frnum)
implicit none
include 'mpif.h'
include 'runhydro.h'
!************************************************************************
!*
!   ritediag calculates a whole boatload of quantities of interest when
!   simulating a binary system.  These quantities include (in order they
!   are computed):
!
!   >  (x,y,z) coordinates of system center of mass
!   >  (R,z,phi) coordinates of max density for each star and
!      the corresponding density values
!   >  (x,y,z) coordinates of L1, L2, L3, outer edge of Roche
!      lobe for each of these points and the corresponding
!      effective potential values at each of these positions
!   >  total volume occupied by each Roche lobe
!   >  total mass of material in star1 (star with largest density
!      maxima).  The material that makes up *1 is determined by
!      whether the fluid has lower energy than the average energy
!      of the fluid that makes up the boundary of the star at a
!      density level of rho_boundary
!   >  total mass of companion star, total mass of envelope material
!      and total mass in the boundary zones
!   >  volume occupied by material in each star that is more
!      dense than the five cutoff values rho[1-5]
!   >  center of mass coordinates for each star 
!   >  three principal momments of inertia in the corotating frame
!   >  gravitational potential energy of the system
!   >  total angular momentum of the system in the inertial frame
!   >  kinetic energy of rotation in the inertial frame
!   >  kinetic energy of motion in radial and vertical directions
!   >  integral of pressure for the virial 
!   >  total mass on the grid, including boundary zones as a check
!   >  coordinates of system's center of mass and each star's center
!      of mass relative to the inertial frame of reference
!   >  z component of spin angular momentum for each star relative
!      to the inertial frame
!   >  integrals of all polynomials in Cartessian coordinates of 
!      order two thruough four times the density field for 
!      calculating the projection of the density field onto the
!      r**l Ylm basis functions.  Done relative to the instantaneous
!      system center of mass
!   >  second time derivatives of Ilm for computing the gravitational
!      strain in Traceless Transverse gauge (h+ and hx).  Done in
!      inertial space.
!
!  Modified 5/18/2000 to add the following functionality:
!
!  => Calculate the total angular momentum about the instantaneous
!     center of mass - specific to the material in each star.
!
!  => Calculate the total mass of fluid component one and two to the
!     left and right of the plane of separation.
!
!  => Calculate total mass to left and right f plane of separation for
!     material with density greater than 1.0e-7.
!*
!************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: rho_diag, potr, vr, vz, vphi
real, dimension(numr_dd,numphi) :: x, y, xin, yin
common /diag/ rho_diag, potr, vr, vz, vphi, x, y, xin, yin

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

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

real, dimension(numr_dd,numz_dd,numphi) :: temp

real, dimension(31) :: mom1, mom2

real, dimension(36) :: diag_sum, diag_summed

real, dimension(4,5) :: lpoints

real, dimension(6) :: ilm

real, dimension(5) :: vol1, vol2

real, dimension(3) :: rholoc1, rholoc2, syscom, star1com,         &
                      star2com, syscomin, star1comin, star2comin, &
                      side, divider
 
real :: factor, rhomax1, rhomax2, sum1, sum2, sum3

real :: roche_vol1, roche_vol2, spin_star1, spin_star2

real :: etot1, etot2, ntot1, ntot2, costemp, sintemp

real :: mass_star1, mass_star2, mass_envelope

real :: mass_boundary, mtot, mtot_check, xx, yy, zz

real :: egrav, total_j, trot, trz, pres, phase

real :: m1_left, m2_left, m1_right, m2_right

real :: mrho_left, mrho_right, binary_j, test

integer, dimension(numr_dd,numz_dd,numphi) :: mask

integer, dimension(3) :: irholoc1, irholoc2

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

integer :: J, K, L

real, parameter :: rho1 = 1.0e-1, rho2 = 1.0e-2, rho3 = 1.0e-3
real, parameter :: rho4 = 1.0e-4, rho5 = 1.0e-5

!*
!************************************************************************
!*
!*  Subroutine Arguments

real :: rho_boundary, q

integer :: frnum

!*
!************************************************************************
!  initialize the local variables
temp = 0.0
diag_sum = 0.0
diag_summed = 0.0
lpoints = 0.0
mom1 = 0.0
mom2 = 0.0
ilm = 0.0
vol1 = 0.0
vol2 = 0.0
rholoc1 = 0.0
rholoc2 = 0.0
syscom = 0.0
star1com = 0.0
star2com = 0.0
syscomin = 0.0
star1comin = 0.0
star2comin = 0.0
side = 0.0
divider = 0.0
factor = 0.0
rhomax1 = 0.0
rhomax2 = 0.0
sum1 = 0.0
sum2 = 0.0
sum3 = 0.0
roche_vol1 = 0.0
roche_vol2 = 0.0
spin_star1 = 0.0
spin_star2 = 0.0 
etot1 = 0.0
ntot1 = 0.0
etot2 = 0.0
ntot2 = 0.0
costemp = 0.0
sintemp = 0.0
mass_star1 = 0.0
mass_star2 = 0.0
mass_envelope = 0.0
mass_boundary = 0.0
mtot = 0.0
mtot_check = 0.0
xx = 0.0
yy = 0.0
zz = 0.0
egrav = 0.0
total_j = 0.0
trot = 0.0
trz = 0.0
pres = 0.0
phase = 0.0
m1_left = 0.0
m2_left = 0.0
m1_right = 0.0
m2_right = 0.0
mrho_left = 0.0
mrho_right = 0.0
binary_j = 0.0
test = 0.0
mask = 0
irholoc1 = 0
irholoc2 = 0
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
! x and y are the cartessian coordiantes of a cell center
do L = philwb, phiupb
   do J = rlwb-1, rupb+1
      x(J,L) = rhf(J) * cos_cc(L)
      y(J,L) = rhf(J) * sin_cc(L)
   enddo
enddo 

!------------------------------------------------------------------------
! xin and yin are the cartessian coordinates of a cell center
! in the inertial frame
phase = 2.0 * pi * time / cirp
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         xin(J,L) = rhf(J) * cos( phi(L) + phase )
         yin(J,L) = rhf(J) * sin( phi(L) + phase ) 
      enddo
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

!------------------------------------------------------------------------
! form the Roche potential
do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         potr(J,K,L) = pot(J,K,L) - 0.5*omega_frame*omega_frame*       &
                                    rhf(J)*rhf(J)
      enddo
   enddo
enddo

!------------------------------------------------------------------------
! form vphi, the cell centered azimuthal velocity, relative
! to the corotating frame
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         vphi(J,K,L) = 0.5*rhfinv(J)*(jn(J,K,L) + jn(J,K,L+1))
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      vphi(J,K,phiupb) = 0.5*rhfinv(J)*(jn(J,K,phiupb) + jn(J,K,philwb))
   enddo
enddo

!------------------------------------------------------------------------
! form vr, the cell centered radial velocity
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         vr(J,K,L) = 0.5*(u(J,K,L)+ u(J+1,K,L))
      enddo
   enddo
enddo

!------------------------------------------------------------------------
! form vz, the cell centered vertical velocity
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         vz(J,K,L) = 0.5*(w(J,K,L) + w(J,K+1,L))
      enddo
   enddo
enddo

!------------------------------------------------------------------------
! find the center of mass for the total system
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         syscom(1) = syscom(1) + rhf(J)*x(J,L)*rho_diag(J,K,L)
         syscom(2) = syscom(2) + rhf(J)*y(J,L)*rho_diag(J,K,L)
         syscom(3) = syscom(3) + rhf(J)*zhf(K)*rho_diag(J,K,L)
         mtot = mtot + rhf(J)*rho_diag(J,K,L)
      enddo
   enddo
enddo
diag_sum(1) = factor*syscom(1)
diag_sum(2) = factor*syscom(2)
diag_sum(3) = factor*syscom(3)
diag_sum(4) = factor*mtot

call mpi_reduce(diag_sum(1:4),diag_summed(1:4),4,REAL_SIZE,          &
                MPI_SUM,root,MPI_COMM_WORLD,ierror)

if( iam_root ) then
   mtot = diag_summed(4)
   syscom(1) = diag_summed(1)/mtot
   syscom(2) = diag_summed(2)/mtot
   syscom(3) = diag_summed(3)/mtot
endif

call mpi_bcast(syscom,3,REAL_SIZE,root,MPI_COMM_WORLD,ierror)   

!------------------------------------------------------------------------------
! calculate the system center of mass coordinates in the inertial frame
costemp = cos(phase)
sintemp = sin(phase)

syscomin(1) = syscom(1)*costemp - syscom(2)*sintemp
syscomin(2) = syscom(1)*sintemp + syscom(2)*costemp
syscomin(3) = syscom(3)

!------------------------------------------------------------------------
! find the global maximum density and the maximum density of
! the companion star
call diag_find_rhomax(rhomax1,rhomax2,rholoc1,rholoc2,                 &
                      irholoc1,irholoc2)


!------------------------------------------------------------------------
! make the mask array that tells us whether a given cell
! is part of star1 or star2 by comparing the total energy
! of the fluid in that cell with the average energy at
! the bounary of either star and some simple geometric
! constraints.
call diag_make_mask(factor,rho_boundary,q,rholoc1,syscom,etot1,         &
                    etot2,ntot1,ntot2,mask,frnum)

!------------------------------------------------------------------------
! this section of dodiag finds the location of the Lagrange
! points along the line of centers and the outer boundary of
! each star's Roche lobe.

call diag_find_lpoints(factor,rholoc1,irholoc1,rholoc2,irholoc2,       &
                       lpoints,roche_vol1,roche_vol2)

!------------------------------------------------------------------------
! calculate the total mass of each star and total mass
! of the envelope
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( mask(J,K,L) < 0 ) then
            sum2 = sum2 + rhf(J)*rho_diag(J,K,L)
         else if( mask(J,K,L) > 0 ) then
            sum1 = sum1 + rhf(J)*rho_diag(J,K,L)
         else
            sum3 = sum3 + rhf(J)*rho_diag(J,K,L)
         endif
      enddo
   enddo
enddo

diag_sum(1) = factor*sum1
diag_sum(2) = factor*sum2
diag_sum(3) = factor*sum3

!------------------------------------------------------------------------
! calculate the total mass in the boundary

if( iam_on_bottom .and. (isym == 1) ) then
   do L = philwb, phiupb
      do K = zlwb-1,zlwb
         do J = rlwb, rupb
            mass_boundary = mass_boundary + rhf(J)*rho(J,K,L)
          enddo
      enddo
   enddo
endif

if( iam_on_top ) then
   do L = philwb, phiupb
      do K = zupb, zupb+1
         do J = rlwb, rupb
            mass_boundary = mass_boundary + rhf(J)*rho(J,K,L)
         enddo
      enddo
   enddo
endif

if( iam_on_edge ) then
   if( iam_on_top ) then
      do L = philwb, phiupb
         do K = zlwb, zupb-1
            do J = rupb, rupb+1
               mass_boundary = mass_boundary + rhf(J)*rho(J,K,L)
            enddo
         enddo
      enddo
   else if( iam_on_bottom .and. (isym ==1) ) then
      do L = philwb, phiupb
         do K = zlwb+1, zupb
            do J = rupb, rupb+1
               mass_boundary = mass_boundary + rhf(J)*rho(J,K,L)
            enddo
         enddo
      enddo	
   else 
      do L = philwb, phiupb
         do K = zlwb, zupb
            do J = rupb, rupb+1
               mass_boundary = mass_boundary + rhf(J)*rho(J,K,L)
            enddo	
         enddo
      enddo 
   endif	
endif 

diag_sum(4) = factor * mass_boundary

!------------------------------------------------------------------------
! calculate volumes occupied by material above a given density 
! threshhold in each Roche lobe
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( mask(J,K,L) < 0 ) then
            if( rho_diag(J,K,L) > rho1 ) then
               vol2(1) = vol2(1) + rhf(J)
               vol2(2) = vol2(2) + rhf(J)
               vol2(3) = vol2(3) + rhf(J)
               vol2(4) = vol2(4) + rhf(J)
               vol2(5) = vol2(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho2 ) then
               vol2(2) = vol2(2) + rhf(J)
               vol2(3) = vol2(3) + rhf(J)
               vol2(4) = vol2(4) + rhf(J)
               vol2(5) = vol2(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho3 ) then
               vol2(3) = vol2(3) + rhf(J)
               vol2(4) = vol2(4) + rhf(J)
               vol2(5) = vol2(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho4 ) then
               vol2(4) = vol2(4) + rhf(J)
               vol2(5) = vol2(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho5 ) then
               vol2(5) = vol2(5) + rhf(J)
            endif
         else if( mask(J,K,L) > 0 ) then
            if( rho_diag(J,K,L) > rho1 ) then
               vol1(1) = vol1(1) + rhf(J)
               vol1(2) = vol1(2) + rhf(J)
               vol1(3) = vol1(3) + rhf(J)
               vol1(4) = vol1(4) + rhf(J)
               vol1(5) = vol1(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho2 ) then
               vol1(2) = vol1(2) + rhf(J)
               vol1(3) = vol1(3) + rhf(J)
               vol1(4) = vol1(4) + rhf(J)
               vol1(5) = vol1(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho3 ) then
               vol1(3) = vol1(3) + rhf(J)
               vol1(4) = vol1(4) + rhf(J)
               vol1(5) = vol1(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho4 ) then
               vol1(4) = vol1(4) + rhf(J)
               vol1(5) = vol1(5) + rhf(J)
               cycle
            else if( rho_diag(J,K,L) > rho5 ) then
               vol1(5) = vol1(5) + rhf(J)
            endif	
         endif	
      enddo
   enddo
enddo

diag_sum(5:9) = factor * vol1
diag_sum(10:14) = factor * vol2

!------------------------------------------------------------------------
! calculate the components of the center of mass of each star
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( mask(J,K,L) < 0 ) then
            star2com(1) = star2com(1) +                        &
                          rhf(J)*x(J,L)*rho_diag(J,K,L)
            star2com(2) = star2com(2) +                        &
                          rhf(J)*y(J,L)*rho_diag(J,K,L)
            star2com(3) = star2com(3) +                        &
                          rhf(J)*zhf(K)*rho_diag(J,K,L)
         else if( mask(J,K,L) > 0 ) then
            star1com(1) = star1com(1) +                        &
                          rhf(J)*x(J,L)*rho_diag(J,K,L)
            star1com(2) = star1com(2) +                        &
                          rhf(J)*y(J,L)*rho_diag(J,K,L)
            star1com(3) = star1com(3) +                        &
                          rhf(J)*zhf(K)*rho_diag(J,K,L) 
         endif
      enddo
   enddo
enddo

diag_sum(15:17) = factor * star1com
      diag_sum(18:20) = factor * star2com

!------------------------------------------------------------------------
! calculate the three principal momments of inertia
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         xx = xx + rhf(J)*x(J,L)*x(J,L)*rho_diag(J,K,L)
         yy = yy + rhf(J)*y(J,L)*y(J,L)*rho_diag(J,K,L)
         zz = zz + rhf(J)*zhf(K)*zhf(K)*rho_diag(J,K,L)
      enddo
   enddo
enddo

diag_sum(21) = factor * xx
diag_sum(22) = factor * yy
diag_sum(23) = factor * zz

!------------------------------------------------------------------------
! calculate the gravitational potential energy
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         egrav = egrav + rhf(J)*pot(J,K,L)*rho_diag(J,K,L)
      enddo
   enddo
enddo

diag_sum(24) = 0.5 * factor * egrav

!------------------------------------------------------------------------
! calculate the total angular momentum relative to the inertial frame
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         temp(J,K,L) = vphi(J,K,L) + rhf(J)*omega_frame
      enddo
   enddo
enddo
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         total_j = total_j + rhf(J)*rhf(J)*rho_diag(J,K,L)*temp(J,K,L)
      enddo
   enddo
enddo

diag_sum(25) = factor * total_j

!------------------------------------------------------------------------
! calculate the total kinetic energy of rotation relative to
! the inertial frame
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         trot = trot + rhf(J)*rho_diag(J,K,L)*temp(J,K,L)*temp(J,K,L)
      enddo
   enddo
enddo
  
diag_sum(26) = 0.5 * factor * trot

!------------------------------------------------------------------------
! kinetic energy due to radial and vertical motion
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         trz = trz + rhf(J)*rho_diag(J,K,L)*(                       &
               vr(J,K,L)*vr(J,K,L) + vz(J,K,L)*vz(J,K,L))
      enddo
   enddo
enddo
 
diag_sum(27) = 0.5 * factor * trz

!------------------------------------------------------------------------
! integrate the virial pressure of the fluid
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         pres = pres + rhf(J) * p(J,K,L)
      enddo
   enddo 
enddo

diag_sum(28) = factor * pres

!------------------------------------------------------------------------
! as a check, calculate the total mass on the grid
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         mtot_check = mtot_check + rhf(J)*rho(J,K,L)
      enddo
   enddo
enddo

diag_sum(29) = factor * mtot_check

!------------------------------------------------------------------------------
!
! calculate the total mass of fluid types 1 and 2 to the left and right
! of the boundary plane that is used to sonctruct the mask array.
! Also calculate the total mass of fluid above a density cutoff of
! rho = 1.0e-7
!
divider(1) = syscom(1)*(1.0-q) + rholoc1(1)*cos(rholoc1(3))*q
divider(2) = syscom(2)*(1.0-q) + rholoc1(1)*sin(rholoc1(3))*q
divider(3) = syscom(3)*(1.0-q) + rholoc1(2)*q 

side = syscom - divider
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         test = side(1)*(x(J,L) - divider(1)) +                      &
                side(2)*(y(J,L) - divider(2)) +                      &
                side(3)*(zhf(K) - divider(3))
         if( test > 0.0 ) then
            ! point lies on *2 side
            m1_left = m1_left + frac2(J,K,L)*rhf(J)*rho_diag(J,K,L)
            m2_left = m2_left + frac1(J,K,L)*rhf(J)*rho_diag(J,K,L)
            if( rho(J,K,L) > 1.0e-7 ) then
               mrho_left = mrho_left + rhf(J)*rho_diag(J,K,L)
            endif
         else
            ! point lies on *1 side
            m1_right = m1_right + frac2(J,K,L)*rhf(J)*rho_diag(J,K,L)
            m2_right = m2_right + frac1(J,K,L)*rhf(J)*rho_diag(J,K,L)
            if( rho(J,K,L) > 1.0e-7 ) then
               mrho_right = mrho_right + rhf(J)*rho_diag(J,K,L)
            endif
         endif
      enddo
   enddo
enddo

diag_sum(30) = factor * m1_left
diag_sum(31) = factor * m2_left
diag_sum(32) = factor * m1_right
diag_sum(33) = factor * m2_right
diag_sum(34) = factor * mrho_left
diag_sum(35) = factor * mrho_right

!------------------------------------------------------------------------------
!
! calculate the total angular momentum of the binary about the 
! instantaneous center of mass
!
! >>NOTE<< this loop expects temp to contain vphi|corotating +
! omega * R
!
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         binary_j = binary_j + rhf(J)*abs(mask(J,K,L))*rho_diag(J,K,L)*   &
                     ((x(J,L) - syscom(1))*                               &
                      (sin_cc(L)*vr(J,K,L) + cos_cc(L)*temp(J,K,L)) -     &
                      (y(J,L) - syscom(2))*                               &
                      (cos_cc(L)*vr(J,K,L) - sin_cc(L)*temp(J,K,L)))
      enddo
   enddo
enddo

diag_sum(36) = factor * binary_j

!------------------------------------------------------------------------
! done with first batch of the work, need the center of
! mass coordinates of both stars before we can continue
call mpi_reduce(diag_sum,diag_summed,36,REAL_SIZE,MPI_SUM,            &
                root,MPI_COMM_WORLD,ierror)

if( iam_root ) then
   mass_star1 = diag_summed(1)
   mass_star2 = diag_summed(2)
   star1com = diag_summed(15:17)/mass_star1
   star2com = diag_summed(18:20)/mass_star2 
endif

call mpi_bcast(star1com,3,REAL_SIZE,root,MPI_COMM_WORLD,ierror)

call mpi_bcast(star2com,3,REAL_SIZE,root,MPI_COMM_WORLD,ierror)

!------------------------------------------------------------------------------
! calculate the center of mass coordinates in the inertial frame for each
! star

star1comin(1) = star1com(1)*costemp - star1com(2)*sintemp
star1comin(2) = star1com(1)*sintemp + star1com(2)*costemp

star2comin(1) = star2com(1)*costemp - star2com(2)*sintemp
star2comin(2) = star2com(1)*sintemp + star2com(2)*costemp

star1comin(3) = star1com(3)
star2comin(3) = star2com(3)

!------------------------------------------------------------------------------
! have the center of mass coordinates for each star, ready to
! calculate the spin angular momentum
call diag_spins(factor,star1com,star2com,spin_star1,spin_star2,mask)
      
!------------------------------------------------------------------------------
! calclate momments of the density field
call diag_momments(factor,star1comin,star2comin,mom1,mom2,mask)

!------------------------------------------------------------------------------
! calculate the second time derritive of the momments of
! inertia that are required to find the strain term in
! the quadropole formula for gravitational radiation
call diag_gwave(factor,ilm)

!------------------------------------------------------------------------------
if( iam_root ) then
   ! write everything out
   call diag_write(mtot, syscom, syscomin, star1comin, star2comin,      &
                   rhomax1, rhomax2, rholoc1, rholoc2, irholoc1,        &
                   irholoc2, etot1, etot2, ntot1, ntot2, lpoints,       &
                   roche_vol1, roche_vol2, diag_summed, spin_star1,     &
                   spin_star2, mom1, mom2, ilm, star1com, star2com)
endif 
          
return
end subroutine dodiag
