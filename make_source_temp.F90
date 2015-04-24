!************************************************************************
!*
!*  MAKE_SOURCE_TEMP
!*
!************************************************************************
subroutine make_source_temp 
implicit none
include 'runhydro.h'
!************************************************************************
!*
!   make source temp sources s and a and puts the results in
!   temps and tempa without disturbing the real s and a.
!   These temp values are then used by source in the right
!   hand side of the lagrangian source equations for momentum
!   density
!
!   use the following equations...
!
!   s  =  s  -  dt ( <rho> dPhieff/dR  - <a>^2/(<rho>R^3) - 2 Omegaframe<a>/R )
!                          ^                   ^                   ^
!                          |                   |                   |
!                   gravity, pressure and   curvature of        Coriolis 
!                      centrifugal       cylindrical coords
!
!   a  =  a  -  dt ( rho dPhieff/dphi  +  2 Omegaframe s R)
!                          ^                      ^ 
!                          |                      |
!                   gravity and pressure       Coriolis
!
!  where,
!
!  Phieff = (pin + 1) p / rho + Phi  - 0.5 Omegaframe^2 R^2
!
!         =           H     +   Phi  - 0.5 Omegaframe^2 R^2
!
!         =  the effective hydrodynamic potential 
!
!         =  the Bernoulli function
!
!    the time centering is described below:
!
!        temps  -  s(n + 1/2S + A)
!        -------------------------  =  f(rho(n+1), phieff(n+1), a(n + 1/2S + A))
!                 0.5 dt(n)
!
!        tempa  -  a(n + 1/2S + A)
!        -------------------------  =  g(rho(n+1), phieff(n+1), s(n + 1/2S + A))
!                 0.5 dt(n)
!
!  04/05/2003 : Changed H (= p/rho) to p
!
!  05/24/2003 : Modified for com motion (only com acceleration in the radial 
!               direction).
!
!  06/09/2003 : Modified the radial and azimuthal equations to include
!               the radial and azimuthal component of the acceleration
!               of the com. Refer to subroutine com_vel_accn.F for 
!               details on how the radial and azimuthal components of
!               aceleration of the com are calculated.
!
!*
!**********************************************************************
!*
!*   Global Variables

real,dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd,numz_dd,numphi) :: temps, tempa
common /source_temp/ temps, tempa

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

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

!*
!************************************************************************
!*
!*   Local Variables

real, dimension(numr_dd,numz_dd,numphi) :: temp, rhjkl, phieff

integer :: J, K, L 

!*
!************************************************************************
!  initialize the local variables

! ! 04/05/2003 : Modified (H = p/rho) to p
! from the effective potential aka the Bernoulli function
do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         phieff(J,K,L) = pot(J,K,L) - 0.5*omega_frame*            &
                         omega_frame*rhf(J)*rhf(J)

!          phieff(J,K,L) = (pin+1.0)*p(J,K,L)/rho(J,K,L) +         &
!                          pot(J,K,L) - 0.5*omega_frame*           &
!                          omega_frame*rhf(J)*rhf(J)
      enddo
   enddo
enddo

!  get rho centered at radial faces
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl(J,K,L) = 0.5 * rinv(J) *                           &
                        ( rho(J,K,L)*(r(J) + 0.25*dr) +           &
                          rho(J-1,K,L)*(r(J) - 0.25*dr) )
      enddo
   enddo
enddo

! and then get A at radial faces, put it in temp
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         temp(J,K,L) = 0.25 * rinv(J) *                &
                      ( (a(J,K,L) + a(J,K,L+1)) *      &
                       (r(J) + 0.25*dr) +              &
                        (a(J-1,K,L) + a(J-1,K,L+1)) *  &
                       (r(J) - 0.25*dr) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      temp(J,K,phiupb) = 0.25 * rinv(J) *                       &
                        ( (a(J,K,phiupb) + a(J,K,philwb)) *     &
                         (r(J) + 0.25*dr) +                     &
                          (a(J-1,K,phiupb) + a(J-1,K,philwb)) * &
                        (r(J) - 0.25*dr) )
   enddo
enddo

! 06/17/2003 :
! 06/09/2003 :
! 05/24/2003 : Modified for com motion (only com acceleration in the radial 
!              direction).
! 04/05/2003 : Changed H (= p/rho) to p

!  source s the radial momentum density to get temps
if (tstep .le. 1000) then
    do L = philwb, phiupb
       do K = zlwb, zupb
          do J = rlwb, rupb
             if( rhjkl(J,K,L) > 0.0 ) then
                temps(J,K,L) = s(J,K,L) - dt*(                       &
                    drinv*( (p(J,K,L)-p(J-1,K,L)) +                  &
                    rhjkl(J,K,L)*(phieff(J,K,L) - phieff(J-1,K,L)) ) &
                    - temp(J,K,L)*temp(J,K,L)*rinv(J)*               &
                      rinv(J)*rinv(J)/rhjkl(J,K,L)                   &
                    - 2.0*omega_frame*temp(J,K,L)*rinv(J) )
             else
                temps(J,K,L) = 0.0
             endif
          enddo
       enddo
    enddo
else
    do L = philwb, phiupb
       do K = zlwb, zupb
          do J = rlwb, rupb
             if( rhjkl(J,K,L) > 0.0 ) then
                temps(J,K,L) = s(J,K,L) - dt*(                       &
                    drinv*( (p(J,K,L)-p(J-1,K,L)) +                  &
                    rhjkl(J,K,L)*(phieff(J,K,L) - phieff(J-1,K,L)) ) &
                    - temp(J,K,L)*temp(J,K,L)*rinv(J)*               &
                      rinv(J)*rinv(J)/rhjkl(J,K,L)                   &
                    - 2.0*omega_frame*temp(J,K,L)*rinv(J)            &
                    + rhjkl(J,K,L)* (cylin_a_com(1,L)                &
                      - omega_frame*omega_frame*R_com                &
                      - 2*omega_frame*cylin_v_com(2,L)) )
             else
                temps(J,K,L) = 0.0
             endif
          enddo
      enddo
   enddo
endif

! get rho centered at the azimuthal faces
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         rhjkl(J,K,L) = 0.5 * ( rho(J,K,L) + rho(J,K,L-1) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      rhjkl(J,K,philwb) = 0.5 * ( rho(J,K,philwb) + rho(J,K,phiupb) )
   enddo
enddo

! and then construct S centered at the azimuthal; faces
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         temp(J,K,L) = 0.25 * rhfinv(J) *                         &
                       ( (s(J+1,K,L) + s(J+1,K,L-1))*             &
                       (rhf(J) + 0.25*dr) +                       &
                       (s(J,K,L) + s(J,K,L-1))*                   &
                       (rhf(J) - 0.25*dr) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      temp(J,K,philwb) = 0.25 * rhfinv(J) *                       &
                         ( (s(J+1,K,philwb) + s(J+1,K,phiupb)) *  &
                         (rhf(J) + 0.25*dr) +                     &
                         (s(J,K,philwb) + s(J,K,phiupb)) *        &
                         (rhf(J) - 0.25*dr) )
   enddo
enddo

! 06/09/2003 : Modified to include azimuthal component of the com 
!              acceleration (cylin_a_com(2))
! 04/05/2003 : Changed H (= p/rho) to p
! apply source terms to A
! apply source terms to a to get tempa

if (tstep .le. 1000) then
    do L = philwb+1, phiupb
       do K = zlwb, zupb
          do J = rlwb, rupb
             tempa(J,K,L) = a(J,K,L) - dt*(                            &
                    dphiinv*( (p(J,K,L)-p(J,K,L-1)) +                  &
                    rhjkl(J,K,L)*(phieff(J,K,L) - phieff(J,K,L-1)) )   &
                    + 2.0*omega_frame*rhf(J)*temp(J,K,L) )
          enddo
       enddo
    enddo
    do K = zlwb, zupb
       do J = rlwb, rupb
          tempa(J,K,philwb) = a(J,K,philwb) - dt*(                     &
                    dphiinv*( (p(J,K,philwb)-p(J,K,phiupb)) +          & 
                    rhjkl(J,K,philwb) *                                &
                    (phieff(J,K,philwb) - phieff(J,K,phiupb)) )        &
                    + 2.0*omega_frame*rhf(J)*temp(J,K,philwb) )
       enddo
    enddo
else
    do L = philwb+1, phiupb
       do K = zlwb, zupb
          do J = rlwb, rupb
             tempa(J,K,L) = a(J,K,L) - dt*(                           &
                    dphiinv*( (p(J,K,L)-p(J,K,L-1)) +                 &
                    rhjkl(J,K,L)*(phieff(J,K,L) - phieff(J,K,L-1)) )  &
                    + 2.0*omega_frame*rhf(J)*temp(J,K,L)              &
                    + rhjkl(J,K,L)*rhf(J)* (cylin_a_com(2,L)          &
                       + 2*omega_frame*cylin_v_com(1,L)) )
          enddo
       enddo
    enddo
    do K = zlwb, zupb
       do J = rlwb, rupb
          tempa(J,K,philwb) = a(J,K,philwb) - dt*(                    &
                    dphiinv*( (p(J,K,philwb)-p(J,K,phiupb)) +         &
                    rhjkl(J,K,philwb) *                               &
                    (phieff(J,K,philwb) - phieff(J,K,phiupb)) )       &
                    + 2.0*omega_frame*rhf(J)*temp(J,K,philwb)         &
                    + rhjkl(J,K,philwb)*rhf(J)*                       &
                      (cylin_a_com(2,philwb)                          &
                      + 2*omega_frame*cylin_v_com(1,philwb)) )
       enddo
    enddo
endif
                  
call make_source_temp_bc

!   need to load guard cell layers with temps and tempa
call comm(temps)
call comm(tempa)

return
end subroutine make_source_temp
