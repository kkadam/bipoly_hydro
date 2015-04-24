!************************************************************************
!*
!*  DIAG_WRITE.1
!*
!************************************************************************
subroutine diag_write(mtot, syscom, syscomin, star1comin,               &
                      star2comin, rhomax1, rhomax2, rholoc1,            &
                      rholoc2, irholoc1, irholoc2, etot1,               &
                      etot2, ntot1, ntot2, lpoints, roche_vol1,         &
                      roche_vol2, diag_sum, spin_star1, spin_star2,     &
                      mom1, mom2, ilm, star1com, star2com)
implicit none
include 'runhydro.h'
!************************************************************************
!*
!   diag_write formats and outputs all the data we have come up
!   with in the dodiag part of the hydrocode package
!
!   the answer key for these variables is as follows:
!
!   mtot: total mass on the grid excluding boundary zones
!
!   syscom: (x,y,z) coordinates of the system center of mass in
!           the corotating frame
!
!   syscomin: (x,y,z) coordinates of the system center of mass in
!             the inertial frame
!
!   star[1-2]comin: (x,y,z) coordinates of the center of mass of the
!                   two stars in the inertial frame
!
!   rhomax[1-2]: maximum density value for each star, rhomax1 is
!                the largest density value on the grid
!
!   rholoc[1-2]: (R,z,Phi) coordinates of the maximum density point
!                for each star
!
!   irholoc[1-2]: (j,k,l) grid coordinates for the density max of
!                 each star
!
!   etot[1-2]: total energy of all cells identified as being on the
!              boundary of the two stars
!
!   ntot[1-2]: total number of cells in the boundary layes for each 
!              star
!
!   lpoints: first index runs through the potential and (x,y,z) coordinates
!            of the point in question.  For the second index:
!            1 => L1 point
!            2 => outer Lagrange point for *1
!            3 => outer boundary of *1's Roche lobe
!            4 => outer Lagrange point for *2
!            5 => outer boundary of *2's Roche lobe
!
!   roche_vol[1-2]: total volume occupied by a stars Roche lobe
!
!   diag_sum(1): total mass of *1
!
!   diag_sum(2): total mass of *2
! 
!   diag_sum(3): total mass of material on the active grid that is in
!                neither star
!
!   diag_sum(4): total mass in the boundary cells
!
!   diag_sum(5-9): volume occupied by material in *1 above rho[1-5]
!                  density level
!
!   diag_sum(10-14): volume occupied by material in *2 above rho[1-5]
!                     density level
!
!   diag_sum(15-17): (x,y,z) coordinates of center of mass for *1
!
!   diag_sum(18-20): (x,y,z) coordinates of center of mass for *2
!
!   diag_sum(21-23): three quadratic momments of the system's density
!                    distribution, relative to corotating frame
!
!   diag_sum(24): gravitational potential energy of the system
!
!   diag_sum(25): total angular momenta of the system, relative to
!                 the inertial frame
!
!   diag_sum(26): kinetic energy of rotation for the system, relative to
!                 the inertial frame
!
!   diag_sum(27): total kinetic energy due to vertical and radial
!                 motion
!
!   diag_sum(28): integral of the pressure over the system
!
!   diag_sum(29): simple total mass, including boundary zones 
!
!   diag_sum(30): frac1 * rho to left of divider
!   diag_sum(31): frac2 * rho to left of divider
!   diag_sum(32): frac1 * rho to right of divider
!   diag_sum(33): frac2 * rho to right of divider
!
!   diag_sum(34-35): mass of material to left and right of divider that
!                    is more dense than 1.0e-7
!
!   diag_sum(36): total angular momentum of material in either star (by energy
!                 determination) in the inertial frame (about a moving z axis).   
!
!   spin[1-2]: z compnent of spin angular momenta for each star, relative
!              to the inertial frame
!
!   mom[1-2]: quadratic, cubic and quartic polynomials in the Cartessian
!             coordinates over the density distribution for each star.
!             origin of the coordinate system is the instantaneous center
!             of mass of the system and uses inertial coordinates
!
!   ilm:  6 components of the second time derivative of the momment of
!         inertia tensor for the system.  Integration is done over inertial
!         space coordinates
!
!
!  Modified 5/20/2000 to write out a new line of output that has 
!  diag_sum(30:36)
!
!*
!************************************************************************
!*
!*   Global variables

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

!*
!************************************************************************
!*
!*   Local Variables


!*
!************************************************************************
!*
!*    Subroutine Arguments

real, dimension(31) :: mom1, mom2

real, dimension(36) :: diag_sum

real, dimension(4,5) :: lpoints

real, dimension(6) :: ilm

real, dimension(3) :: syscom, syscomin, star1comin, star2comin,       &
                      rholoc1, rholoc2, star1com, star2com

real :: mtot, rhomax1, rhomax2, etot1, etot2, ntot1, ntot2,           &
        roche_vol1, roche_vol2, spin_star1, spin_star2

integer, dimension(3) :: irholoc1, irholoc2

!*
!************************************************************************ 
!format statements for the output
  100  format(i7,' timestep ',6e16.8)

  110  format(' angular momentum ',6e16.8)

  120  format(' masses ',7e16.8)

  130  format(' density ',8e16.8,6i4)

  140  format(' position ',18e16.8)

  150  format(' volume1 ',5e16.8)

  160  format(' volume2 ',5e16.8)

  170  format(' gravrad ',6e16.8)

  180  format(' Roche ',22e16.8)

  190  format(' momments1 ',31e16.8)

  200  format(' momments2 ',31e16.8)

  210  format(' boundary ',6e16.8)

  220  format(' newmass ',7e16.8)

  230  format(' com_coords ', 6e16.8)

  240  format(' com_vel_xyz ', 3e16.8)

  250  format(' com_cylin_vel ', 768e16.8)

  260  format(' com_acc_xyz ', 3e16.8)

  270  format(' com_cylin_acc ', 768e16.8)

!  280  format(' source_terms ', 3e16.8)

! write everything out 
write(6,100) tstep, time/cirp, dt, diag_sum(24), diag_sum(26),       &
             diag_sum(27), diag_sum(28)

write(6,110) diag_sum(25), spin_star1, spin_star2, diag_sum(21),     &
             diag_sum(22), diag_sum(23)

write(6,120) diag_sum(1), diag_sum(2), diag_sum(3), diag_sum(4),     &
             diag_sum(1)+diag_sum(2)+diag_sum(3)+diag_sum(4),        &
             mtot, diag_sum(29)

write(6,130) rhomax1, rhomax2, rholoc1, rholoc2, irholoc1, irholoc2

write(6,140) syscom, star1com, star2com, syscomin, star1comin,       &
             star2comin

write(6,150) diag_sum(5:9)

write(6,160) diag_sum(10:14)

write(6,170) ilm

write(6,180) roche_vol1, roche_vol2, lpoints

write(6,190) mom1

write(6,200) mom2

write(6,210) etot1, ntot1, etot2, ntot2, etot1/ntot1, etot2/ntot2

write(6,220) diag_sum(30:36)

write(6,230) delt(1), phi_com, com(1,1), com(1,2), com(1,3), R_com

write(6,240) v_com

!write(6,250) cylin_v_com

write(6,260) a_com

!write(6,270) cylin_a_com

!write(6,280) -omega_frame**2*R_com, -2*omega_frame*cylin_v_com(2), 2*omega_frame*cylin_v_com(1)
          
return
end subroutine diag_write
