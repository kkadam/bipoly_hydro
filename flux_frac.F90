!**************************************************************************
!*
!*  FLUX_FRAC
!*
!**************************************************************************
subroutine flux_frac(rho_tmp, slope_r, slope_z, slope_phi, varea_r, varea_z, &
                varea_phi, deltat)
implicit none
include 'runhydro.h'
!**************************************************************************
!
!  Modified 5/17/2000 to keep track of the mass fraction during the
!   advection of mass density.  A very limited way to implement this
!   but it fits my needs, I only need to distinguish between two fluid
!   components and I simply asume that they are evenly mixed in every
!   cell at all times.
!
!**************************************************************************
!*
!*  Subroutine Arguments
      
real rho_tmp(numr_dd,numz_dd,numphi)
real slope_r(numr_dd,numz_dd,numphi)
real slope_z(numr_dd,numz_dd,numphi)
real slope_phi(numr_dd,numz_dd,numphi)
real varea_r(numr_dd,numz_dd,numphi)
real varea_z(numr_dd,numz_dd,numphi)
real varea_phi(numr_dd,numz_dd,numphi)
real deltat

!*
!**************************************************************************
!*
!*

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

real, dimension(numr_dd) :: vol_cc, volinv_cc, ar_cc, az_cc
real :: aphi_cc
common /geometry_cc/ vol_cc, volinv_cc, ar_cc, az_cc, aphi_cc

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv  

!*
!**************************************************************************
!*
!*  Local Variables

real, dimension(numr_dd,numz_dd,numphi) :: dm
real, dimension(numr_dd,numz_dd,numphi) :: dm1, dm2
real :: qdot, mass_res_limit
integer ::  J, K, L

!*
!**************************************************************************
!  initialize the local variables
!
!  the mass resolution limit is set to be the floor value for the density
!  that will be followed times the smallest cell volume times 0.01.
!  the factor of 0.01 being pretty much an arbitrary choice but it
!  does have the implcation that the code will lose track of mass type in a
!  cell if it is emptied further than this ammount during a timestep
!
qdot = 0.0
mass_res_limit = 5.0e-3*densmin*dr*dr*dz*dphi
dm = 0.0
dm1 = 0.0
dm2 = 0.0

!  change in cell's content of mass due to advection in radial direction 
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb-1, rupb
         if( varea_r(J,K,L) > 0.0 ) then
            qdot = (rho_tmp(J,K,L) + 0.5*slope_r(J,K,L))*             &
                   varea_r(J,K,L)
            dm1(J,K,L)   = dm1(J,K,L)   - frac1(J,K,L) * qdot * deltat
            dm1(J+1,K,L) = dm1(J+1,K,L) + frac1(J,K,L) * qdot * deltat
            dm2(J,K,L)   = dm2(J,K,L)   - frac2(J,K,L) * qdot * deltat
            dm2(J+1,K,L) = dm2(J+1,K,L) + frac2(J,K,L) * qdot * deltat
         else
            qdot = (rho_tmp(J+1,K,L) - 0.5*slope_r(J+1,K,L))*         &
                   varea_r(J,K,L)
            dm1(J,K,L)   = dm1(J,K,L)   - frac1(J+1,K,L) * qdot * deltat
            dm1(J+1,K,L) = dm1(J+1,K,L) + frac1(J+1,K,L) * qdot * deltat
            dm2(J,K,L)   = dm2(J,K,L)   - frac2(J+1,K,L) * qdot * deltat
            dm2(J+1,K,L) = dm2(J+1,K,L) + frac2(J+1,K,L) * qdot * deltat
         endif
         dm(J,K,L)   = dm(J,K,L)   - qdot * deltat
         dm(J+1,K,L) = dm(J+1,K,L) + qdot * deltat
      enddo
   enddo
enddo

!  same for vertical direction
do L = philwb, phiupb
   do K = zlwb-1, zupb
      do J = rlwb, rupb
         if( varea_z(J,K,L) > 0.0 ) then
            qdot = (rho_tmp(J,K,L) + 0.5*slope_z(J,K,L))*             &
                   varea_z(J,K,L)
            dm1(J,K,L)   = dm1(J,K,L)   - frac1(J,K,L) * qdot * deltat
            dm1(J,K+1,L) = dm1(J,K+1,L) + frac1(J,K,L) * qdot * deltat
            dm2(J,K,L)   = dm2(J,K,L)   - frac2(J,K,L) * qdot * deltat
            dm2(J,K+1,L) = dm2(J,K+1,L) + frac2(J,K,L) * qdot * deltat
         else
            qdot = (rho_tmp(J,K+1,L) - 0.5*slope_z(J,K+1,L))*         &
                   varea_z(J,K,L)
            dm1(J,K,L)   = dm1(J,K,L)   - frac1(J,K+1,L) * qdot * deltat
            dm1(J,K+1,L) = dm1(J,K+1,L) + frac1(J,K+1,L) * qdot * deltat
            dm2(J,K,L)   = dm2(J,K,L)   - frac2(J,K+1,L) * qdot * deltat
            dm2(J,K+1,L) = dm2(J,K+1,L) + frac2(J,K+1,L) * qdot * deltat
         endif
         dm(J,K,L)   = dm(J,K,L)   - qdot * deltat
         dm(J,K+1,L) = dm(J,K+1,L) + qdot * deltat
      enddo
   enddo
enddo
     
!  same for azimuthal direction
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( varea_phi(J,K,L) > 0.0 ) then
            qdot = (rho_tmp(J,K,L) + 0.5*slope_phi(J,K,L))*            &
                   varea_phi(J,K,L)
            dm1(J,K,L)   = dm1(J,K,L)   - frac1(J,K,L) * qdot * deltat
            dm1(J,K,L+1) = dm1(J,K,L+1) + frac1(J,K,L) * qdot * deltat
            dm2(J,K,L)   = dm2(J,K,L)   - frac2(J,K,L) * qdot * deltat
            dm2(J,K,L+1) = dm2(J,K,L+1) + frac2(J,K,L) * qdot * deltat
         else
            qdot = (rho_tmp(J,K,L+1) - 0.5*slope_phi(J,K,L+1))*        &
                           varea_phi(J,K,L)
            dm1(J,K,L)   = dm1(J,K,L)   - frac1(J,K,L+1) * qdot * deltat
            dm1(J,K,L+1) = dm1(J,K,L+1) + frac1(J,K,L+1) * qdot * deltat
            dm2(J,K,L)   = dm2(J,K,L)   - frac2(J,K,L+1) * qdot * deltat
            dm2(J,K,L+1) = dm2(J,K,L+1) + frac2(J,K,L+1) * qdot * deltat
         endif
         dm(J,K,L)   = dm(J,K,L)   - qdot * deltat
         dm(J,K,L+1) = dm(J,K,L+1) + qdot * deltat
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      if( varea_phi(J,K,numphi) > 0.0 ) then
         qdot = (rho_tmp(J,K,numphi) + 0.5*slope_phi(J,K,numphi))*     &
                varea_phi(J,K,numphi)
         dm1(J,K,numphi) = dm1(J,K,numphi) - frac1(J,K,numphi) * &
                                                 qdot * deltat
         dm1(J,K,1)      = dm1(J,K,1)      + frac1(J,K,numphi) * &
                                                 qdot * deltat 
         dm2(J,K,numphi) = dm2(J,K,numphi) - frac2(J,K,numphi) * &
                                                 qdot * deltat
         dm2(J,K,1)      = dm2(J,K,1)      + frac2(J,K,numphi) * &
                                                 qdot * deltat
      else
         qdot = (rho_tmp(J,K,1) - 0.5*slope_phi(J,K,1))*               &
                varea_phi(J,K,numphi)
         dm1(J,K,numphi) = dm1(J,K,numphi) - frac1(J,K,1) *      &
                                                 qdot * deltat       
         dm1(J,K,1)      = dm1(J,K,1)      + frac1(J,K,1) *      &
                                                 qdot * deltat
         dm2(J,K,numphi) = dm2(J,K,numphi) - frac2(J,K,1) *      &
                                                 qdot * deltat
         dm2(J,K,1)      = dm2(J,K,1)      + frac2(J,K,1) *      &
                                                 qdot * deltat
      endif
      dm(J,K,numphi) = dm(J,K,numphi) - qdot * deltat
      dm(J,K,1)      = dm(J,K,1)      + qdot * deltat
   enddo
enddo

 
! calculate the new mass fraction
!
!  [m1/mtot] = [( m1 + dm1 ) / m + dm ]
!
!  [m2/mtot] = [( m2 + dm2 ) / m + dm ]
!
do L = philwb, phiupb
   do K = zlwb, zupb
      do  J = rlwb, rupb

          if( (rho_tmp(J,K,L)*vol_cc(J) + dm(J,K,L)) < mass_res_limit ) then

             frac1(J,K,L) = 0.0
           
             frac2(J,K,L) = 0.0

          else
    
             if( (frac1(J,K,L)*rho_tmp(J,K,L)*vol_cc(J) + dm1(J,K,L)) < 0.0 ) then
 
                frac1(J,K,L) = 0.0

             else

                frac1(J,K,L) = ( frac1(J,K,L) * rho_tmp(J,K,L) * vol_cc(J) + dm1(J,K,L) ) / &
                               ( rho_tmp(J,K,L) * vol_cc(J) + dm(J,K,L) )

             endif

             if( (frac2(J,K,L)*rho_tmp(J,K,L)*vol_cc(J) + dm2(J,K,L)) < 0.0 ) then

                frac2(J,K,L) = 0.0

             else

               frac2(J,K,L) = ( frac2(J,K,L) * rho_tmp(J,K,L) * vol_cc(J) + dm2(J,K,L) ) / &
                              ( rho_tmp(J,K,L) * vol_cc(J) + dm(J,K,L) )

            endif
         endif
      enddo
   enddo
enddo


return
end subroutine flux_frac
