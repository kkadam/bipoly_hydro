!***********************************************************************
!*
!*  SCFIN_FRAC
!*
!***********************************************************************
subroutine scfin
implicit none
include 'runhydro.h'
!***********************************************************************
!*
!  scfin reads in the three dimensional arrays for the density,
!  anglular momentum density and radial momentum density.  Then
!  we set the internal energy density array from the mass density
!  and the appropriate polytropic constant for each star.
!
!  Assume that an external program generates the initial conditions
!  for all data on all processors.  This means that guard cell
!  values are assumed to be valid in the input file and appropriate
!  boundary conditions are enforced as well.
!
!  >>NOTE<< for isym = 1, the input files have to extend the entire
!           range of z values.  Communication to replicate data
!           from upper half space to lower half space would be
!           too big a pain in the neck
!
!  5/17/2000 modified to set up the fluid's mass fraction array 
!  for an initial model.  The fluid is all star1 if it is left
!  of x = 0.1 initially... hard coded and stupid.
!*
!***********************************************************************
!*
!*  Global Variable

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

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

real :: kappac1, kappac2, rho_c1, rho_c2, np1, np2   !bipoly 
common /bipoly/ kappac1, kappac2, rho_c1, rho_c2, np1, np2

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

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
!***********************************************************************
!*
!*  Local Variable

real, dimension(numr_dd,numz_dd,numphi) :: eps

real :: temp 
 
integer :: record_length, three_numphi_by_four

integer :: J, K, L

!*
!***********************************************************************
!  initialize the local variables
eps = 0.0
temp = 0.0
three_numphi_by_four = 3 * numphi_by_four
inquire(iolength=record_length) rho

open(unit=21, file='input/density', &
     form='unformatted', status='old', access='direct', &
     recl=record_length)

open(unit=22, file='input/ang_mom', &
     form='unformatted', status='old', access='direct', &
     recl=record_length)

open(unit=23, file='input/rad_mom', &
     form='unformatted', status='old', access='direct', &
     recl=record_length)
       
read(21,rec=iam+1) rho

read(22,rec=iam+1) a

read(23,rec=iam+1) s

close(21)
close(22)
close(23)

! 10/08/2003 : Modified x to -0.2
! set up the mass fraction array
!
! if fluid is of density greater then 10 times the vacuum level 
! and is to the left of x = -0.2 then call it fluid of type 1
!
do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         frac1(J,K,L) = 0.0
         frac2(J,K,L) = 0.0
         if( rho(J,K,L) > 10.0 * densmin ) then
            if( rhf(J) * cos_cc(L) < 0.15 ) then
               frac2(J,K,L) = 1.0
            else
               frac1(J,K,L) = 1.0
            endif
         endif 
      enddo
   enddo
enddo


! given the density generate the internal energy per unit mass
!
!  p  =  kappa * rho ** gamma  =  (gamma - 1) * rho * eps
!
!  so that
! 
! eps = kappa * rho**(1/n) / (gamma - 1)
!
temp = gamma - 1.0

!do L = philwb, phiupb    !old
!   do K = zlwb-1, zupb+1
!      do J = rlwb-1, rupb+1
!         if( rhf(J) * cos_cc(L) >= 0.15 ) then
!           eps(J,K,L) = (kappa1 * rho(J,K,L)**(1.0/pin))/temp
!         else
!           eps(J,K,L) = (kappa2 * rho(J,K,L)**(1.0/pin))/temp
!         endif
!      enddo
!   enddo
!enddo
!temp = 1.0/pin


   do L = philwb, phiupb            !current
      do K = zlwb-1, zupb+1
         do J = rlwb-1, rupb+1
            if( rhf(J) * cos_cc(L) >= 0.15 ) then
               if (rho(j,k,l).gt.rho_c1) then
                  eps(J,K,L) = (kappac1 * rho(J,K,L)**(1.0/np1))/temp
               else
                  eps(J,K,L) = (kappa1 * rho(J,K,L)**(1.0/np2))/temp
               endif
            else
               if (rho(j,k,l).gt.rho_c2) then
                  eps(J,K,L) = (kappac2 * rho(J,K,L)**(1.0/np1))/temp
               else
                  eps(J,K,L) = (kappa2 * rho(J,K,L)**(1.0/np2))/temp
               endif
            endif
         enddo
      enddo
   enddo


! do L = philwb, phiupb            !bipoly
!      do K = zlwb-1, zupb+1
!         do J = rlwb-1, rupb+1
!            if( rhf(J) * cos_cc(L) >= 0.15 ) then
!               if (rho(j,k,l).gt.rho_c1) then
!                  eps(J,K,L) = (kappac1 * rho(J,K,L)**(1.0/np1))/temp1
!               else
!                  eps(J,K,L) = (kappa1 * rho(J,K,L)**(1.0/np2))/temp2
!               endif
!            else
!               if (rho(j,k,l).gt.rho_c2) then
!                  eps(J,K,L) = (kappac2 * rho(J,K,L)**(1.0/np1))/temp1
!               else
!                  eps(J,K,L) = (kappa2 * rho(J,K,L)**(1.0/np2))/temp2
!               endif
!            endif
!         enddo
!      enddo
!   enddo


if( iam_root ) then
open(unit=76, file='output/test.out', &
     form='formatted',status='unknown', position='append')
write(76,*) "temp", temp
write(76,*) "rho_c1", rho_c1
write(76,*) "rho_c2",rho_c2
write(76,*) "kappa1",kappa1
write(76,*) "kappa2",kappa2
write(76,*) "kappac1",kappac1
write(76,*) "kappac2",kappac2
write(76,*) "np1",np1
write(76,*) "np2",np2
!write(76,*) "temp", temp, "temp1",temp1, "temp2",temp2
!write(76,*) "eps(J,K,L) = (kappac1 * 0.9**(1.0/np1))/temp1 = ",(kappac1 * 0.9**(1.0/np1))/temp1
!write(76,*) "eps(J,K,L) = (1.6e-2 * 0.9**(1.0/pin))/temp = ", (1.6e-2 * 0.9**(1.0/pin))/temp


close(76)
endif




! now convert eps to be the entropy tracer, tau
! 
!  tau  =  (eps * rho) ** (1/gamma)
!
! and tau is related to the entropy of the fluid as follows
!
!  s  =  c_p * ln( tau / rho ) + constant
!
!  for adiabatic flow tau has no source term in the eulerian
!  fluid equations (just like rho).  this is not true if
!  the effects of viscosity (artificial or real) are included. 
tau = (eps*rho)**gammainv

return
end subroutine scfin
