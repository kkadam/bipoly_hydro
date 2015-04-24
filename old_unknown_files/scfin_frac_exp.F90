!***********************************************************************
!*
!*  SCFIN_FRAC
!*
!***********************************************************************
subroutine scfin
implicit none
include 'runhydro.h'
include 'mpif.h'
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

open(unit=21,file='density_130_98_256_fc',                           &
     form='unformatted',status='old',access='direct',                  &
     recl=record_length)

open(unit=22,file='a_130_98_256_fc',                                 &
     form='unformatted',status='old',access='direct',                  &
     recl=record_length)

open(unit=23,file='s_130_98_256_fc',                                 &
     form='unformatted',status='old',access='direct',                  &
     recl=record_length)
       
read(21,rec=iam+1) rho

read(22,rec=iam+1) a

read(23,rec=iam+1) s

close(21)
close(22)
close(23)

! set up the mass fraction array
!
! if fluid is of density greater then 10 times the vacuum level 
! and is to the left of x = 0.1 then call it fluid of type 1
!
do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         frac1(J,K,L) = 0.0
         frac2(J,K,L) = 0.0
         if( rho(J,K,L) > 10.0 * densmin ) then
            if( rhf(J) * cos_cc(L) < 0.1 ) then
               frac1(J,K,L) = 1.0
            else
               frac2(J,K,L) = 1.0
            endif
         endif 
      enddo
   enddo
enddo

! calculate the center of mass for each body
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( x >= 0.0 ) then
            sum_it(1) = sum_it(1) + rhf(J)*rho(J,K,L)
            sum_it(2) = sum_it(2) + rhf(J)*rhf(J)*cos_cc(L)rho(J,K,L)
            sum_it(3) = sum_it(3) + rhf(J)*rhf(J)*sin_cc(L)rho(J,K,L)
            sum_it(4) = sum_it(4) + rhf(J)*zhf(K)*rho(J,K,L)
         else
            sum_it(5) = sum_it(5) + rhf(J)*rho(J,K,L)
            sum_it(6) = sum_it(6) + rhf(J)*rhf(J)*cos_cc(L)*rho(J,K,L)
            sum_it(7) = sum_it(7) + rhf(J)*rhf(J)*sin_cc(L)*rho(J,K,L)
            sum_it(8) = sum_it(8) + rhf(J)*zhf(K)*rho(J,K,L)
         endif
      enddo
   enddo
enddo

call mpi_reduce(sum_it,summed_it,8,REAL_SIZE,MPI_SUM,root,   &
                MPI_COMM_WORLD,ierror)

if( iam_root ) then
   summed_it(2:4) = summed_it(2:4) / summed_it(1)
   summed_it(6:8) = summed_it(6:8) / summed_it(5)
   write(*,*) 'center of mass calculation in scfin: ',summed_it
endif

call mpi_bcast(summed_it,8,REAL_SIZE,root,MPI_COMM_WORLD,ierror)

x1 = summed_it(2)
y1 = summed_it(3)
z1 = summed_it(4)

x2 = summed_it(6)
y2 = summed_it(7)
z2 = summed_it(8)


! given the density generate the internal energy per unit mass
!
!  p  =  kappa * rho ** gamma  =  (gamma - 1) * rho * eps
!
!  so that
! 
! eps = kappa * rho**(1/n) / (gamma - 1)
!
temp = gamma - 1.0

do L = philwb, phiupb
   do K = zlwb-1, zupb+1
      do J = rlwb-1, rupb+1
         if( rhf(J) * cos_cc(L) >= 0.0 ) then
           eps(J,K,L) = (kappa1 * rho(J,K,L)**(1.0/pin))/temp
         else
           eps(J,K,L) = (kappa2 * rho(J,K,L)**(1.0/pin))/temp
         endif
      enddo
   enddo
enddo

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
