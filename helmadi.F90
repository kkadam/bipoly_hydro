!********************************************************************
!*
!*  HELMADI
!*
!********************************************************************
subroutine helmadi(nsteps)
implicit none
include 'runhydro.h'
!********************************************************************
!*
!  helmadi solves Poisson's equation for the self-gravity potential
!  First, Fourier transform the density and initial potential to
!  decouple the azimuthal direction.  Then use ADI method to 
!  iteratively solve for the potential and inverse Fourier transform
!  the potential to get back to physical space.
! 
!  Initially implemented by Howard Cohl in hpf.  See Howie's
!  theis for original hpf source code.  For a discussion
!  of the ADI method to solve Poisson's equation in cylindrical
!  coordinates see Black and Bodenheimer, ApJ vol 199, p 619, 1975
!*
!********************************************************************
!*
!*  Subroutine Arguemtns

integer :: nsteps

!*
!********************************************************************
!*
!*  Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: potp, rhop
common /potarrays/ potp, rhop

real, dimension(numr) :: ar, cr, alphar
real, dimension(numr,numphi_dd) :: brb
common /ADI_R_sweep/ ar, cr, alphar, brb

real :: az, cz
real, dimension(numr_dd_z) :: alphaz, betaz
real, dimension(numz) :: bzb
real, dimension(numr_dd_z,numphi_dd) :: elambdazb
common /ADI_Z_sweep/ az, cz, alphaz, betaz, bzb, elambdazb

real :: gamma, piinv, four_pi
common /pot_constants/ gamma, piinv, four_pi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g, r_g, rhfinv_g, rinv_g, zhf_g

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
!********************************************************************
!*
!*  Local Variables
  
real, dimension(numr_dd,numz_dd,numphi) :: ffrho, ffpot

real, dimension(numr,numz_dd,numphi_dd) :: knownr, rhor, potr

real, dimension(numr_dd_z,numz,numphi_dd) :: knownz, rhoz, potz

real, dimension(numz) :: bz

real, dimension(numr,numphi_dd) :: br

real, dimension(numr_dd_z,numphi_dd) :: elambdaz

integer :: I, J, K, L, index

real :: alph, dtt, dtt_mtwogamma

real, dimension(nsteps) :: dt

!*
!********************************************************************
! initialize the local variables
ffrho =  0.0
ffpot = 0.0
knownr = 0.0
rhor = 0.0
potr = 0.0
knownz = 0.0
rhoz = 0.0
potz = 0.0
bz = 0.0
br = 0.0
elambdaz= 0.0
index = 0
alph = 0.0
dtt = 0.0
dtt_mtwogamma = 0.0 
dt = 0.0

!  setup the pseudo timesteps for the ADI iterations
!  the values for the timestep have been chosen to
!  optimize convergence, see Black and Bodenheimer
dt(nsteps) = 4.0*(rinv_g(numr)*rinv_g(numr))
alph = (r_g(numr)*rinv_g(3))**(2.0/(nsteps-1.0))
do I = 2, nsteps
   index = nsteps + 1 - I 
   dt(index) = alph*dt(index+1)
enddo 

!  fourier transform the density and re-order the output
call realft(rhop,numr_dd,numz_dd,numphi,+1)
do K = 1, numz_dd
   do J = 1, numr_dd
      ffrho(J,K,1) = rhop(J,K,1)
      ffrho(J,K,numphi_by_two+1) = rhop(J,K,2)
   enddo
enddo
index = 3
do L = 2, numphi_by_two
   do K = 1, numz_dd
      do J = 1, numr_dd
         ffrho(J,K,L) = rhop(J,K,index)
      enddo
   enddo
   index = index + 2
enddo
index = 4
do L = numphi_by_two+2, numphi
   do K = 1, numz_dd
      do J = 1, numr_dd
         ffrho(J,K,L) = - rhop(J,K,index)
      enddo
   enddo
   index = index + 2
enddo

!  fourier transform the potential and re-order
call realft(potp,numr_dd,numz_dd,numphi,+1)
do K = 1, numz_dd
   do J = 1, numr_dd
      ffpot(J,K,1) = potp(J,K,1)
      ffpot(J,K,numphi_by_two+1) = potp(J,K,2)
   enddo
enddo
index = 3
do L = 2, numphi_by_two
   do K = 1, numz_dd
      do J = 1, numr_dd
         ffpot(J,K,L) = potp(J,K,index)
      enddo
   enddo
   index = index + 2
enddo
index = 4
do L = numphi_by_two+2, numphi
   do K = 1, numz_dd
      do J = 1, numr_dd
         ffpot(J,K,L) = - potp(J,K,index)
      enddo
   enddo
   index = index + 2
enddo

!  make sure all pes are starting with correct guard cell values
call comm(ffpot)
call comm(ffrho)
 
!  swap arrays around to get data aligned for ADI iteration
call swap_phi_r(ffpot,potr)
call swap_phi_r(ffrho,rhor)
call swap_r_z(rhor,rhoz)

!  ADI iteration cycle
do I = 1, nsteps

   dtt = dt(I)

   dtt_mtwogamma = dtt - 2.0*gamma

   ! Radial ADI sweep

   do L = 1, numphi_dd
      do J = 1, numr
         br(J,L) = brb(J,L) + dtt
      enddo
   enddo

   knownr = 0.0

   if( iam_on_bottom .and. isym /= 1 ) then
      !  the conditional looks a little stange, here is the deal...
      !  in this data decomposition, the vertical index is block
      !  distributed across numz_procs, for the global K index of
      !  two need to make a special case if isym is 2 or 3 which
      !  is that you don't include the potential at K = 1 in
      !  knownr with K = 2 and the weighting of potr with K = 2
      !  is different.  Do this to impose symmetry between solution
      !  at K = 1 and K = 2 when solution at K = 2 is not known
      !  initially
      do L = 1, numphi_dd
         do K = 3, numz_dd-1
            do J = 2, numr-1
               knownr(J,K,L) = -four_pi * rhor(J,K,L) +                &
                                dtt_mtwogamma*potr(J,K,L) +            &
                                gamma*(potr(J,K+1,L) + potr(J,K-1,L))
            enddo
         enddo
      enddo
      ! the special case
      K = 2
      do L = 1, numphi_dd
         do J = 2, numr-1
            knownr(J,K,L) = - four_pi * rhor(J,K,L) +                  &
                              (dtt-gamma)*potr(J,K,L) +                &
                              gamma*potr(J,K+1,L)
         enddo
      enddo
   else
      do L = 1, numphi_dd
         do K = 2, numz_dd-1
            do J = 2, numr-1
               knownr(J,K,L) = - four_pi * rhor(J,K,L) +               &
                                 dtt_mtwogamma*potr(J,K,L) +           &
                                 gamma*(potr(J,K+1,L) + potr(J,K-1,L))
            enddo
         enddo
      enddo
   endif

   !  add in the boundary potential on side to knownr
   do L = 1, numphi_dd
      do K = 2, numz_dd-1
         knownr(numr-1,K,L) = knownr(numr-1,K,L) -        &
                     alphar(numr-1)*potr(numr,K,L)
      enddo
   enddo

   !  solve the system of equations
   call tridagr(ar,br,cr,knownr,potr)

   ! communicate guard cell values
   call comm_r_sweep(potr)

   !  go from having all J values in local memory to having
   !  all K values in local memory to do the vertical ADI sweep
   call swap_r_z(potr,potz)

   
   !  Vertical ADI sweep
   
   do K = 1, numz
      bz(K) = bzb(K) + dtt
   enddo

   do L = 1, numphi_dd
      do J = 1, numr_dd_z
         elambdaz(J,L) = elambdazb(J,L) + dtt
      enddo
   enddo

   knownz = 0.0

   if( iam_on_bottom ) then
      ! now the special case is independent of isym and exists for
      ! the global radial index J = 2.  We have swapped having all
      ! J values in local memory to having all K values in local 
      ! memory with J block distributed across numz_procs so pe's
      ! on the bottom of the pe grid hold the special case radial
      ! index.  The special case is that knownz with J of 2 doesn't
      ! include the influence of potz with J = 1.  Do this to impose
      ! consistent boundary condition across the z axis.
      do L = 1, numphi_dd
         do K = 2, numz-1
            do J = 3, numr_dd_z-1
               knownz(J,K,L) = - four_pi * rhoz(J,K,L) +        &
                                 elambdaz(J,L)*potz(J,K,L) -    &
                                 alphaz(J)*potz(J+1,K,L) -      &
                                 betaz(J)*potz(J-1,K,L)
            enddo
         enddo
      enddo
      ! the sepcial case
      J = 2
      do L = 1, numphi_dd
         do K = 2, numz-1
            knownz(J,K,L) = - four_pi * rhoz(J,K,L) +           &
                              elambdaz(J,L)*potz(J,K,L) -       &
                              alphaz(J)*potz(J+1,K,L)
         enddo
      enddo
   else
      do L = 1, numphi_dd
         do K = 2, numz-1
            do J = 2, numr_dd_z-1
               knownz(J,K,L) = - four_pi * rhoz(J,K,L) +        &
                                 elambdaz(J,L)*potz(J,K,L) -    &
                                 alphaz(J)*potz(J+1,K,L) -      &
                                 betaz(J)*potz(J-1,K,L)
            enddo
         enddo
      enddo
   endif

   ! add boundary potential at top and (if isym = 1)
   ! bottom of the grid to knownz
   do L = 1, numphi_dd
      do J = 2, numr_dd_z-1
         knownz(J,numz-1,L) = knownz(J,numz-1,L) +     &
                  gamma * potz(J,numz,L)
      enddo
   enddo

   if( isym == 1 ) then
      do L = 1, numphi_dd
         do J = 2, numr_dd_z-1
            knownz(J,2,L) = knownz(J,2,L) + gamma*potz(J,1,L)
         enddo
      enddo
   endif
        
   ! solve the system of equations
   call tridagz(az,bz,cz,knownz,potz)

   ! communicate guard cell values
   call comm_z_sweep(potz)

   ! rearrange the data so that K is block distributed
   ! and J is in local memory to get ready for the next
   ! iteration or to rearrange the data so that J is
   ! block distributed and l is in local memory after
   ! leaving ADI iteration loop
   call swap_z_r(potz,potr)

enddo	! END of ADI cycle

call swap_r_phi(potr,ffpot)

! inverse fourier transform the potential
do K = 1, numz_dd
   do J = 1, numr_dd
      potp(J,K,1) = ffpot(J,K,1)
      potp(J,K,2) = ffpot(J,K,numphi_by_two+1)
   enddo
enddo
index = 3
do L = 2, numphi_by_two
   do K = 1, numz_dd
      do J = 1, numr_dd
         potp(J,K,index) = ffpot(J,K,L)
      enddo
   enddo
   index = index + 2
enddo
index = 4
do L = numphi_by_two+2, numphi
   do K = 1, numz_dd
      do J = 1, numr_dd
         potp(J,K,index) = - ffpot(J,K,L)
      enddo
   enddo
   index = index + 2
enddo 
call realft(potp,numr_dd,numz_dd,numphi,-1) ! this call dies...
! Fourier transform normalization
do L = 1, numphi
   do K = 1, numz_dd
      do J = 1, numr_dd
         potp(J,K,L) = 2.0 * numphiinv * potp(J,K,L)
      enddo
   enddo
enddo

! impose boundary conditions on the potential across the
! z axis and across the equatorial plane
if( iam_on_axis ) then
   if( isym == 3 ) then
      potp(1,:,:) = potp(2,:,:)
   else
      do L= 1, numphi_by_two
         potp(1,:,L) = potp(2,:,L+numphi_by_two)
         potp(1,:,L+numphi_by_two) = potp(2,:,L)
      enddo
   endif
endif

if( iam_on_bottom ) then
   if( isym == 2 .or. isym == 3 ) then
      potp(:,1,:) = potp(:,2,:)
   endif
endif

! fill in guard cells with the potential from neighbors
call comm(potp)
 
return
end subroutine helmadi
