!*********************************************************************
!*
!*      ADVECT_FRAC
!*
!*********************************************************************
subroutine advect(advect_time)
implicit none
include 'runhydro.h'
!*********************************************************************
!*
!   Advect takes conserved quantities and, depending on the flow
!   field, moves (that is advects) material from one cell to 
!   neighboring cells.
!
!   Consider the generic advection equation cast in conservative form
!
!      d Q                    ( Q Area(i) v(i) )
!      ---  =  SUM(i = 1, 6) ------------------- 
!      d t                         Volume 
!                           
!    > where d is acutally a partial derivative
!    > Q is the density of some conserved quantity in the cell 
!      of interest
!    > Area v is a vector that gives the volume swept out
!      at the interfaces of the cell by the flow in unit time
!    > Volume is the volume of the cell in quesiton 
!
!    calculate the right hand side by the following algorithm,
!    
!    --> calculate slopes of each conserved density in each cell
!        in each of the three dimensions
!
!    --> calculate product of advection areas and velocities
!
!    --> using the slopes form a linear approximation to Q in 
!        each cell to calculate the value of Q at the interface
!        between cells (where the advection velocities are
!        defined) and depending on the sign of the advection velocity
!        add or sutract the terms of the above sum times dt
!
!    the algorithm is implemented to be conservative so that a conserved
!    quantity that leaves one cell has to appear in another cell to 
!    machine accuracy per timestep.
!
!    >>NOTE<< Nothing in the hydrocode algorithm depends on having 
!             tau(n + 1/2A), therefore on the first call to advect
!             in a timestep cycle tau is not altered
!
!    Modified 5/17/2000 to make a call to a seperate routine to
!     make a different subroutine call on the second call to
!     advect for the advection of the mass density.
!     The new call keeps track of the mass fraction in each cell
!     due to one of two fluid components.  Only need to change this
!     on the second call in a timestep as it is the only time the
!     change in the density field is remembered.
!
!*
!*********************************************************************
!*
!*  Subroutine Arguments

integer :: advect_time

!*
!*********************************************************************
!*
!*  Global variables

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd,numz_dd,numphi) :: s1, t1, a1, rho1
common /save/ s1, t1, a1, rho1

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr_dd) :: vol_cc, volinv_cc, ar_cc, az_cc
real :: aphi_cc
common /geometry_cc/ vol_cc, volinv_cc, ar_cc, az_cc, aphi_cc

real, dimension(numr_dd) :: vol_fc_r, volinv_fc_r, ar_fc_r, az_fc_r
real :: aphi_fc_r
common /geometry_fc_r/ vol_fc_r, volinv_fc_r, ar_fc_r, az_fc_r, &
                       aphi_fc_r

real, dimension(numr_dd) :: vol_fc_z, volinv_fc_z, ar_fc_z, az_fc_z
real :: aphi_fc_z
common /geometry_fc_z/ vol_fc_z, volinv_fc_z, ar_fc_z, az_fc_z, &
                       aphi_fc_z

real, dimension(numr_dd) :: vol_fc_phi, volinv_fc_phi, ar_fc_phi,  &
                            az_fc_phi
real :: aphi_fc_phi
common /geometry_fc_phi/ vol_fc_phi, volinv_fc_phi, ar_fc_phi,     &
                         az_fc_phi, aphi_fc_phi

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

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
!*********************************************************************
!*
!*  Local variables

integer :: J, K, L 

real, dimension(numr_dd,numz_dd,numphi) :: retflx, rettot

real, dimension(numr_dd,numz_dd,numphi) :: transu, transv, transw

real, dimension(numr_dd,numz_dd,numphi) :: qj, qk, ql

!*
!*********************************************************************
!   initialize the local variables
retflx = 0.0
rettot = 0.0
transu = 0.0
transw = 0.0
transv = 0.0
qj = 0.0
qk = 0.0
ql = 0.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Advect the Radial Momeuntum Density, S
!

! transu is the velocity at the forward face of the
! cell times the area of that face
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb-1, rupb
         transu(J,K,L) = 0.5 * ar_fc_r(J) *           &
	                 ( u(J,K,L) + u(J+1,K,L) )
      enddo
   enddo
enddo

! transw is the velocity at the top face of the cell
! times the area of the face
do L = philwb, phiupb
   do K = zlwb-1, zupb
      do J = rlwb, rupb
         transw(J,K,L) = 0.5 * az_fc_r(J) *           &
	                 ( w(J,K+1,L) + w(J-1,K+1,L) )
      enddo
   enddo
enddo

! transv is the velocity at the left face of the cell
! times the area of the face
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         transv(J,K,L) = 0.5 * aphi_fc_r *              &
	                 ( rhfinv(J) * jn(J,K,L+1) +    &
			   rhfinv(J-1) * jn(J-1,K,L+1) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      transv(J,K,phiupb) = 0.5 * aphi_fc_r *              &
                          ( rhfinv(J) * jn(J,K,philwb) +  &
			    rhfinv(J-1) * jn(J-1,K,philwb) )
   enddo
enddo

! if you are on the edge of the grid make sure that matter
! can only flow outwards
if( iam_on_edge ) then
   do L = philwb, phiupb
      do K = zlwb, zupb
	 transu(rupb,K,L) = 0.5 * ( abs(transu(rupb,K,L))  &
	                              + transu(rupb,K,L) )
      enddo
   enddo
endif

! if you are on the top of the grid make sure that material
! can only flow upwards off the gird
if( iam_on_top ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
	 transw(J,zupb,L) = 0.5 * ( abs(transw(J,zupb,L)) + &
	                                transw(J,zupb,L) )
      enddo
   enddo
endif

! if you are on the bottom of the grid (and no symmetry) make
! sure material can only flow downwards off the grid
if( iam_on_bottom .and. isym == 1 ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
	 transw(J,zlwb-1,L) = 0.5 * ( abs(transw(J,zlwb-1,L)) -  &
	                                  transw(J,zlwb-1,L) )
      enddo
   enddo
endif

! at the axis (J = 2) for the variable S transv is not defined
if( iam_on_axis ) then
   do L = philwb, phiupb
      do K = zlwb, zupb
         transv(rlwb,K,L) = 0.0
      enddo
   enddo
endif

!   advect s
if( advect_time == 1 ) then
   call slope(s, qj, qk, ql)
   call flux(s,qj,transu,dt,volinv_fc_r,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(s,qk,transw,dt,volinv_fc_r,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(s,ql,transv,dt,volinv_fc_r,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    s(J,K,L) = s(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
else
   call slope(s1, qj, qk, ql)
   call flux(s1,qj,transu,dt,volinv_fc_r,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(s1,qk,transw,dt,volinv_fc_r,2,retflx)
   do L =  philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
            rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(s1,ql,transv,dt,volinv_fc_r,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    s(J,K,L) = s1(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Advect the Vertical Momentum Density, T
!

! transu is the velocity at the forward cell face times
! the area of the face
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb-1, rupb
         transu(J,K,L) = 0.5 * ar_fc_z(J) *                     &
	                 ( u(J+1,K,L) + u(J+1,K-1,L) )
      enddo
   enddo
enddo

! transw is the velocity at the top face of the cell times
! the area of that face
do L = philwb, phiupb
   do K = zlwb-1, zupb
      do J = rlwb, rupb
         transw(J,K,L) = 0.5 * az_fc_z(J) *                     &
	                 ( w(J,K,L) + w(J,K+1,L) )
      enddo
   enddo
enddo

! transv is the velocity at the left face of the cell times
! the area of the face
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         transv(J,K,L) = 0.5 * aphi_fc_z * rhfinv(J) *          &
	                 ( jn(J,K,L+1) + jn(J,K-1,L+1) ) 
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      transv(J,K,phiupb) = 0.5 * aphi_fc_z * rhfinv(J) *        &
                           ( jn(J,K,philwb) + jn(J,K-1,philwb) )
   enddo
enddo

! if on the edge of the grid allow the matter to only flow outwards
! off the grid
if( iam_on_edge ) then
   do L = philwb, phiupb
      do K = zlwb, zupb
	 transu(rupb,K,L) = 0.5 * ( abs(transu(rupb,K,L)) +     &
	                                transu(rupb,K,L) )
      enddo
   enddo   
endif

! if on the top of the grid allow matter to only flow upwards
! off the top of the grid
if( iam_on_top ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
	 transw(J,zupb,L) = 0.5 * ( abs(transw(J,zupb,L)) +     &
	                                transw(J,zupb,L) )
      enddo
   enddo
endif

! if on the bottom of the grid and running with no symmetries
! only allow matter to flow downwards off the grid
if( iam_on_bottom .and. isym == 1 ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
         transw(J,zlwb-1,L) = 0.5 * ( abs(transw(J,zlwb-1,L)) -   &
	                                  transw(J,zlwb-1,L) ) 
      enddo
   enddo
endif

if( advect_time == 1 ) then
   call slope(t, qj, qk, ql)
   call flux(t,qj,transu,dt,volinv_fc_z,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(t,qk,transw,dt,volinv_fc_Z,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(t,ql,transv,dt,volinv_fc_z,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    t(J,K,L) = t(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
else
   call slope(t1, qj, qk, ql)
   call flux(t1,qj,transu,dt,volinv_fc_z,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(t1,qk,transw,dt,volinv_fc_z,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(t1,ql,transv,dt,volinv_fc_z,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    t(J,K,L) = t1(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Advect the Angular Momentum Density, A
!

! transu is the velocity at the forward cell face times the
! area of the face
do L = philwb+1, phiupb
   do K = zlwb, zupb
      do J = rlwb-1, rupb
         transu(J,K,L) = 0.5 * ar_fc_phi(J) *                &
	                 ( u(J+1,K,L) + u(J+1,K,L-1) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb-1, rupb
      transu(J,K,philwb) = 0.5 * ar_fc_phi(J) *               &
                           ( u(J+1,K,philwb) + u(J+1,K,phiupb) )
   enddo
enddo

! transw is the velocity at the top cell face times the area
! of the face
do L = philwb+1, phiupb
   do K = zlwb-1, zupb
      do J = rlwb, rupb
         transw(J,K,L) = 0.5 * az_fc_phi(J) *                 &
	                 ( w(J,K+1,L) + w(J,K+1,L-1) )
      enddo
   enddo
enddo
do K = zlwb-1, zupb
   do J = rlwb, rupb
      transw(J,K,philwb) = 0.5 * az_fc_phi(J) *               &
                           ( w(J,K+1,philwb) + w(J,K+1,phiupb) )
   enddo
enddo

! transv is the velocity at the left face of the cell times the
! area of that face
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         transv(J,K,L) = 0.5 * aphi_fc_phi * rhfinv(J) *       &
	                 ( jn(J,K,L) + jn(J,K,L+1) )
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      transv(J,K,phiupb) = 0.5 * aphi_fc_phi * rhfinv(J) *     &
                           ( jn(J,K,phiupb) + jn(J,K,philwb) )
   enddo
enddo

! if on the edge of the grid only allow the matter to flow
! outwards off the grid
if( iam_on_edge ) then
   do L = philwb, phiupb
      do K = zlwb, zupb
	 transu(rupb,K,L) = 0.5 * ( abs(transu(rupb,K,L)) +    &
	                                transu(rupb,K,L) )
      enddo
   enddo
endif

! if on the top of the grid allow matter to only flow upwards off
! the grid
if( iam_on_top ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
	 transw(J,zupb,L) = 0.5 * ( abs(transw(J,zupb,L)) +    &
	                                transw(J,zupb,L) )
      enddo
   enddo
endif

! if on the bottom of the grid and if running with no symmetries
! then only allow matter to flow downwards off the grid
if( iam_on_bottom .and. isym == 1 ) then
   do L = philwb, phiupb
      do J= rlwb, rupb
         transw(J,zlwb-1,L) = 0.5 * ( abs(transw(J,zlwb-1,L)) - &
	                                  transw(J,zlwb-1,L) )
      enddo
   enddo
endif

!   advect a:
if( advect_time == 1 ) then
   call slope(a, qj, qk, ql)
   call flux(a,qj,transu,dt,volinv_fc_phi,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(a,qk,transw,dt,volinv_fc_phi,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(a,ql,transv,dt,volinv_fc_phi,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    a(J,K,L) = a(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
else
   call slope(a1, qj, qk, ql)
   call flux(a1,qj,transu,dt,volinv_fc_phi,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(a1,qk,transw,dt,volinv_fc_phi,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(a1,ql,transv,dt,volinv_fc_phi,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    a(J,K,L) = a1(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Advect the Mass Density, rho and if this is the second
!  call to advect in a timestep cycle advect the entropy
!  tracer Tau as well
!

! transu is the velocity at the front face of a cell times
! the area of the face
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb-1, rupb
         transu(J,K,L) = ar_cc(J) * u(J+1,K,L)
      enddo
   enddo
enddo

! transw is the velocity at the top face of the cell times
! the area of that face
do L = philwb, phiupb
   do K = zlwb-1, zupb
      do J = rlwb, rupb
         transw(J,K,L) = az_cc(J) * w(J,K+1,L)
      enddo
   enddo
enddo

! transv is the velocity at the left face of the cell times
! the area of that face
do L = philwb, phiupb-1
   do K = zlwb, zupb
      do J = rlwb, rupb
         transv(J,K,L) = aphi_cc * rhfinv(J) * jn(J,K,L+1)
      enddo
   enddo
enddo
do K = zlwb, zupb
   do J = rlwb, rupb
      transv(J,K,phiupb) = aphi_cc * rhfinv(J) * jn(J,K,philwb)
   enddo
enddo

! if on the edge of the grid only allow the matter to flow
! outwards off the grid
if( iam_on_edge ) then
   do L = philwb, phiupb
      do K = zlwb, zupb
         transu(rupb,K,L) = 0.5 * ( abs(transu(rupb,K,L)) +     &
                                        transu(rupb,K,L) )
      enddo
   enddo
endif

! if on the top of the grid allow the matter to only flow
! upwards off the grid
if( iam_on_top ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
	 transw(J,zupb,L) = 0.5 * ( abs(transw(J,zupb,L)) +      &
	                                transw(J,zupb,L) )
      enddo
   enddo
endif

! if on the bottom of the grid and running with no symmetries
! then only allow the matter to flow downwards off the grid
if( iam_on_bottom .and. isym == 1 ) then
   do L = philwb, phiupb
      do J = rlwb, rupb
         transw(J,zlwb-1,L) = 0.5 * ( abs(transw(J,zlwb-1,L)) -  &
	                                  transw(J,zlwb-1,L) )
      enddo
   enddo
endif

!   advect rho:
if( advect_time == 1 ) then
   call slope(rho, qj, qk, ql)
   call flux(rho,qj,transu,dt,volinv_cc,1,retflx)
   do L= philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(rho,qk,transw,dt,volinv_cc,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(rho,ql,transv,dt,volinv_cc,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rho(J,K,L) = rho(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
else
   call slope(rho1, qj, qk, ql)
   ! keep track of fluid types
   call flux_frac(rho1,qj,qk,ql,transu,transw,transv,dt) 
   ! actually update rho due to advection
   call flux(rho1,qj,transu,dt,volinv_cc,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(rho1,qk,transw,dt,volinv_cc,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(rho1,ql,transv,dt,volinv_cc,3,retflx)
   do L= philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rho(J,K,L) = rho1(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
endif

!   advect tau if second call to advcet in timestep cycle:
if( advect_time == 2 ) then
   call slope(tau, qj, qk, ql)
   call flux(tau,qj,transu,dt,volinv_cc,1,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
            rettot(J,K,L) = retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(tau,qk,transw,dt,volinv_cc,2,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    rettot(J,K,L) = rettot(J,K,L) + retflx(J,K,L)
	 enddo
      enddo
   enddo
   call flux(tau,ql,transv,dt,volinv_cc,3,retflx)
   do L = philwb, phiupb
      do K = zlwb, zupb
         do J = rlwb, rupb
	    tau(J,K,L) = tau(J,K,L) + rettot(J,K,L) + retflx(J,K,L)
         enddo	
      enddo
   enddo
endif

call advect_bc

!   communicate guard cell values of s, t, a, rho, tau
call comm(s)
call comm(t)
call comm(rho)
call comm(frac1)
call comm(frac2)
call comm(a)
if( advect_time == 2 ) then
   call comm(tau)
endif

return
end subroutine advect
