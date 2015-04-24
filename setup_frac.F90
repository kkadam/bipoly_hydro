!********************************************************************************
!*
!*  SETUP_FRAC
!*
!********************************************************************************
subroutine setup(tstart, tstop, do_diag, intrvl, frnum, rho_boundary, q)
implicit none
include 'runhydro.h'
include 'mpif.h' 
include 'setup_formats.h'
!********************************************************************************
!*
!   setup does a whole bunch of things, highlights are
!   listed below:
!
!   --> read input file, fort.7, and broadcast the values to all pe's
!
!   --> initialize coordinate differentials and arrays for the volumes,
!       face areas and radii for the cell centered grid and vertex
!       centered grid
!
!   --> if an initial model call scfin3d (model_type = 0) which reads in
!       arrays for rho, s, and a.  Then set values for densmin and taumin
!       that will be used throughout the evolution and also initialize
!       the frame number counter to 1
!
!   --> if a continuation model (model_type = 1) read in the fort.12
!       and fort.13 files 
!
!   --> if a cooked up set of distributions use model_type = 2 which
!       causes debugin to be called to generate the initial conditions
!
!  5/17/2000  added support for reading in a continution file with the
!  mass fraction array present
!*
!********************************************************************************
!*
!*  Subroutine Arguments

integer :: tstart, tstop, do_diag, frnum

real :: intrvl, rho_boundary, q

!*
!********************************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

real, dimension(numr_dd,numz_dd,numphi) :: potp, rhop
common /potarrays/ potp, rhop

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd,numz_dd,numphi) :: temps, tempa
common /source_temp/ temps, tempa

real :: pi, grav, cfl_factor, den_cutoff, vlim_factor, viscosity
common /constants/ pi, grav, cfl_factor, den_cutoff, vlim_factor,    &
                   viscosity

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

real :: pin, gamma, kappa1, kappa2, gammainv 
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: kappac1, kappac2, rho_c1, rho_c2, np1, np2   !bipoly 
common /bipoly/ kappac1, kappac2, rho_c1, rho_c2, np1, np2 

integer :: isoadi, call_pot, zero_out
common /flags/ isoadi, call_pot, zero_out

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g,r_g,rhfinv_g,rinv_g,zhf_g

real, dimension(numr_dd) ::  vol_cc, volinv_cc, ar_cc, az_cc
real :: aphi_cc
common /geometry_cc/ vol_cc, volinv_cc, ar_cc, az_cc, aphi_cc

real, dimension(numr_dd) :: vol_fc_r, volinv_fc_r, ar_fc_r, az_fc_r
real :: aphi_fc_r
common /geometry_fc_r/ vol_fc_r, volinv_fc_r, ar_fc_r, az_fc_r, aphi_fc_r

real, dimension(numr_dd) :: vol_fc_z, volinv_fc_z, ar_fc_z, az_fc_z
real :: aphi_fc_z
common /geometry_fc_z/ vol_fc_z, volinv_fc_z, ar_fc_z, az_fc_z, aphi_fc_z

real, dimension(numr_dd) :: vol_fc_phi, volinv_fc_phi, ar_fc_phi, az_fc_phi
real :: aphi_fc_phi
common /geometry_fc_phi/ vol_fc_phi, volinv_fc_phi, ar_fc_phi, az_fc_phi, &
                         aphi_fc_phi

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

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
!*********************************************************************************
!*
!*  Local variables

integer :: model_type, J, K, L

#ifdef SHORT
integer*8 :: ierror
#else
integer :: ierror
#endif

real :: x, xinv, offset

integer :: record_length

integer, dimension(8) :: init_int

real, dimension(20) :: init_real           !bipoly

character(len=150) :: conts_template
character(len=150) :: conts8_file

!*
!*********************************************************************************
!*
!*  Expected format of input file called fort.7
!*
!  isym				model_type
!  1 => no symmetry		0 => an initital binary model
!  2 => equatorial symmetry	1 => continuation of evolution
!  3 => pi symmetry		2 => debugging model
!
!  tstart		tstop			do_diag
!  starting timestep #  ending timestep #	call diagnostic program
!						every do_diag timesteps
!
!  isoadi		call_pot		zero_out
!  3 => polytropic	0 => not self-		0 => don't zero out at axis
!  evolution		gravitating		# => zero out variables for
!                       1 => solve for		inner zero_out zones in j 
!                       self-gravity		for all k and l
!                       potential of fluid
!  pin                  gamma            
!  polytropic           polytropic      
!  index                exponent        
!                                                            
!
! kappa1                        kappa2                        !bipoly 
! polytropic constant for       polytropic constant for 
! 0 <= phi < pi/2               pi/2 <= phi < 3 pi/2
! 3 pi/2 <= phi < 0       
!
!                               
! kappac1                       kappac2 
! polytropic constant for       polytropic constant for 
! the core                      the core
! 0 <= phi < pi/2               pi/2 <= phi < 3 pi/2
! 3 pi/2 <= phi < 0     
!
!
! rho_c1                        rho_c2
! threshold density above       threshold density above
! which the material            which the material
! belongs to the core for       belongs to the core for
! star1                         star2
!       
!
!  np1                          np2 
!  Structural polytropic        Structural polytropic
!  indx of both the core        indx of both the envelopes
!
!               
!  omgfrm               intrvl                  scfomega
!  angular frequency    number of frames to     angular frequency of
!  of rotating frame    output per orbit        initial model
!  omgfrm = 0 =>
!  inertial frame
!  calculation
!
!  vmax                         constp
!  maximum allowed velocity     to simulate a uniform
!  for material above           pressure background, set
!  den_cutoff                   p to constp if p is less
!                               than constp
!
!  densmin                      taumin
!  floor value for density      floor value for tau = (rho eps)**1/gamma
!  shouldn't change the values of taumin and densmin during a run
!
!
!
!  boundary_condition(1)	boundary_consition(2)	boundary_condition(3)
!  bottom of grid		side of grid		top of grid
!
!  for all boundary condition entries a value of: 
!                                        1 => wall boundary
!  					 2 => free boundary
!					 3 => dirichlet boundary condition
!					      or the free flow off grid
!					      condition
!
!  rho_boundary                 q
!  value of the density         a parameter that divides the grid between
!  that marks the transition    the two stars.  It appears in the equation
!  from the star to the         of the line connecting the center of mass
!  envelope                     of the system and the maximum density point
!                               on the grid.  0 <= q <= 1
!
!
!  viscosity
!  isotropic coeeficient of the 
!  artificial viscosity
!
!*
!*********************************************************************************      
!  initialize the local variables
model_type = 0
ierror = 0
x = 0.0
xinv = 0.0
offset = 0.0
record_length = 0
init_int = 0
init_real = 0.0

!conts_template = 'input/conts/fort.12.1180000.'
 conts_template = 'input/conts/fort.12..'

!  Have root read the fort.7 file into arrays and then broadcast
!  those arrays to all other pe's
if( iam_root ) then
   open(unit=20, file='run/fort.7', &
        form='formatted' ,status='old')
   read(20,*) init_int(1:2)
   read(20,*) init_int(3:5) 
   read(20,*) init_int(6:8)

   read(20,*) init_real(1:2)          !bipoly
   read(20,*) init_real(3:4)
   read(20,*) init_real(15:16)
   read(20,*) init_real(17:18)
   read(20,*) init_real(19:20)

   read(20,*) init_real(5:7)
   read(20,*) init_real(8:9)
   read(20,*) init_real(10:11)
   read(20,*) boundary_condition(1:3) 
   read(20,*) init_real(12:14)
   close(20)
endif

call mpi_bcast(init_int,8,INT_SIZE,root,MPI_COMM_WORLD,ierror)
call mpi_bcast(init_real,20,REAL_SIZE,root,MPI_COMM_WORLD,ierror)   !bipoly
call mpi_bcast(boundary_condition,3,INT_SIZE,root,MPI_COMM_WORLD,ierror)


isym = init_int(1)
model_type = init_int(2)
tstart = init_int(3)
tstop = init_int(4)
do_diag = init_int(5)
isoadi = init_int(6)
call_pot = init_int(7)
zero_out = init_int(8)

pin = init_real(1)                !bipoly
gamma = init_real(2)
kappa1 = init_real(3)
kappa2 = init_real(4)
kappac1 = init_real(15)
kappac2 = init_real(16)
rho_c1 = init_real(17)
rho_c2 = init_real(18)
np1 = init_real(19)
np2 = init_real(20)

omega_frame = init_real(5)
intrvl = init_real(6)
scf_omega = init_real(7)
vmax = init_real(8)
constp = init_real(9)
densmin = init_real(10)
taumin = init_real(11)
rho_boundary = init_real(12)
q = init_real(13)
viscosity = init_real(14)

!  initialize the run parameters
pi = acos(-1.0)
grav = 1.0
cirp = 2.0 * pi / scf_omega
cfl_factor = 0.5
vlim_factor = 0.5
den_cutoff = 5.0e3 * densmin
gammainv = 1.0 / gamma

!  if running with pi or equaorial symmetry then the boundary
!  condition at the bottom of the grid has to be a wall
!  condition
if( isym == 2 .or. isym == 3 ) then
   boundary_condition(1) = 1
endif

! 07/24/2003 : Modified dr to reflect the numr used by Patrick to
!              create the q = 1.323 initial model. He used numr = 128
!              in the scf code which implies numr = 130 in the
!              calculation of dr.

!  setup coordinates, all pe's use their own portion of the domain
!  decomposition except for the global ( _g ) arrays which span the
!  entire grid

!  set the coordinate differentials, numr, numz and numphi come
!  from the runhydro.h header file
! >>>NOTE<<< the - 3 part comes from getting agreement from scf
!            code and hydrocode.  if numr  = 128 in the scf code
!            that translates to numr = 130 in the hydrocode
dr = 1.0 / (oldnumr - deltar_parameter)
dz = dr
dphi = 2.0 * pi * numphiinv
drinv = 1.0 / dr
dzinv = 1.0 / dz
dphiinv = 1.0 / dphi

!  define r array on every processor, use temp here to avoid 
!  coercions from integer to floating types
x = 1.0
offset = column_num*(numr_dd-2)*dr
do J = rlwb-1,rupb+1
   r(J) = offset + (x - 2.0)*dr
   x = x + 1.0
enddo

!  now define rhf from r
x = 0.5 * dr 
do J = rlwb-1,rupb+1
   rhf(J) = r(J) + x
enddo

!  and define the inverses of both arrays
where(r /= 0.0) 
  rinv = 1.0/r
elsewhere
  rinv = 0.0
endwhere

rhfinv = 1.0/rhf 

! setup the local zhf array
x = 1.0
if( isym /= 1 ) then
   offset = row_num * (numz_dd-2) * dz
   do K = zlwb-1, zupb+1
      zhf(K) = offset + (x-1.5)*dz
      x = x + 1.0
   enddo
else
   offset = (row_num*(numz_dd-2) - numz/2)*dz
   do K = zlwb-1,zupb+1
      zhf(K) = offset + (x - 0.5)*dz
      x = x + 1.0
   enddo 
endif

! set up the azimuthal angle
x = 0.0
do L = philwb, phiupb
   phi(L) = x * dphi
   x = x + 1.0
enddo

! global radius array
x = 1.0
do J = 1, numr
   r_g(J) = (x-2.0)*dr
   x = x + 1.0
enddo

! global rhf array
x = 0.5*dr
do J = 1, numr
   rhf_g(J) = r_g(J) + x
enddo

! define the inverse arrays for r and rhf
where( r_g /= 0.0 ) 
   rinv_g = 1.0/r_g
elsewhere
   rinv_g = 0.0
endwhere

rhfinv_g = 1.0/rhf_g

! setup the global zhf array
x = 1.0
if( isym /= 1 ) then
   do K = 1, numz
      zhf_g(K) = (x-1.5)*dz
      x = x + 1.0
   enddo
else
   offset = - (numz/2) * dz
   do K = 1, numz
      zhf_g(K) = offset + (x-0.5)*dz
      x = x + 1.0
   enddo
endif

! trigonemtric arrays for cell centered and vertex centered grid
x = 0.5*dphi
do L = 1, numphi
   cos_cc(L) = cos(phi(L))
   sin_cc(L) = sin(phi(L))
   cos_vc(L) = cos(phi(L)-x)
   sin_vc(L) = sin(phi(L)-x)
enddo

!  volume and area values for cell-centered grid locations

x = dr * dz * dphi
xinv = drinv * dzinv * dphiinv
!  vol_cc is the volume of a cell centered grid cell
vol_cc = x * rhf

!  volinv_cc is the inverse of this volume
volinv_cc = xinv * rhfinv

!  ar_cc is the area of the front face perpendicular
!  to the radial unit vector
do J = rlwb-1,rupb
   ar_cc(J) = r(J+1) * dphi * dz
enddo

!  az_cc is the area of cell face perpendicular to vertical
!  unit vecotr
az_cc = dr * dphi * rhf

!  aphi_cc is the area of cell face perpendicular to 
!  the azimuthal unit vector, just a scalar
aphi_cc = dr * dz

!  volume and area arrays for radial face-centered variables
vol_fc_r = x * r

volinv_fc_r = xinv * rinv

ar_fc_r = dphi * dz * rhf

az_fc_r = dr * dphi * r

aphi_fc_r = dr * dz

! volume and area arrays for vertical face-centered variables
vol_fc_z = x * rhf

volinv_fc_z = xinv * rhfinv

do J = rlwb-1, rupb
   ar_fc_z(J) = r(J+1) * dz * dphi
enddo

az_fc_z = dr * dphi * rhf

aphi_fc_z = dr * dz

! volume and area arrays for azimuthal face-centered variables
vol_fc_phi = x * rhf

volinv_fc_phi = xinv * rhfinv

do J = rlwb-1, rupb
   ar_fc_phi(J) = r(J+1) * dz * dphi
enddo

az_fc_phi = dr * dphi * rhf

aphi_fc_phi = dr * dz

if( iam_on_axis ) then
   vol_cc(rlwb-1) = vol_cc(rlwb)
   volinv_cc(rlwb-1) = volinv_cc(rlwb)
   vol_fc_r(rlwb-1) = vol_fc_r(rlwb+1)
   volinv_fc_r(rlwb-1) = volinv_fc_r(rlwb+1)
   vol_fc_z(rlwb-1) = vol_fc_z(rlwb)
   volinv_fc_z(rlwb-1) = volinv_fc_z(rlwb)
   vol_fc_phi(rlwb-1) = vol_fc_phi(rlwb)
   volinv_fc_phi(rlwb-1) = volinv_fc_phi(rlwb)
endif

if( iam_root ) then
   ! write out information about the model we are
   ! about to evolve
   if( model_type == 0 ) then
      write(6,110)
   else if( model_type == 1 ) then
      write(6,120)
   endif

   if( isym == 1 ) then
      write(6,130) tstart, tstop
   else if( isym == 2 ) then
      write(6,140) tstart, tstop
   else if( isym == 3 ) then
      write(6,150) tstart, tstop
   endif

   write(6,160) numr, numz, numphi
 
   if( boundary_condition(1) == 1 ) then
      write(6,*) 'Boundary Condition at Bottom of Grid is WALL'
   else if( boundary_condition(1) == 2 ) then
      write(6,*) 'Boundary Condition at Bottom of Grid is FREE'
   else if( boundary_condition(1) == 3 ) then
      write(6,*) 'Boundary Condition at Bottom of Grid is DIRICHLET'
   endif
   write(6,*)

   if( boundary_condition(2) == 1 ) then
      write(6,*) 'Boundary Condition at Side of Grid is WALL'
   else if( boundary_condition(2) == 2 ) then
      write(6,*) 'Boundary Condition at Side of Grid is FREE'
   else if( boundary_condition(2) == 3 ) then
      write(6,*) 'Boundary Condition at Side of Grid is DIRICHLET'
   endif
   write(6,*)

   if( boundary_condition(3) == 1 ) then
      write(6,*) 'Boundary Condition at Top of Grid is WALL'
   else if( boundary_condition(3) == 2 ) then
      write(6,*) 'Boundary Condition at Top of Grid is FREE'
   else if( boundary_condition(3) == 3 ) then
      write(6,*) 'Boundary Condition at Top of Grid is DIRICHLET'
   endif
   write(6,*)

   write(6,*) 'Zero out from j = 1 to j = ',zero_out
   write(6,*)

   write(6,170) pin, gamma, kappa1, kappa2

   write(6,180) densmin, taumin, constp
 
   write(6,190) cfl_factor

   write(6,200) den_cutoff, vlim_factor, vmax

   write(6,210) intrvl, cirp, do_diag

   write(6,220) rho_boundary, q

   write(6,225) viscosity

endif

if( model_type == 0 ) then
   ! perfrom initialization specific to a 3-d  scf model
   ! scfin reads in the following three arrays:
   !  mass density
   !  angular momentum density
   !  radial momentum density
   ! and then constructs the internal energy per unit mass
   ! velocities are initialized with the call to vel
   call scfin

   ! some housekeeping initializations
   frnum = 1000
   time = 0.0

else if(  model_type == 1 )  then
   ! perform initialization specfic to a continuation of
   ! a previous  evolution

! 10/02/2002 : Modified the write statements for fort.8

! create the filenames for the files every processor is going to read in
    if (iam .lt. 10) then
        write(conts8_file,'(a,i1)') trim(conts_template),iam
    elseif ( (iam .ge. 10) .and. (iam .le. 99) ) then
        write(conts8_file,'(a,i2)') trim(conts_template),iam
    else
        write(conts8_file,'(a,i3)') trim(conts_template),iam
    endif


!    if (iam < 10) then
!        write(conts8_file,'(a,i1)') trim(conts_template),iam
!    else
!        write(conts8_file,'(a,i2)') trim(conts_template),iam
!    endif

    open(unit=8,file=trim(conts8_file),form='unformatted', status='old')
    read(8) s, t, a, rho, tau, pot, tempa, frac1, frac2
    close(8)
   
   if( iam_root ) then
     !open(unit=9, file='input/conts/fort.13.1180000', &
 open(unit=9, file='input/conts/fort.13.', &
           form='unformatted', status='old')
      read(9) time, frnum, phi_com, R_com, com, v_com, a_com, cylin_a_com,  &
              cylin_v_com, delt

! 02/08/2003 : To check the values being read in
      write(6,*) 'fort.13 values -'
      write(6,*) time, frnum, phi_com, R_com, com, v_com,             &
                 cylin_a_com, a_com, cylin_v_com, delt
      close(9)
   endif

   call mpi_bcast(time, 1, REAL_SIZE, root, MPI_COMM_WORLD, ierror)
   call mpi_bcast(frnum, 1, INT_SIZE, root, MPI_COMM_WORLD, ierror)

else if( model_type == 2 ) then
   call debugin
   frnum = 1000
   time = 0.0
else
   ! unsupported model type
   if( iam_root ) then
      write(6,100) model_type
   endif
   stop 
endif

!  initialize the pressure
call state  

!   if self gravitating system do setup for the
!   potential solver and get the initial potential
if( call_pot == 1 ) then
   if( iam_root ) then
      write(6,230) 
   endif
   call potsetup
   call potential_solver
endif

!   fill in the velocity arrays so delta can use them
!   to find the initial timestep
call vel

return
end subroutine setup
