!*********************************************************************
!*
!*  INITIALIZE
!*
!*********************************************************************
subroutine initialize
implicit none
include 'runhydro.h'
include 'pot.h'
!*********************************************************************
!*
!  initialize does initialization for all data in commmon blocks
!  in hydrocode.
!*
!*********************************************************************
!*
!*  Global Variables
  
real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

real, dimension(numr_dd,numz_dd,numphi) :: s, t, a
common /kinematic/ s, t, a

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

real :: phi_com, R_com, R_com_inv
real, dimension(3,3) :: com
real, dimension(3) :: v_com, a_com, delt
real, dimension(3,numphi) :: cylin_v_com, cylin_a_com
common /centerofmass/ phi_com, R_com, R_com_inv, com, v_com,    &
                      a_com, delt, cylin_v_com, cylin_a_com

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd,numz_dd,numphi) :: potp, rhop
common /potarrays/ potp, rhop

real, dimension(numr_dd,numz_dd,numphi) :: s1, t1, a1, rho1
common /save/ s1, t1, a1, rho1

real, dimension(numr_dd,numz_dd,numphi) :: temps, tempa
common /source_temp/ temps, tempa

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g, r_g, rhfinv_g, rinv_g, zhf_g

real, dimension(numr_dd) :: vol_cc, volinv_cc, ar_cc, az_cc
real :: aphi_cc
common /geometry_cc/ vol_cc, volinv_cc, ar_cc, az_cc, aphi_cc

real, dimension(numr_dd) :: vol_fc_r, volinv_fc_r, ar_fc_r, az_fc_r
real :: aphi_fc_r
common /geometry_fc_r/ vol_fc_r, volinv_fc_r, ar_fc_r, az_fc_r, aphi_fc_r

real, dimension(numr_dd) :: vol_fc_z, volinv_fc_z, ar_fc_z, az_fc_z
real :: aphi_fc_z
common /geometry_fc_z/ vol_fc_z, volinv_fc_z, ar_fc_z, az_fc_z, aphi_fc_z

real, dimension(numr_dd) :: vol_fc_phi, volinv_fc_phi, ar_fc_phi,   &
                            az_fc_phi
real :: aphi_fc_phi
common /geometry_fc_phi/ vol_fc_phi, volinv_fc_phi, ar_fc_phi,     &
                         az_fc_phi, aphi_fc_phi

real, dimension(numphi) :: cos_cc, sin_cc, cos_vc, sin_vc
common /trig/ cos_cc, sin_cc, cos_vc, sin_vc

real, dimension(numr_dd,numz_dd,numr,mmax) :: tmr
real, dimension(numr_dd,numz_dd,numr,mmax) :: bmr
real, dimension(numr_dd,numz_dd,numz,mmax) :: smz
common /green_functions/ tmr, bmr, smz

real, dimension(numphi,mmax) :: bes_cos, bes_sin
common /bessel_trig/ bes_cos, bes_sin

real, dimension(numr) :: ar, cr, alphar
real, dimension(numr, numphi_dd) ::  brb
common /ADI_R_sweep/ ar, cr, alphar, brb

real :: az, cz
real, dimension(numr_dd_z) :: alphaz, betaz
real, dimension(numz) :: bzb
real, dimension(numr_dd_z,numphi_dd) :: elambdazb
common /ADI_Z_sweep/ az, cz, alphaz, betaz, bzb, elambdazb

real :: dr, dz, dphi, drinv, dzinv, dphiinv
common /coord_differentials/ dr, dz, dphi, drinv, dzinv, dphiinv

real :: dt, time, dt_visc
integer :: tstep
common /timestep/ dt, time, dt_visc, tstep

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

real :: kappac1, kappac2, rho_c1, rho_c2, np1, np2   !bipoly 
common /bipoly/ kappac1, kappac2, rho_c1, rho_c2, np1, np2

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp

integer :: isoadi, call_pot, zero_out
common /flags/ isoadi, call_pot, zero_out

real :: pi, grav, cfl_factor, den_cutoff, vlim_factor, viscosity
common /constants/ pi, grav, cfl_factor, den_cutoff, vlim_factor,  &
                   viscosity

real :: omega_frame, cirp, scf_omega
common /rotframe/ omega_frame, cirp, scf_omega

real :: gam, piinv, four_pi
common /pot_constants/ gam, piinv, four_pi

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
pot = 0.0
rho = 0.0

frac1 = 0.0
frac2 = 0.0

s = 0.0
t = 0.0
a = 0.0

u = 0.0
w = 0.0
jn = 0.0

! 05/24/2003 : Modified for com motion
phi_com = 0.0
R_com = 0.0
com = 0.0
v_com = 0.0
a_com = 0.0
cylin_a_com = 0.0
cylin_v_com = 0.0
delt = 0.0

p = 0.0
tau = 0.0

potp = 0.0
rhop = 0.0

s1 = 0.0
t1 = 0.0
a1 = 0.0
rho1 = 0.0

temps = 0.0
tempa = 0.0

rhf = 0.0
r = 0.0
rhfinv = 0.0
rinv = 0.0
zhf = 0.0
phi = 0.0

rhf_g = 0.0
r_g = 0.0
rhfinv_g = 0.0
rinv_g = 0.0
zhf_g = 0.0

vol_cc = 0.0
volinv_cc = 0.0
ar_cc = 0.0
az_cc = 0.0
aphi_cc = 0.0

vol_fc_r = 0.0
volinv_fc_r = 0.0
ar_fc_r = 0.0
az_fc_r = 0.0
aphi_fc_r = 0.0

vol_fc_z = 0.0
volinv_fc_z = 0.0
ar_fc_z = 0.0
az_fc_z = 0.0
aphi_fc_z = 0.0

vol_fc_phi = 0.0
volinv_fc_phi = 0.0
ar_fc_phi = 0.0
az_fc_phi = 0.0
aphi_fc_phi = 0.0

cos_cc = 0.0
sin_cc = 0.0
cos_vc = 0.0
sin_vc = 0.0

tmr = 0.0
bmr = 0.0
smz = 0.0

bes_cos = 0.0
bes_sin = 0.0

ar = 0.0
cr = 0.0
alphar = 0.0
brb = 0.0

az = 0.0
cz = 0.0
alphaz = 0.0
betaz = 0.0
bzb = 0.0
elambdazb = 0.0

dr = 0.0
dz = 0.0
dphi = 0.0
drinv = 0.0
dzinv = 0.0
dphiinv = 0.0

dt = 0.0
time = 0.0
tstep = 0
dt_visc = 0.0

pin = 0.0
gamma = 0.0
kappa1 = 0.0
kappa2 = 0.0
gammainv = 0.0

kappac1 = 0.0      !bipoly
kappac2 = 0.0
rho_c1 = 0.0
rho_c2 = 0.0
np1 = 0.0
np2 = 0.0

densmin = 0.0
taumin = 0.0
vmax = 0.0
constp = 0.0

isoadi = 0
call_pot = 0
zero_out = 0

pi = 0.0
grav = 0.0
cfl_factor = 0.0
den_cutoff = 0.0
vlim_factor = 0.0

omega_frame = 0.0
cirp = 0.0
scf_omega = 0.0

gam = 0.0
piinv = 0.0
four_pi = 0.0

isym = 0
boundary_condition = 0

iam_on_top = .false.
iam_on_bottom = .false.
iam_on_axis = .false.
iam_on_edge = .false.
iam_root = .false.
iam = 0
numprocs = 0
down_neighbor = -1
up_neighbor = -1
in_neighbor = -1
out_neighbor = -1
root = 0
column_num = 0
row_num = 0
pe_grid = 0 
REAL_SIZE = 0
INT_SIZE = 0

return
end subroutine initialize
