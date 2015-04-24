!**************************************************************
!*
!*  STATE
!*
!**************************************************************
subroutine state
implicit none
include 'runhydro.h'
!**************************************************************
!*
!    Using an ideal gas equation of state, assign the
!    pressure as a function of the density and internal
!    energy per unit mass
!
!      p  =  (gamma - 1) tau ** gamma
!
!    where tau = (rho * eps) ** 1 / gamma
!*
!**************************************************************
!*
!*   Global variables

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real :: pin, gamma, kappa1, kappa2, gammainv
common /polytrope/ pin, gamma, kappa1, kappa2, gammainv

!*
!**************************************************************
!*
!*   Local variables

!*
!**************************************************************

p = (gamma - 1.0) * tau**gamma

return
end subroutine state
