!*************************************************************************
!*
!*  OUTPUT_FILES
!*
!*************************************************************************
subroutine output_files(frnum,frnum_init)
implicit none
include 'runhydro.h'
!*************************************************************************
!*
!*  Global Variables

real, dimension(numr_dd,numz_dd,numphi) :: pot, rho
common /poisson/ pot, rho

real, dimension(numr_dd,numz_dd,numphi) :: frac1, frac2
common /fluid_frac/ frac1, frac2

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real, dimension(numr_dd,numz_dd,numphi) :: u, w, jn
common /velocities/ u, w, jn

!*
!************************************************************************* 
!*
!*   Local variables

integer :: record_number

!*
!*************************************************************************
!*
!*   Subroutine Arguments

integer :: frnum, frnum_init

!*
!*************************************************************************
! Initialize local variables

record_number = frnum - frnum_init + 1

call output_kernel(rho, record_number, 50 )

call output_kernel(frac1, record_number, 51 )

call output_kernel(frac2, record_number, 52 )

call output_kernel(tau, record_number, 53 )

call output_kernel(pot, record_number, 54 )

call output_kernel(u, record_number, 55 )

call output_kernel(w, record_number, 56 )

call output_kernel(jn, record_number, 57 )

return
end subroutine output_files
