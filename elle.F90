function elle(ak)
implicit none
!   From Section 6.11 page 261 of Numerical Recipes in Fortran,
!   2nd edition.  elle returns Legendre's elliptic integral
!   of the second kind, E(phi, kappa) using Carlson's
!   functions Rd and Rf.  See source for rf.f and rd.f.
!
!   For the hydrocode implementation only:
!   All evaluations of E(phi, kappa) are done for phi = pi / 2,
!   avoid unnecessary work by taking advantage of this fact
!
real :: elle,ak
!    USES rd,rf
real :: q,rd,rf

q = (1.0 - ak)*(1.0 + ak)
elle= rf(0.0,q,1.0) -                                         &
      ak*ak*rd(0.0,q,1.0)*0.3333333333333333

return
end function elle
!  (C) Copr. 1986-92 Numerical Recipes Software .
