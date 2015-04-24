function gammln(xx)
implicit none
!  gammln retruns the natural logarithm of Gamma(xx).  From
!  Section 6.1, p 207 of Numerical Recipes in Fortran, 2nd ed.
!  xx > 0
!  On the T3E the default kind for real is 64 bit (double precision)
!  therefor declare floating variables in gammln to be real not double
!  In the Numerical Recipes version, they had used double precision for
!  float variables private to the gammln function
real :: gammln,xx
integer :: j
real :: ser,stp,tmp,x,y,cof(6)
SAVE cof,stp
data cof,stp/76.18009172947146e0,-86.50532032941677e0,             &
24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2,     &
-.5395239384953e-5,2.5066282746310005e0/
x=xx
y=x
tmp=x+5.5e0
tmp=(x+0.5e0)*log(tmp)-tmp
ser=1.000000000190015e0
do j= 1, 6
  y=y+1.e0
  ser=ser+cof(j)/y
enddo
gammln=tmp+log(stp*ser/x)

return
end function gammln
!  (C) Copr. 1986-92 Numerical Recipes Software .
