function rf(x,y,z)
implicit none
!  From Section 6.11, page 257 of Numerical Recipes in Fortran,
!  2nd edition.  rf computes Carlson's elliptic integral of
!  the first kind, Rf(x,y,z).  x, y, and z must be non-negative
!  and at most one argument may be zero.
!  For single precision:
!    ERRTOL = 0.08
!    TINY = 5 * Machine Underflow Limit = 6.0e-38
!    BIG = 0.2 * Machine Overflow limit = 6.0e37
!  For double precision:
!    ERRTOL = 0.0025
!    TINY = 2.0e-307
!    BIG = 2.0e307
!
!  Note that with Carlson's alogorithm, the error in the
!  return value scales as ERRTOL to the 6th power
real :: rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
#ifdef SHORT
parameter (ERRTOL=0.08,TINY=6.0e-38,BIG=6.0e37)
#else
parameter (ERRTOL=0.0025,TINY=2.0e-307,BIG=2.0e307)
#endif
parameter (THIRD=1./3.,C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
real :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,z).gt.BIG) then
   write(6,*)   'invalid arguments in rf',x,y,z
   STOP
endif 

xt=x
yt=y
zt=z
1     continue
   sqrtx=sqrt(xt)
   sqrty=sqrt(yt)
   sqrtz=sqrt(zt)
   alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
   xt=.25*(xt+alamb)
   yt=.25*(yt+alamb)
   zt=.25*(zt+alamb)
   ave=THIRD*(xt+yt+zt)
  delx=(ave-xt)/ave
   dely=(ave-yt)/ave
   delz=(ave-zt)/ave
if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
e2=delx*dely-delz**2
e3=delx*dely*delz
rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)

return
end function rf
!  (C) Copr. 1986-92 Numerical Recipes Software .
