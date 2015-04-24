function rd(x,y,z)
implicit none
!  From Section 6.11, page 257 of Numerical Recipes in Fortran, 2nd ed.
!  Compute Carlson's ellitpic integral of the second kind, Rd(x, y, z).
!  x and y must be nonnegative and at most one of them can be zero 
!  z  must be positive
!  TINY must be at least twice (Machine overflow limit) ^ -2/3
!  BIG must be at most 0.1 ERRTOL (Machine underflow limit) ^ -2/3  
!  single precision:
!    ERRTOL =  0.05
!    TINY = 1.0e-25
!    BIG = 4.5e21
!  double precision:
!    ERRTOL: 0.0015
!    TINY = 1.0e-204
!    BIG = 1.0e201

real :: rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6

#ifdef SHORT
parameter (ERRTOL=0.05,TINY=1.0e-25,BIG=4.5e21)
#else
parameter (ERRTOL=0.0015,TINY=1.0e-204,BIG=1.0e201)
#endif
parameter (C1=3./14.,C2=1./6.,C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
real :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,         &
        sqrtz,sum,xt,yt,zt
if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG) then
   write(6,*) 'Invalid arguments in rd',x,y,z
   stop 
endif
xt=x
yt=y
zt=z
sum=0.
fac=1.
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=.2*(xt+yt+3.*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
ea=delx*dely
eb=delz*delz
ec=ea-eb
ed=ea-6.*eb
ee=ed+ec+ec
rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*      &
      ec+delz*C4*ea)))/(ave*sqrt(ave))

return
end function rd
!  (C) Copr. 1986-92 Numerical Recipes Software .
