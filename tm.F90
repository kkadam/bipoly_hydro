!****************************************************************
!*
!*  TM
!*
!****************************************************************
subroutine tm(tmr)
implicit none
include 'runhydro.h'
include 'pot.h'
!****************************************************************
!*
!   tm tabulates the Green function that gets convolved with
!   the density field to evaluate the potential at the top 
!   of the grid.  Initially implemented by Howard Cohl, 
!   see his thesis for discussion and original hpf source code.
!   Also, see comment section of sm to see description of
!   parameters that determine the accuracy of the boundary
!   solution calculated with this Green function.
!*
!****************************************************************
!*
!*  Subroutine Arguments

       real, dimension(numr_dd,numz_dd,numr,mmax) :: tmr

!*
!****************************************************************
!*
!* Global Variables

real, dimension(numr_dd) :: rhf, r, rhfinv, rinv
real, dimension(numz_dd) :: zhf
real, dimension(numphi) :: phi
common /grid/ rhf, r, rhfinv, rinv, zhf, phi

real, dimension(numr) :: rhf_g, r_g, rhfinv_g, rinv_g
real, dimension(numz) :: zhf_g
common /global_grid/ rhf_g, r_g, rhfinv_g, rinv_g, zhf_g

real :: gamma, piinv, four_pi
common /pot_constants/ gamma, piinv, four_pi

integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition 

!*
!****************************************************************
!*
!*  Local Variables

real, dimension(mmax) :: qp, qm

real, dimension(mmax) :: nu, coefh

real, dimension(hypr_upr_bnd+1,mmax) :: dcoefh

real :: zB, mm, index, ir_index, sum, gammln_half
 
real :: aa, bb, cc, coef, xm, xp, mum, mup, lam, lap

real :: Kmum, Kmup, Emum, Emup, coefh_m

real :: elle, ellf, gammln

integer :: jl ! radial index for local portion of pe's data
integer :: kl ! vertical index for local portion of pe's data
integer :: jg ! global radial index
integer :: m  ! azimuthal mode number      
integer :: ir ! hypergeometric series loop index
 
!*
!*
!****************************************************************
!  initialize the local variables
qp = 0.0  
qm = 0.0  
nu = 0.0  
coefh = 0.0  
dcoefh = 0.0  
mm = 0.0  
index = 0.0
ir_index = 0.0  
sum = 0.0  
aa = 0.0  
bb = 0.0  
cc = 0.0  
coef = 0.0
xm = 0.0  
xp = 0.0  
mum = 0.0 
mup = 0.0  
lam = 0.0  
lap = 0.0  
Kmum = 0.0  
Kmup = 0.0
Emum = 0.0  
Emup = 0.0  
coefh_m = 0.0  
       
zB = zhf_g(numz)

gammln_half = gammln(0.5)

index = 3.0
if( isym == 3 ) then
   do m = 3, mmax
      nu(m) = 2.0*index - 4.5
      index = index + 1.0
   enddo
else
   do m = 3, mmax
      nu(m) = index - 2.5
      index = index + 1.0
   enddo
endif

index = 1.0   
if( isym == 3 ) then
   do m = 1, mmax
      mm = 2.0*(index - 1.0)
      aa = 0.5*(mm + 1.5)
      bb = 0.5*(mm + 0.5)
      cc = mm + 1.0
      coefh(m) = gammln_half + gammln(mm+0.5) -        &
                 gammln(aa) - gammln(bb)
      ir_index = 0.0
      do ir = 0, hypr_upr_bnd
         dcoefh(ir+1,m) = gammln(aa+ir_index) +        &
                          gammln(bb+ir_index) -        &
                          gammln(cc+ir_index) -        &
                          gammln(ir_index+1.0)
         ir_index = ir_index + 1.0
      enddo
      index = index + 1.0
   enddo
else
   do m = 1, mmax
      mm = index - 1.0
      aa = 0.5*(mm + 1.5)
      bb = 0.5*(mm + 0.5)
      cc = mm + 1.0
      coefh(m) = gammln_half + gammln(mm+0.5) -        &
                 gammln(aa) - gammln(bb)
      ir_index = 0.0
      do ir = 0, hypr_upr_bnd
         dcoefh(ir+1,m) = gammln(aa+ir_index) +        &
                          gammln(bb+ir_index) -        &
                          gammln(cc+ir_index) -        &
                          gammln(ir_index+1.0)
         ir_index = ir_index + 1.0
      enddo
      index = index + 1.0
   enddo
endif

if( isym == 1 ) then          
!  set up tmr for case of no assumed symmetry
   do jl = rlwb, rupb
      do kl = zlwb, zupb
         do jg = 2, numr-1
            coef = sqrt(rhf(jl)*rhfinv_g(jg))*         &
                   piinv
            xm = 0.5*rhfinv(jl)*rhfinv_g(jg)*          &
                 ( (zB-zhf(kl))*(zB-zhf(kl)) +         &
                   rhf_g(jg)*rhf_g(jg) +               &
                   rhf(jl)*rhf(jl) )
            if( xm < hypr_cutoff ) then
               mum = sqrt(2.0 / (1.0 + xm))
               lam = sqrt(2.0*(1.0 + xm))
               Kmum = ellf(mum)
               Emum = elle(mum)
               qm(1) = Kmum * mum
               qm(2) = xm * mum * Kmum - lam * Emum
               do m = 3, mmax
                  qm(m) = (2.0*nu(m)+1.0)/(nu(m)+1.0)*xm*qm(m-1) -      &
                          nu(m)/(nu(m)+1.0)*qm(m-2)
               enddo
            else
               index = 0.0
               do m = 1, mmax
                  coefh_m = exp( coefh(m) - (index + 0.5)*              &
                                alog(2.0*xm) )
                  sum = 0.0 
                  ir_index = 0.0
                  do ir = 0, hypr_upr_bnd
                     sum = sum + exp( dcoefh(ir+1,m) -                  &
                                  2.0*ir_index*alog(xm) )
                     ir_index = ir_index + 1.0
                  enddo
                  qm(m) = coefh_m * sum
                  index = index + 1.0
               enddo
            endif
            do m = 1, mmax
               tmr(jl,kl,jg,m) = coef * qm(m)
            enddo
         enddo
      enddo
   enddo
else if( isym == 2) then
!  set up tmr for case of equatorial symmetry
   do jl = rlwb, rupb
      do kl = zlwb, zupb
         do jg = 2, numr-1
            coef = sqrt(rhf(jl)*rhfinv_g(jg)) * piinv
            xp = 0.5*rhfinv(jl)*rhfinv_g(jg)*           &
                  ( (zB+zhf(kl))*(zB+zhf(kl)) +         &
                    rhf(jl)*rhf(jl) +                   &
                    rhf_g(jg)*rhf_g(jg) )
            xm = 0.5*rhfinv(jl)*rhfinv_g(jg) *          &
                 ( (zB-zhf(kl))*(zB-zhf(kl)) +          &
                   rhf(jl)*rhf(jl) +                    &
                   rhf_g(jg)*rhf_g(jg) )
            if( xp < hypr_cutoff ) then
               mup = sqrt(2.0/(1.0+xp))
               lap = sqrt(2.0*(1.0+xp))
               Kmup = ellf(mup)
               Emup = elle(mup)
               qp(1) = Kmup*mup
               qp(2) = xp*mup*Kmup - lap*Emup
               do m = 3, mmax
                  qp(m) = (2.0*nu(m)+1.0)/(nu(m)+1.0)*xp*qp(m-1) -   &
                          nu(m)/(nu(m)+1.0)*qp(m-2)
               enddo
            else
               index = 0.0
               do m = 1, mmax
                  coefh_m = exp(coefh(m) - (index+0.5)*      &
                              alog(2.0*xp))
                  sum = 0.0 
                  ir_index = 0.0
                  do ir = 0, hypr_upr_bnd
                     sum = sum + exp(dcoefh(ir+1,m) -        &
                                  2.0*ir_index*alog(xp))
                     ir_index = ir_index + 1.0
                  enddo
                  qp(m) = coefh_m * sum
                  index = index + 1.0
               enddo
            endif
            if( xm < hypr_cutoff ) then
               mum = sqrt(2.0/(1.0+xm))
               lam = sqrt(2.0*(1.0+xm))
               Kmum = ellf(mum)
               Emum = elle(mum)
               qm(1) = Kmum*mum
               qm(2) = xm*mum*Kmum - lam*Emum
               do m = 3, mmax
                  qm(m) = (2.0*nu(m)+1.0)/(nu(m)+1.0)*xm*qm(m-1) -     &
                          nu(m)/(nu(m)+1.0)*qm(m-2)
               enddo
            else
               index = 0.0
               do m = 1, mmax
                  coefh_m = exp(coefh(m) - (index+0.5)*       &
                              alog(2.0*xm))
                  ir_index = 0.0
                  sum = 0.0
                  do ir = 0, hypr_upr_bnd
                     sum = sum + exp(dcoefh(ir+1,m) -         &
                                  2.0*ir_index*alog(xm))
                     ir_index = ir_index + 1.0
                  enddo
                  qm(m) = coefh_m * sum
                  index = index + 1.0
               enddo
            endif
            do m = 1, mmax
               tmr(jl,kl,jg,m) = coef*(qp(m)+qm(m))
            enddo
         enddo
      enddo
   enddo
else if( isym == 3 ) then
!  set up tmr for the case of pi + equatorial symmetry
   do jl = rlwb, rupb
      do kl = zlwb, zupb
         do jg = 2, numr-1
            coef = sqrt(rhf(jl)*rhfinv_g(jg))* piinv
            xp = 0.5*rhfinv(jl)*rhfinv_g(jg)*       &
                 ( (zB+zhf(kl))*(zB+zhf(kl)) +      &
                   rhf(jl)*rhf(jl) +                &
                   rhf_g(jg)*rhf_g(jg) )
            xm = 0.5*rhfinv(jl)*rhfinv_g(jg)*       &
                 ( (zB-zhf(kl))*(zB-zhf(kl)) +      &
                   rhf(jl)*rhf(jl) +                &
                   rhf_g(jg)*rhf_g(jg) )
            if( xp < hypr_cutoff ) then
               mup = sqrt(2.0/(1.0+xp))
               lap = sqrt(2.0*(1.0+xp))
               Kmup = ellf(mup)
               Emup = elle(mup)
               qp(1) = Kmup*mup
               qp(2) = (four_thirds*xp*xp-one_third)*mup*Kmup -             &
                       four_thirds*xp*lap*Emup
               do m = 3, mmax
                  qp(m) = qp(m-1)*((2.0*nu(m)+3.0)*(2.0*nu(m)+1.0)*xp*xp/   &
                          ((nu(m)+2.0)*(nu(m)+1.0)) - (2.0*nu(m)+3.0)*      &
                          nu(m)*nu(m)/                                      &
                          ((2.0*nu(m)-1.0)*(nu(m)+2.0)*(nu(m)+1.0)) -       &
                          (nu(m)+1.0)/(nu(m)+2.0)) -                        &
                          qp(m-2)*(2.0*nu(m)+3.0)*(nu(m)-1.0)*nu(m)/        &
                          ((2.0*nu(m)-1.0)*(nu(m)+2.0)*(nu(m)+1.0))
               enddo
            else
               index = 1.0
               do m = 1, mmax
                  mm = 2.0*(index-1.0)
                  coefh_m = exp(coefh(m) - (mm+0.5)*        &
                              alog(2.0*xp))
                  sum = 0.0
                  ir_index = 0.0
                  do ir = 0, hypr_upr_bnd
                     sum = sum + exp(dcoefh(ir+1,m) -       &
                                  2.0*ir_index*alog(xp))
                     ir_index = ir_index + 1.0
                  enddo
                  qp(m) = coefh_m * sum
                  index = index + 1.0
               enddo
            endif
            if( xm < hypr_cutoff ) then
               mum = sqrt(2.0/(1.0+xm))
               lam = sqrt(2.0*(1.0+xm))
               Kmum = ellf(mum)
               Emum = elle(mum)
               qm(1) = Kmum*mum
               qm(2) = (four_thirds*xm*xm-one_third)*mum*Kmum -             &
                       four_thirds*xm*lam*Emum
               do m = 3, mmax
                  qm(m) = qm(m-1)*((2.0*nu(m)+3.0)*(2.0*nu(m)+1.0)*xm*xm/   &
                          ((nu(m)+2.0)*(nu(m)+1.0)) - (2.0*nu(m)+3.0)*      &
                          nu(m)*nu(m)/                                      &
                          ((2.0*nu(m)-1.0)*(nu(m)+2.0)*(nu(m)+1.0)) -       &
                          (nu(m)+1.0)/(nu(m)+2.0)) -                        &
                          qm(m-2)*(2.0*nu(m)+3.0)*(nu(m)-1.0)*nu(m)/        &
                          ((2.0*nu(m)-1.0)*(nu(m)+2.0)*(nu(m)+1.0))
               enddo
            else
               index = 1.0
               do m = 1, mmax
                  mm = 2.0*(index-1.0)
                  coefh_m = exp(coefh(m) - (mm+0.5)*        &
                              alog(2.0*xm))
                  sum = 0.0
                  ir_index = 0.0
                  do ir = 0, hypr_upr_bnd
                     sum = sum + exp(dcoefh(ir+1,m) -       &
                                  2.0*ir_index*alog(xm))
                     ir_index = ir_index + 1.0
                  enddo
                  qm(m) = coefh_m * sum
                  index = index + 1.0
               enddo
            endif
            do m = 1, mmax
               tmr(jl,kl,jg,m) = coef * (qp(m) + qm(m))
            enddo
         enddo
      enddo
   enddo
endif

return
end subroutine tm
