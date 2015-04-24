!*********************************************************************
!*
!*      ADVECT_BC
!*
!*********************************************************************
subroutine advect_bc
implicit none
include 'runhydro.h'
!*********************************************************************
!*
!
!  Apply the boundary conditions after advection
!
!*
!*********************************************************************
!*
!*  Subroutine Arguments
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

real, dimension(numr_dd,numz_dd,numphi) :: p, tau
common /thermo/ p, tau

real :: densmin, taumin, vmax, constp
common /limits/ densmin, taumin, vmax, constp
 
integer :: isym
integer, dimension(3) :: boundary_condition
common /boundary_conditions/ isym, boundary_condition

logical :: iam_on_top, iam_on_bottom, iam_on_axis,                    &
           iam_on_edge, iam_root
integer :: column_num, row_num
#ifdef SHORT
integer*8 :: iam, down_neighbor, up_neighbor,                         &
             in_neighbor, out_neighbor, root,                         &
             REAL_SIZE, INT_SIZE, numprocs
#else
integer :: iam, down_neighbor, up_neighbor,                           &
           in_neighbor, out_neighbor, root,                           &
           REAL_SIZE, INT_SIZE, numprocs
#endif
integer, dimension(numr_procs,numz_procs) :: pe_grid
common /processor_grid/ iam, numprocs, iam_on_top,                    &
                        iam_on_bottom, iam_on_axis,                   &
                        iam_on_edge, down_neighbor,                   &
                        up_neighbor, in_neighbor,                     &
                        out_neighbor, root, column_num,               &
                        row_num, pe_grid, iam_root,                   &
                        REAL_SIZE, INT_SIZE
!*
!*********************************************************************
!*
!*  Local variables

integer :: J, K, L 

real :: densmininv, sum, sum1, sum2

!*
!*********************************************************************
!   initialize the local variables
densmininv = 1.0 / densmin

!   set floor values of rho and tau
do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( rho(J,K,L) <= densmin ) then
            frac1(J,K,L) = frac1(J,K,L)*rho(J,K,L)*densmininv
            frac2(J,K,L) = frac2(J,K,L)*rho(J,K,L)*densmininv
            rho(J,K,L) = densmin
         endif
      enddo
   enddo
enddo

do L = philwb, phiupb
   do K = zlwb, zupb
      do J = rlwb, rupb
         if( tau(J,K,L) <= taumin ) tau(J,K,L) = taumin
      enddo
   enddo
enddo

!   apply boundary conditions to conserved quantities
if( iam_on_axis ) then

   ! s not defined on the axis
   do L = philwb, phiupb
      do K = zlwb, zupb
         s(rlwb,K,L) = 0.0
      enddo
   enddo

   do K = zlwb, zupb
      sum = 0.0
      sum1 = 0.0
      sum2 = 0.0
      do L = philwb, phiupb
         sum = sum + rho(rlwb,K,L)
	 sum1 = sum1 + rho(rlwb,K,L)*frac1(rlwb,K,L)
	 sum2 = sum2 + rho(rlwb,K,L)*frac2(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      sum1 = sum1 * numphiinv
      sum2 = sum2 * numphiinv
      do L = philwb, phiupb
         rho(rlwb,K,L) = sum
	 frac1(rlwb,K,L) = sum1 / rho(rlwb,K,L)
	 frac2(rlwb,K,L) = sum2 / rho(rlwb,K,L)
      enddo
      sum = 0.0
      do L = philwb, phiupb
         sum = sum + a(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         a(rlwb,K,L) = sum
      enddo
      sum = 0.0
      do L = philwb, phiupb
         sum = sum + tau(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         tau(rlwb,K,L) = sum
      enddo
      sum = 0.0
      do L = philwb, phiupb
         sum = sum + t(rlwb,K,L)
      enddo
      sum = sum * numphiinv
      do L = philwb, phiupb
         t(rlwb,K,L) = sum
      enddo
   enddo

   if( isym == 3 ) then
      do L = philwb, phiupb
         do K = zlwb, zupb
            rho(rlwb-1,K,L)   = rho(rlwb,K,L)
            frac1(rlwb-1,K,L) = frac1(rlwb,K,L)
            frac2(rlwb-1,K,L) = frac2(rlwb,K,L)
            tau(rlwb-1,K,L)   = tau(rlwb,K,L)
            a(rlwb-1,K,L)     =   a(rlwb,K,L)
            s(rlwb-1,K,L)     = - s(rlwb+1,K,L)
            t(rlwb-1,K,L)     =   t(rlwb,K,L)
         enddo
      enddo
   else
      do L = 1, numphi_by_two
         do K = zlwb, zupb
            rho(rlwb-1,K,L)               = rho(rlwb,K,L+numphi_by_two)
            rho(rlwb-1,K,L+numphi_by_two) = rho(rlwb,K,L)
	 enddo
      enddo

      do L = 1, numphi_by_two
         do K = zlwb, zupb
            frac1(rlwb-1,K,L)               = frac1(rlwb,K,L+numphi_by_two)
            frac1(rlwb-1,K,L+numphi_by_two) = frac1(rlwb,K,L)
	 enddo
      enddo

      do L = 1, numphi_by_two
         do K = zlwb, zupb
            frac2(rlwb-1,K,L)               = frac2(rlwb,K,L+numphi_by_two)
            frac2(rlwb-1,K,L+numphi_by_two) = frac2(rlwb,K,L)
	 enddo
      enddo

      do L = 1, numphi_by_two
         do K = zlwb, zupb
            tau(rlwb-1,K,L)               = tau(rlwb,K,L+numphi_by_two)
            tau(rlwb-1,K,L+numphi_by_two) = tau(rlwb,K,L)
	 enddo
      enddo

      do L = 1, numphi_by_two
         do K = zlwb, zupb
            a(rlwb-1,K,L)               = a(rlwb,K,L+numphi_by_two)
            a(rlwb-1,K,L+numphi_by_two) = a(rlwb,K,L)
	 enddo
      enddo

      do L = 1, numphi_by_two
         do K = zlwb, zupb
            s(rlwb-1,K,L)               = - s(rlwb+1,K,L+numphi_by_two)
            s(rlwb-1,K,L+numphi_by_two) = - s(rlwb+1,K,L)
	 enddo
      enddo

      do L = 1, numphi_by_two
         do K = zlwb, zupb
            t(rlwb-1,K,L)               = t(rlwb,K,L+numphi_by_two)
            t(rlwb-1,K,L+numphi_by_two) = t(rlwb,K,L)
	 enddo
      enddo
            
   endif
endif

if( iam_on_bottom ) then
   if( isym == 1 ) then
      if( boundary_condition(1) == 1 ) then
         ! wall boundary condition
	 do L = philwb, phiupb
	    do J = rlwb, rupb
               rho(J,zlwb-1,L)   = rho(J,zlwb,L)
               frac1(J,zlwb-1,L) = frac1(J,zlwb,L)
               frac2(J,zlwb-1,L) = frac2(J,zlwb,L)
               tau(J,zlwb-1,L)   = tau(J,zlwb,L)
               a(J,zlwb-1,L)     =   a(J,zlwb,L)
               s(J,zlwb-1,L)     =   s(J,zlwb,L)
               t(J,zlwb-1,L)     = - t(J,zlwb+1,L)
               t(J,zlwb,L)       =   0.0
	    enddo
	 enddo
      else if( boundary_condition(1) == 2 ) then
         ! free boundary condition
	 do L = philwb, phiupb
	    do J = rlwb, rupb
               rho(J,zlwb-1,L)   = rho(J,zlwb,L)
               frac1(J,zlwb-1,L) = frac1(J,zlwb,L)
               frac2(J,zlwb-1,L) = frac2(J,zlwb,L)
               tau(J,zlwb-1,L)   = tau(J,zlwb,L)
               a(J,zlwb-1,L)     =   a(J,zlwb,L)
               s(J,zlwb-1,L)     =   s(J,zlwb,L)
               t(J,zlwb-1,L)     =   t(J,zlwb,L)
	    enddo
	 enddo
      else 
         !dirichlet bc
	 do L = philwb, phiupb
	    do J = rlwb, rupb
	       s(J,zlwb-1,L) = 0.0
	       a(J,zlwb-1,L) = 0.0
	       t(J,zlwb,L) = 0.5*(abs(t(J,zlwb,L)) -   &
                                      t(J,zlwb,L)  )
	       t(J,zlwb-1,L) = 0.0
	    enddo
	 enddo
      endif
   else
      ! if isym /= 1 enforce wall boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            rho(J,zlwb-1,L)   =   rho(J,zlwb,L)
            frac1(J,zlwb-1,L) = frac1(J,zlwb,L)
            frac2(J,zlwb-1,L) = frac2(J,zlwb,L)
            tau(J,zlwb-1,L)   =   tau(J,zlwb,L)
            a(J,zlwb-1,L)     =     a(J,zlwb,L)
            s(J,zlwb-1,L)     =     s(J,zlwb,L)
            t(J,zlwb-1,L)     =  -  t(J,zlwb+1,L)
            t(J,zlwb,L)       =      0.0
         enddo
      enddo
   endif
endif 

if( iam_on_edge ) then
   if( boundary_condition(2) == 1 ) then
      ! wall boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            rho(rupb+1,K,L)   =   rho(rupb,K,L)
            frac1(rupb+1,K,L) = frac1(rupb,K,L)
            frac2(rupb+1,K,L) = frac2(rupb,K,L)
            tau(rupb+1,K,L)   =   tau(rupb,K,L)
            a(rupb+1,K,L)     =     a(rupb,K,L)
            s(rupb+1,K,L)     =  -  s(rupb-1,K,L)
            s(rupb,K,L)       =     0.0
            t(rupb+1,K,L)     =     t(rupb,K,L)
         enddo
      enddo
   else if( boundary_condition(2) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            rho(rupb+1,K,L)   =   rho(rupb,K,L)
            frac1(rupb+1,K,L) = frac1(rupb,K,L)
            frac2(rupb+1,K,L) = frac2(rupb,K,L)
            tau(rupb+1,K,L)   =   tau(rupb,K,L)
            a(rupb+1,K,L)     =     a(rupb,K,L)
            s(rupb+1,K,L)     =     s(rupb,K,L)
            t(rupb+1,K,L)     =     t(rupb,K,L)
         enddo
      enddo
   else	
      ! dirichlet boundary condition
      do L = philwb, phiupb
         do K = zlwb, zupb
            s(rupb,K,L)   = 0.5*(abs(s(rupb,K,L)) +         &
                                     s(rupb,K,L) )
            s(rupb+1,K,L) = 0.0
            a(rupb+1,K,L) = 0.0
            t(rupb+1,K,L) = 0.0
         enddo
      enddo
   endif
endif

if( iam_on_top ) then
   if( boundary_condition(3) == 1 ) then
      ! wall boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            rho(J,zupb+1,L)   =   rho(J,zupb,L)
            frac1(J,zupb+1,L) = frac1(J,zupb,L)
            frac2(J,zupb+1,L) = frac2(J,zupb,L)
            tau(J,zupb+1,L)   =   tau(J,zupb,L)
            a(J,zupb+1,L)     =     a(J,zupb,L)
            s(J,zupb+1,L)     =     s(J,zupb,L)
            t(J,zupb+1,L)     =  -  t(J,zupb-1,L)
            t(J,zupb,L)       =     0.0
         enddo
      enddo
   else if( boundary_condition(3) == 2 ) then
      ! free boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            rho(J,zupb+1,L)   =   rho(J,zupb,L)
            frac1(J,zupb+1,L) = frac1(J,zupb,L)
            frac2(J,zupb+1,L) = frac2(J,zupb,L)
            tau(J,zupb+1,L)   =   tau(J,zupb,L)
            a(J,zupb+1,L)     =     a(J,zupb,L)
            s(J,zupb+1,L)     =     s(J,zupb,L)
            t(J,zupb+1,L)     =     t(J,zupb,L)
         enddo
      enddo
   else	
      ! dirichlet boundary condition
      do L = philwb, phiupb
         do J = rlwb, rupb
            s(J,zupb+1,L) = 0.0
            a(J,zupb+1,L) = 0.0
            t(J,zupb,L)   = 0.5*(abs(t(J,zupb,L)) +          &
                                     t(J,zupb,L))
            t(J,zupb+1,L) = 0.0
         enddo
      enddo
   endif
endif

return
end subroutine advect_bc
