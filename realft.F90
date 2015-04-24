!-------------------------------------------------------------
!     REALFT
!-------------------------------------------------------------

      subroutine realft(data,nx,ny,nz,isign)

!-------------------------------------------------------------
implicit none
integer :: nx, ny, nz, isign

real,dimension(nx,ny,nz)::data

real,dimension(nx,ny)::h1i,h1r,h2i,h2r

#ifdef SHORT
real*8 :: theta, wi, wpi, wpr, wr, wtemp
#else
real :: theta, wi, wpi, wpr, wr, wtemp
#endif

real :: c1, c2, wis, wrs

integer :: i, i1, i2, i3, i4

integer :: n2p3

#ifdef SHORT
theta = 3.141592653589793 / dble(nz/2)
#else
theta = 3.141592653589793 / (nz/2)
#endif

c1=0.5
if( isign == 1 ) then
   c2=-0.5
   call four1(data,nx,ny,nz/2,+1)
else
   c2=0.5
   theta=-theta
endif
wpr=-2.0*sin(0.5*theta)**2
wpi=sin(theta)
wr=1.0+wpr
wi=wpi
n2p3=nz+3
do i=2,nz/4
   i1=2*i-1
   i2=i1+1
   i3=n2p3-i2
   i4=i3+1
   wrs=wr
   wis=wi
   h1r=c1*(data(:,:,i1)+data(:,:,i3))
   h1i=c1*(data(:,:,i2)-data(:,:,i4))
   h2r=-c2*(data(:,:,i2)+data(:,:,i4))
   h2i=c2*(data(:,:,i1)-data(:,:,i3))
   data(:,:,i1)=h1r+wrs*h2r-wis*h2i
   data(:,:,i2)=h1i+wrs*h2i+wis*h2r
   data(:,:,i3)=h1r-wrs*h2r+wis*h2i
   data(:,:,i4)=-h1i+wrs*h2i+wis*h2r
   wtemp=wr
   wr=wr*wpr-wi*wpi+wr
   wi=wi*wpr+wtemp*wpi+wi
enddo
if( isign == 1 ) then
   h1r(:,:)=data(:,:,1)
   data(:,:,1)=h1r(:,:)+data(:,:,2)
   data(:,:,2)=h1r(:,:)-data(:,:,2)
else
   h1r(:,:)=data(:,:,1)
   data(:,:,1)=c1*(h1r(:,:)+data(:,:,2))
   data(:,:,2)=c1*(h1r(:,:)-data(:,:,2))
   call four1(data,nx,ny,nz/2,-1)
endif

return
end subroutine realft

!-------------------------------------------------------------
!     FOUR1
!-------------------------------------------------------------

      subroutine four1(data,nx,ny,nnz,isign)

!-------------------------------------------------------------
implicit none

real,dimension(nx,ny,2*nnz)::data

real,dimension(nx,ny)::tempi,tempr

#ifdef SHORT
real*8 :: theta,wi,wpi,wpr,wr,wtemp
#else
real::theta,wi,wpi,wpr,wr,wtemp
#endif

integer:: nx,ny,nnz,isign,i,istep,j,m,mmax,nz

nz=2*nnz
j=1
do i=1,nz,2
   if( j > i )then
      tempr=data(:,:,j)
      tempi=data(:,:,j+1)
      data(:,:,j)=data(:,:,i)
      data(:,:,j+1)=data(:,:,i+1)
      data(:,:,i)=tempr
      data(:,:,i+1)=tempi
   endif
   m=nz/2

 1      continue
   if( (m >= 2) .and. (j > m) ) then
      j=j-m
      m=m/2
      goto 1
   endif
   j=j+m
enddo
mmax=2

 2    continue
if( nz > mmax) then
   istep=2*mmax
   theta=6.28318530717959/(isign*mmax)
   wpr=-2.*sin(0.5*theta)**2
   wpi=sin(theta)
   wr=1.
   wi=0.
   do m=1,mmax,2
      do i=m,nz,istep
         j=i+mmax
         tempr=wr*data(:,:,j)-wi*data(:,:,j+1)
         tempi=wr*data(:,:,j+1)+wi*data(:,:,j)
         data(:,:,j)=data(:,:,i)-tempr
         data(:,:,j+1)=data(:,:,i+1)-tempi
         data(:,:,i)=data(:,:,i)+tempr
         data(:,:,i+1)=data(:,:,i+1)+tempi
      enddo
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
   enddo
   mmax=istep
   goto 2
endif

return
end subroutine four1
