! choose a direction for the next atom
! randomly choose theta and phi

subroutine pickdir(i,r,N,pos,nang)
   implicit none
   integer, intent(in) :: i,N,nang
   real*8, intent(inout) :: r(N,3),pos(nang,3)
   real*8 theta,phi,pi
   integer j
   
   ! define pi
   pi = 4.d0 * datan(1.d0)
   
   !random number
   do j=1,nang
      call random_number(theta)
      call random_number(phi)
      !print *, theta,phi
   
      ! sample theta from cos(theta)
      ! theta between 0 and pi
      ! phi between 0 and 2pi
      theta = 2*theta - 1 ! between 1 and -1
      theta = acos(theta) ! becomes an angle
      phi = phi*2.d0*pi
      ! define array of new possible positions
      pos(j,1) = r(i-1,1)+sin(theta)*cos(phi)
      pos(j,2) = r(i-1,2)+sin(theta)*sin(phi)
      pos(j,3) = r(i-1,3)+cos(theta)
    enddo
   
end subroutine pickdir
