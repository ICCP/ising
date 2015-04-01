   subroutine calcpot(i,V,Vtot,N,r)
   implicit none
   integer, intent(in) :: i,N
   real*8, intent(inout) :: Vtot
   real*8, intent(inout) :: V(N,N)
   real*8, intent(in) :: r(N,3)
   integer j,k
   real*8 distsq
   
   ! calculate pair-wise potentials
   if (i==3) then
      ! calculate the potential between atoms 1 and 2
      distsq = (r(1,1)-r(2,1))**2 + (r(1,2)-r(2,2))**2 + (r(1,3)-r(2,3))**2
      V(1,2) = 4*((1.d0/distsq)**6 - (1.d0/distsq)**3)
      V(2,1) = V(1,2)
      !print *, distsq,V(1,2),"1,2"
      ! should be zero unless we change sigma
   else 
      ! calculate the potential between the newest atom and the previous ones
      do j=1,i-2
         distsq = (r(j,1)-r(i-1,1))**2 + (r(j,2)-r(i-1,2))**2 + (r(j,3)-r(i-1,3))**2
	 V(i-1,j) = 4*((1.d0/distsq)**6 - (1.d0/distsq)**3)
	 V(j,i-1) = V(i-1,j)
      enddo 
   end if 
   
   ! calculate the total potential energy of the system
   ! this is for the old configuration - before the ith atom is added
   if (i==3) then
      Vtot = V(1,2)
   else
      do j=1,i-2
         Vtot = Vtot + V(i-1,j)
      enddo 
   end if 
   
   end subroutine calcpot
