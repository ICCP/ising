   subroutine pickconfig(i,r,N,T,pos,Vtot,nang)
   implicit none
   integer, intent(in) :: i,N,nang
   real*8, intent(in) :: T,pos(nang,3),Vtot
   real*8, intent(inout) :: r(N,3)
   real*8 Z,Vtemp(20),dtemp,rand,prob
   integer j,k,count,pickk,l
   
   ! calculate potential for each possible configuration
   ! for the ith atom
   do k=1,nang
      Vtemp(k) = 0
      do j=1,i-1
         dtemp = 0
         dtemp = (r(j,1)-pos(k,1))**2 + (r(j,2)-pos(k,2))**2 + (r(j,3)-pos(k,3))**2
         Vtemp(k) = Vtemp(k) + 4*((1.d0/dtemp)**6 - (1.d0/dtemp)**3)
      enddo 
      Vtemp(k) = Vtemp(k) !+ Vtot
   enddo 
   
   ! calculate partition function
   Z = 0
   do l=1,nang
      Z = Z + exp(-Vtemp(l)/T)
   enddo 
   
   ! generate random number 
   call random_number(rand)
   ! pick new position based on random number
   count = 1
   prob = 0
   do while (count < nang+1)
      prob = prob + exp(-Vtemp(count)/T)/Z
      if (rand <= prob) then
         pickk = count
	 count = 1000
      else
         count = count + 1
      end if 
   enddo 
   
   ! define the new position
   r(i,1) = pos(pickk,1)
   r(i,2) = pos(pickk,2)
   r(i,3) = pos(pickk,3)   
   
   end subroutine pickconfig
