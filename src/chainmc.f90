! main program for advanced monte carlo
! creates a chain polymer
! in 3 dimensions

program chainmc
   implicit none
   integer N,i,nangles,j,k
   real*8 T,pos,Vtot
   real*8, allocatable, dimension(:,:) :: r
   real*8, allocatable, dimension(:,:) :: V
   real*8, allocatable, dimension(:,:) :: temppos
   
   ! seed random number generator
   !call random_seed()
   call init_random_seed()
   
   ! number atoms in the chain
   N=10000
   
   ! Temperature
   T = 2.d0
   
   ! allocate memory for the positions in the chain
   allocate (r(N,3))
   ! allocate memory for the potentials between pairs
   allocate (V(N,N))
   ! initilize to zero
   r=0.d0
   V=0.d0
   
   ! set positions of the first two atoms (0,0,0) and (1,0,0)
   r(1,1)=0
   r(1,2)=0
   r(1,3)=0
   r(2,1)=1
   r(2,2)=0
   r(2,3)=0
   
   ! add the rest of the atoms
   ! define a number of sample points
   nangles = 90
   do i=3,N
      ! calculate the energy of the previous configuration
      call calcpot(i,V,Vtot,N,r)
      
      ! pick new direction, this will fill positions (r)
      ! first pick possible directions
      allocate (temppos(nangles,3))
      call pickdir(i,r,N,temppos,nangles)
      
      ! then calculate the energy of each configuration
      ! use metropolis algorithm to decide "best" configuration
      call pickconfig(i,r,N,T,temppos,Vtot,nangles)
      deallocate(temppos)
   enddo
   
   ! check the lengths of all adjacent segements
   do i=1,N-1
      pos = 0
      pos = (r(i,1)-r(i+1,1))**2 + (r(i,2)-r(i+1,2))**2 + (r(i,3)-r(i+1,3))**2
      pos = sqrt(pos)
   enddo 
   
   ! calculate the size of the system
   call calcsize(r,N)
   
   ! prepare to run program many times to get many chains
   ! for some sort of error
   deallocate(r)
   deallocate(V)


end program chainmc
