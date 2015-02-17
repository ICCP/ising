program ising
  
  use global

  implicit none
  
  !External Functions
  real(8),external :: energy
  integer,external :: hfield

  !Variables for main routine
  integer :: ii
  integer :: tot_energy

  !output files
  open (unit = 1, file = 'energy.out',status = 'unknown')
  open (unit = 2, file = 'hfield.out',status = 'unknown')

  !initializing variables
  ii = 0                 !Iteration Variable
  energy_sum = 0         !Energy Sum
  
  !Settings
  latticeSize = 100      !Lattice Size
  iters = 10000000       !Number of iterations the program runs through
  rand_loc = 2           !Init Random Location Variable
  j = 1.d0               !magnetic coeffient
  k = 1                  !Nomarlized Boltzman Constant
  temp = 2               !Tempereature of simulation
  B = 1/(k*temp)         !constant of temp and boltman

  allocate(lattice(0:latticeSize + 1,0:latticeSize + 1))
  lattice = 0
  lattice(1:latticeSize, 1:latticeSize) = 1
  
  !main program 
  tot_energy = energy(lattice,latticeSize)!Calculate Chagnge total energy
  do ii = 0,iters
      !generate a random location 
      Call random_number(rand_loc)
      rand_loc = rand_loc*(latticeSize)
      rand_loc_int = ceiling(rand_loc)
      
      call del_energy !calculate change in energy
      call flip_bit   !Perform metropolis test
      
      !Write results to file
      write (1,*),ii, energy_sum 
      write (2,*),ii, sum(lattice)/float(latticeSize**2)

  end do
end program 
!------------------------------------------------------------------------------
subroutine del_energy!{{{

  use global, only: del_e

  integer :: spin,h1,h2,h3,h4
      spin = lattice(rand_loc_int(1),rand_loc_int(2))
      h1 = lattice(rand_loc_int(1)+1,rand_loc_int(2)+1)
      h2 = lattice(rand_loc_int(1)+1,rand_loc_int(2)-1)
      h3 = lattice(rand_loc_int(1)-1,rand_loc_int(2)+1)
      h4 = lattice(rand_loc_int(1)-1,rand_loc_int(2)-1)
      del_e = 2*spin*(h1 + h2 + h3 + h4)
end subroutine !}}}
!------------------------------------------------------------------------------
subroutine flip_bit !{{{
  use global
  real(KIND=8) :: boltzmann

    
  integer :: spin
  spin = 0
  spin = lattice(rand_loc_int(1),rand_loc_int(2))
!   if (del_e .lt. 0) then  
!       spin=lattice(rand_loc_int(1),rand_loc_int(2))
!       lattice(rand_loc_int(1),rand_loc_int(2))=-spin
!       energy_sum=energy_sum+del_e
    call random_number(rand_temp)
    if (exp(-B*del_e) .ge. rand_temp) then
      lattice(rand_loc_int(1),rand_loc_int(2))=-spin
      energy_sum=energy_sum+del_e
    end if 
end subroutine !}}}
!------------------------------------------------------------------------------
real(8) function energy(lattice1,latticeSize1) !{{{
 
  integer :: latticeSize1
  integer,dimension(latticeSize1,latticeSize1) :: lattice1
  
  integer :: ii,jj
  integer :: spin,h1,h2,h3,h4
  
  energy=0
  
  do ii = 2,latticeSize1 - 1
      do jj = 2,latticeSize1 - 1
          spin = lattice1(ii,jj)
          h1 = lattice1(ii + 1,jj + 1)
          h2 = lattice1(ii + 1,jj - 1)
          h3 = lattice1(ii - 1,jj + 1)
          h4 = lattice1(ii - 1,jj - 1)
          energy = energy+spin*(h1+h2*h3+h4)
      end do
  end do

end function !}}}
!------------------------------------------------------------------------------
