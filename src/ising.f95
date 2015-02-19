program ising
  
  use global

  implicit none
  
  !External Functions
  real(8),external :: energy
  integer,external :: hfield

  !Variables for main routine
  integer :: ii,jj
  integer :: tot_energy
  real(8) :: energy_avg
  real(8) :: hfeild_avg
  real(8) :: temp

  !output files
  open (unit = 1, file = 'energy.out',     status = 'unknown')
  open (unit = 2, file = 'hfield.out',     status = 'unknown')
  open (unit = 3, file = 'ising.inp',      status = 'unknown')
  open (unit = 4, file = 'hfield_avg.out', status = 'unknown')
  open (unit = 5, file = 'energy_avg.out', status = 'unknown')

  !initializing variables
  ii = 0                 !Iteration Variable
  jj = 0                 !Averages Variables
  energy_sum = 0         !Energy Sum
  rand_loc = 2           !Init Random Location Variable

  call ReadInput
  
  allocate(lattice(0:latticeSize + 1,0:latticeSize + 1))
  lattice = 0
  lattice(1:latticeSize, 1:latticeSize) = 1
  
  !main program 
  
  do jj = 0,temp_points
      print*, 'Current Temp Point:',jj
      hfeild_avg=0
      energy_avg=0
      tot_energy=0
      temp=min_temp+((max_temp-min_temp)/temp_points)*jj
      B = 1/(k*temp)         !constant of temp and boltman
      tot_energy = energy(lattice,latticeSize)!Calculate Chagnge total energy
      do ii = 0,iters
          !generate a random location 
          Call random_number(rand_loc)
          rand_loc = rand_loc*(latticeSize)
          rand_loc_int = ceiling(rand_loc)
          
          call del_energy !calculate change in energy
          call flip_bit   !Perform metropolis test
          
          !Write results to file
          if (ii > iters-1) then   
            write (1,*),ii, energy_sum 
            write (2,*),ii, sum(lattice)/float(latticeSize**2)
          end if
          if (ii > iters-t_avg) then
              hfeild_avg=hfeild_avg+(1/(float(ii)+1))*(sum(lattice)/float(latticeSize**2)-hfeild_avg)
              energy_avg=energy_avg+(1/(float(ii)+1))*(tot_energy/float(latticeSize*2)-energy_avg)
          end if 
      end do
      
      !writing averages to file
      write(4,*),temp,hfeild_avg
      write(5,*),temp,energy_avg

  end do 
end program 

!------------------------------------------------------------------------------
subroutine del_energy!{{{

  use global, only: del_e,rand_loc_int,lattice

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
  
  use global, only: lattice,rand_loc_int,rand_temp,energy_sum,del_e,B
  
  integer :: spin
  
  spin = 0
  spin = lattice(rand_loc_int(1),rand_loc_int(2))
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
subroutine ReadInput !{{{
    
    use global, only: latticeSize,iters,k,min_temp,max_temp,temp_points, &
                      t_avg,k,j  

    read(3,*) latticeSize
    read(3,*) iters
    read(3,*) min_temp,max_temp,temp_points
    read(3,*) t_avg
    read(3,*) k
    read(3,*) j
  

    close(3)

end subroutine !}}}
!------------------------------------------------------------------------------
