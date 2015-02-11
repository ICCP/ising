program ising
  
  use global

  implicit none
  
  !External Functions
  real(8),external :: energy
  integer,external :: hfield

  integer :: ii
  integer :: tot_energy
  open (unit = 1, file = 'energy.out',status = 'unknown')
  open (unit = 2, file = 'hfield.out',status = 'unknown')

  latticeSize=10  !Lattice Size
  ii=0            !Iteration Variable
  iters=1000      !Number of iterations the program runs through
  rand_loc=2      !Init Random Location Variable
  n1=1            !Random Number Generator Seed lattice
  n2=2            !Random Number Generator Seed temp
  j=-1.d0         !magnetic coeffient
  k=1             !Nomarlized Boltzman Constant
  temp= .5          !Tempereature of simulation
  B=1/(k*temp)    !constant of temp and boltman
  energy_sum=0    !Energy Sum

  allocate(lattice(latticeSize,latticeSize))
  lattice=1
  lattice(:,1)=0
  lattice(:,latticeSize)=0
  lattice(1,:)=0
  lattice(latticeSize,:)=0
   
  tot_energy=energy(lattice,latticeSize)
  do ii=0,iters
      Call random_seed(n1)
      Call random_number(rand_loc)
      rand_loc=rand_loc*(latticeSize-2)
      rand_loc_int=int(rand_loc)+2
      call del_energy
      call flip_bit
      write (1,*),ii,energy_sum 
      write (2,*),ii,hfield(lattice,latticeSize)

  end do
end program 

subroutine del_energy

  use global

  integer :: spin,h1,h2,h3,h4
      spin=lattice(rand_loc_int(1),rand_loc_int(2))
      h1=lattice(rand_loc_int(1)+1,rand_loc_int(2)+1)
      h2=lattice(rand_loc_int(1)+1,rand_loc_int(2)-1)
      h3=lattice(rand_loc_int(1)-1,rand_loc_int(2)+1)
      h4=lattice(rand_loc_int(1)-1,rand_loc_int(2)-1)
      h(1)=j*spin*(h1+h2+h3+h4)
      print*,'h(1):',h(1)
      if(spin==-1) then
          spin=1
          h(2)=j*spin*(h1+h2+h3+h4)
      else
          spin=-1
          h(2)=j*spin*(h1+h2+h3+h4)
      end if
      print*,'h(2)',h(2)
      del_e=h(1)-h(2)
      print*,'del_e',del_e
end subroutine

subroutine flip_bit

  use global
    
  integer :: spin
  spin = 0
    if (del_e .lt. 0) then  
        spin=lattice(rand_loc_int(1),rand_loc_int(2))
        lattice(rand_loc_int(1),rand_loc_int(2))=-spin
        energy_sum=energy_sum+del_e
    else
        temp_prob=exp(-B*del_e)
        call random_seed(n2)
        call random_number(rand_temp)
        if (rand_temp .lt. temp_prob) then
            lattice(rand_loc_int(1),rand_loc_int(2))=-spin
            energy_sum=energy_sum+del_e
        end if 
    end if 
end subroutine 

real(8) function energy(lattice1,latticeSize1)
 
  integer :: ii,jj
  integer :: spin,h1,h2,h3,h4
  integer :: latticeSize1
  integer,dimension(latticeSize1,latticeSize1) :: lattice1

  do ii=2,latticeSize1-1
      do jj=2,latticeSize1-1
          spin=lattice1(ii,jj)
          h1=lattice1(ii+1,jj+1)
          h2=lattice1(ii+1,jj-1)
          h3=lattice1(ii-1,jj+1)
          h4=lattice1(ii-1,jj-1)
          energy=energy+spin*(h1+h2*h3+h4)
      end do
  end do

  print*,'energy_sum',energy
end function

integer function hfield(lattice1,latticeSize1)
  
  integer :: ii,jj !local indexers
  integer :: latticeSize1
  integer,dimension(latticeSize1,latticeSize1) :: lattice1
  hfield=0
  do ii=1,latticeSize1
    do jj=1,latticeSize1
       hfield=hfield+lattice1(ii,jj)
    end do
  end do

end function hfield
