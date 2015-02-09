program ising
  
  use global

  implicit none
  
  latticeSize=10  !Lattice Size
  ii=0            !Iteration Variable
  iters=3         !Number of iterations the program runs through
  rand_loc=2      !Init Random Location Variable
  n=1             !Random Number Generator Seed
  j=-1.d0         !magnetic coeffient
  spin=0          !spin initation

  allocate(lattice(latticeSize,latticeSize))
  lattice=1

  do ii=0,iters
      Call random_seed(n)
      Call random_number(rand_loc)
      rand_loc=rand_loc*latticeSize
      rand_loc_int=int(rand_loc)+1
      print*, rand_loc_int
  end do

  call del_energy(lattice)
  print*, del_e
end program 

subroutine del_energy

  use global
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
