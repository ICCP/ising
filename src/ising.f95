program ising
  
  use global

  implicit none
  
  latticeSize=10  !Lattice Size
  ii=0            !Iteration Variable
  iters=3         !Number of iterations the program runs through
  rand_loc=2      !Init Random Location Variable
  n1=1            !Random Number Generator Seed lattice
  n2=2            !Random Number Generator Seed temp
  j=-1.d0         !magnetic coeffient
  spin=0          !spin initation
  k=1             !Nomarlized Boltzman Constant
  temp=1          !Tempereature of simulation
  B=1/(k*temp)    !constant of temp and boltman

  allocate(lattice(latticeSize,latticeSize))
  lattice=1
  lattice(:,1)=0
  lattice(:,latticeSize)=0
  lattice(1,:)=0
  lattice(latticeSize,:)=0

  do ii=0,iters
      Call random_seed(n1)
      Call random_number(rand_loc)
      rand_loc=rand_loc*(latticeSize-2)
      rand_loc_int=int(rand_loc)+2
      call del_energy
      call flip_bit
  end do
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

subroutine flip_bit

  use global
    
    if (del_e .lt. 0) then  
        lattice(rand_loc_int(1),rand_loc_int(2))=-spin
        print*,'flip neg_temp'
    else
        temp_prob=exp(-B*del_e)
        call random_seed(n2)
        call random_number(rand_temp)
        if (rand_temp .lt. temp_prob) then
            print*,'rand_temp',rand_temp,'=<','temp_prob',temp_prob
            lattice(rand_loc_int(1),rand_loc_int(2))=-spin
            print*,'flip based on prob'
        end if 
    end if 
end subroutine 

