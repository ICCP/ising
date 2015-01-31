!program to impliment the Markov method to calculate critical 
!temperature for the 2D ising model
!intially with a finite boundary (not periodic)
!initially a ferro magnet Eij=-JSiSj with J>0

      
      function spinflip(L) result(flip)
      !chooses a random spin to flip in the LxL grid
      !one random number for the row and one for the column
      integer*8, intent(in) :: L ! input
      !real*8, intent(in)     :: iter !input
      integer*8 :: flip,p !output
      real*8 :: ran
      call random_number(ran)
       do p=1,L
         if (ran <= p*0.01) then
           flip = p
           exit
         end if
       enddo
      end function spinflip
 
      
      program isingmodel
      
      !define variables and function
      implicit none
      !variables from input file
      real*8 T
      integer ttot,fac
      character (len=10) Tfile   
      !variables defined in the program   
      integer i,j,sold,L,p,c,d,time
      real*8 num,Eold,Enew,stot,mtot
      integer, dimension(:,:),allocatable::spins
      real*8, dimension(:),allocatable::m
      integer spinflip
      
      !read in the parameters from the file given in the command line
      read(*,*) T
      read(*,*) ttot
      read(*,*) Tfile
      
      L=100 !lattice size
      !spin flip initializations (c,d) is the spin to flip
      c=1
      d=1
            
      !make an LXL array of spins, all initilized to 1
      allocate (spins(L,L))
      spins = 1
      !count how many iterations, beginning with 1
      time=1
      !total magnetization
      mtot=0
      !array of m values
      allocate (m(ttot))
      !print *, T, ttot, Tfile
      
      !seed the random number generator
      call random_seed()
          
      !repeat the Monte Carlo ttot times (read in from input file)	    
      do time=1,ttot
      Eold=0 !start with all energies equal to zero
      Enew=0
      !calculate initial energy of the system
      do i=1,L
        do j=1,L
          if ((i==L .and. j==L)) then
            Eold=Eold
          else if (i==L) then
            Eold = Eold - spins(i,j)*spins(i,j+1)
          else if (j==L) then
            Eold = Eold - spins(i,j)*spins(i+1,j)
          else
            Eold = Eold - spins(i,j)*spins(i+1,j) 
     & - spins(i,j)*spins(i,j+1)
          end if
        enddo
      enddo 
      !print *, Eold
      
      !choose which spin to flip
      c=spinflip(L)
      d=spinflip(L)
      !print *, c,d
      !flip the spin (1 -> -1 or -1 -> 1)
      if (spins(c,d)==1) then
        spins(c,d)=-1  !(1 -> -1)
        sold=1 !keep track of the old spin in case spin is flipped back
      else 
        spins(c,d)=1  !(-1 -> 1)
        sold=-1 !old spin
      end if
            
      !calculate new energy
      do i=1,L
        do j=1,L
          if ((i==L .and. j==L)) then
            Enew=Enew
           !print *, i,j,Enew
          else if (i==L) then
            Enew = Enew - spins(i,j)*spins(i,j+1)
            !print *, i,j,Enew
          else if (j==L) then
            Enew = Enew - spins(i,j)*spins(i+1,j)
            !print *, i,j,Enew
          else
            Enew = Enew - spins(i,j)*spins(i+1,j) 
     & - spins(i,j)*spins(i,j+1)
            !print *, i,j,Enew
          end if
        enddo
      enddo 
      !print *, Enew
      
      !calculate a random number to check if spin stays flipped
      call random_number(num)
      
      !check whether to keep spin flip
      if (Enew-Eold < 0) then
        Eold=Enew ! keep spin flip
        !print *, "spin flipped"
      else 
        if ( exp(-(Enew-Eold)/T) < num) then
          Enew=0
          spins(c,d)=sold !change the new value back to old value
          !print *, "spin did not flip"
        else 
          Eold=Enew !keep new spin
          !print *, "spin flipped"
        end if
      end if 
      !print *, Eold
      
      !calculate observables
      
      !calculate the total spin of the iteration
      do i=1,L
        do j=1,L
          stot=stot+spins(i,j)
        enddo
      enddo 
      
      !calculate the average magnetization for the given iteration
      m(time)=stot/(L*L) !for one iteration so no time average
      
      !write out the average magnetization for each iteration
      open(40,file="mtot.txt")
      write(40,*) time,m(time)
      stot=0 !reset stot for each iteration (new m)
      enddo
      close(40)
      
      !sample nines times to create a mean and std deviation
        !calculated in stats.f
      open(50,file=Tfile)
      fac=0.1*ttot  !averages over the last 10th of the number of iterations
        	    !in nine blocks
      do i=1,9
        do p=i*fac+1,(i+1)*fac 
          mtot=mtot+m(p)  !total magnetization in given iteration block
        enddo
        mtot=mtot/fac  !average magnetization
        write(50,*) mtot  !print to file
        mtot=0  !reset magnetization to zero
      enddo
      close(50)
      
      end program isingmodel 
