! program to calculate the mean and std deviation for the temperature
!sampling for the ising model

      program stats

      !defome variables
      implicit none
      integer i,p
      real*8 mval(10),mtot,mavg,var,std_dev
      character (len=10) Tfile
      real*8 T,iter
      
      !read temp, # of iterations and file name from input file
      read(*,*) T
      read(*,*) iter
      read(*,*) Tfile

      !read the 9 magnetization values calcuated in ising.f
      open(30,file=Tfile,status='old',action='read')
      do i=1,9
        read(30,*) mval(i)
	!print *, mval(i)
      enddo 
!20    continue

      ! calculate the average m for the 9 trials
      mtot=0
      do p=1,9
        mtot=mtot+mval(p)
	!print *, mval(p),mtot
      enddo 
      mavg=mtot/9d0
      !print *, mtot,mavg
      
      !calcuate the standard deviation - and use as error
      var=0
      do p=1,9
        var=var+(mavg-mval(p))**2
	!print *, var
      enddo 
      std_dev=var**(0.5000000)
      print *, mavg,std_dev
      
      !print 
      open(31,file='mag_temp.txt')
      write(31,*) T,mavg,std_dev
      close(31)
      
      end program stats
