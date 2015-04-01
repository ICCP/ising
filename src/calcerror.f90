   ! calculate the mean and error for each temperature/chain length set
   
   program error
   implicit none
   integer i,nruns,io
   real*8 logn,total,mean,variance,err
   real*8, allocatable :: logrg(:)
   
   ! read in file until eof to get number of data points
   nruns = 0
   do
      read(5,*,iostat=io) 
      if (io/=0) exit
      nruns = nruns + 1
   enddo 
   
   rewind(5)
   
   allocate(logrg(nruns))
   
   ! read file to get logR and logN
   total=0
   do i=1,nruns
      read(5,*) logrg(i),logn
      total = total + logrg(i)
   enddo  
   
   ! calculate mean
   mean = total/real(nruns)
   
   ! calculate variance
   total=0
   do i=1,nruns
      variance = (mean - logrg(i))**2
      total = total + variance
   enddo
   
   ! calculate the error
   err = sqrt(total/real(nruns))
   
   print *, mean, err
   
   end program error
