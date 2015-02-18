program isingtest

IMPLICIT NONE

integer :: i, j, flipSpins, k, Tstep, everyHundred
integer :: size, i_swap_value, j_swap_value
integer :: avgCount, mCount
integer, parameter :: out_unit=20
integer, parameter :: out_unit2=10
integer, dimension(100,100) :: spins

real, dimension(3,7000000) :: m_T
real, dimension(20) :: m_bins
real :: e_new, e_old, delta_e
real :: m, rand, kB, beta, g, T
real :: avgM, m_avg, stdev, diff
real :: diffSq, diffSqSum

character(len=10) ::  out_file
character(len=7) ::  out_file_base
character(len=3) ::  out_file_no

size=100

kB=1

out_file_base="Results"

open (unit=out_unit2,file="avg_m_values",action="write",status="replace")

do Tstep=1,2

   k=1

   T=Tstep
   beta=1/(kB*T)
   write(out_file_no,"(F3.1)") T
   out_file=out_file_base//out_file_no
   open (unit=out_unit,file=out_file,action="write",status="replace")
   do i=1,size
      do j=1,size
         spins(i,j)=1
      end do
   end do
   
   m_Avg=0
   everyHundred=1
   avgCount=1
   mCount=1 
   avgM=0   
   diffSqSum=0 
   do flipSpins=1,5000000
  
   rand = get_random_number()

   rand = get_random_number()
   i_swap_value = get_random_whole_number(rand,size)

   rand = get_random_number()
   j_swap_value = get_random_whole_number(rand,size)

   spins(i_swap_value, j_swap_value)=spins(i_swap_value, j_swap_value)*(-1)

   e_new=0
   e_old=0

   if (i_swap_value>1) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((i_swap_value-1),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((i_swap_value-1),j_swap_value)*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((size),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((size),j_swap_value)*(-1))
   end if
   if (j_swap_value>1) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value-1)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value-1))*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(size)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(size))*(-1))
   end if
   if (i_swap_value<size) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((i_swap_value+1),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((i_swap_value+1),j_swap_value)*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(1,j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(1,j_swap_value)*(-1))
   end if
   if (j_swap_value<size) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value+1)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value+1))*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,1))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,1)*(-1))
   end if

   e_new=(-1)*e_new
   e_old=(-1)*e_old

   delta_e=e_new-e_old


   if (delta_e>0) then
      rand=get_random_number()
      g=exp((-beta*delta_e))
      if(g<=rand) then
         spins(i_swap_value,j_swap_value)=spins(i_swap_value,j_swap_value)*(-1)
      end if
   end if
   
   if (everyHundred == 100) then
       m=0
       do i=1,size
           do j=1,size
               m=m+spins(i,j) 
           end do
       end do

       m=m/(size*size)
	   if (m<0) then
	       m=(-1)*m
	   end if
       m_T(2,k)=m
       m_T(1,k)=flipSpins
       m_T(1,k)=m_T(1,k)/1000000.0 
       write (out_unit,*) m_T(1,k),m_T(2,k)
	   
	   if (flipSpins>2999999) then
	       m_avg=m_avg+m
		   if (avgCount==1000) then
			   m_avg=m_avg/avgCount
			   avgCount=0
			   m_bins(mCount)=m_avg
			   mCount=mCount+1
		   end if
		   avgCount=avgCount+1
	   end if
	   
   end if
	   
   k=k+1
   if (everyHundred==100) then
       everyHundred=1
   else   
       everyHundred=everyHundred+1
   end if
   
    !Iterations end
   end do
   mCount=mCount-1
   do i=1,mCount
       avgM=avgM+m_bins(i)
   end do
   avgM=avgM/mCount
   
   do i=1,mCount
       diff=avgM-m_bins(i)
	   diffSq=diff**2
	   diffSqSum=diffSqSum+diffSq
   end do
       stdev=sqrt((diffSqSum/mCount))
	  
   write (out_unit2,*) T, avgM, stdev
   close(out_unit)
   
!Temp end
end do

do Tstep=1,9

   k=1

   T=2+(Tstep/10.0)
   beta=1/(kB*T)
   write(out_file_no,"(F3.1)") T
   out_file=out_file_base//out_file_no
   open (unit=out_unit,file=out_file,action="write",status="replace")
   do i=1,size
      do j=1,size
         spins(i,j)=1
      end do
   end do
   
   m_Avg=0
   everyHundred=1
   avgCount=1
   mCount=1 
   avgM=0   
   diffSqSum=0 
   do flipSpins=1,5000000
  
   rand = get_random_number()

   rand = get_random_number()
   i_swap_value = get_random_whole_number(rand,size)

   rand = get_random_number()
   j_swap_value = get_random_whole_number(rand,size)

   spins(i_swap_value, j_swap_value)=spins(i_swap_value, j_swap_value)*(-1)

   e_new=0
   e_old=0

   if (i_swap_value>1) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((i_swap_value-1),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((i_swap_value-1),j_swap_value)*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((size),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((size),j_swap_value)*(-1))
   end if
   if (j_swap_value>1) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value-1)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value-1))*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(size)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(size))*(-1))
   end if
   if (i_swap_value<size) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((i_swap_value+1),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((i_swap_value+1),j_swap_value)*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(1,j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(1,j_swap_value)*(-1))
   end if
   if (j_swap_value<size) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value+1)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value+1))*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,1))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,1)*(-1))
   end if

   e_new=(-1)*e_new
   e_old=(-1)*e_old

   delta_e=e_new-e_old


   if (delta_e>0) then
      rand=get_random_number()
      g=exp((-beta*delta_e))
      if(g<=rand) then
         spins(i_swap_value,j_swap_value)=spins(i_swap_value,j_swap_value)*(-1)
      end if
   end if
   
   if (everyHundred == 100) then
       m=0
       do i=1,size
           do j=1,size
               m=m+spins(i,j) 
           end do
       end do

       m=m/(size*size)
	   if (m<0) then
	       m=(-1)*m
	   end if
       m_T(2,k)=m
       m_T(1,k)=flipSpins
       m_T(1,k)= m_T(1,k)/1000000
       write (out_unit,*) m_T(1,k),m_T(2,k)
	   
	   if (flipSpins>2999999) then
	       m_avg=m_avg+m
		   if (avgCount==1000) then
			   m_avg=m_avg/avgCount
			   avgCount=0
			   m_bins(mCount)=m_avg
			   mCount=mCount+1
		   end if
		   avgCount=avgCount+1
	   end if
	   
   end if
	   
   k=k+1
   if (everyHundred==100) then
       everyHundred=1
   else   
       everyHundred=everyHundred+1
   end if
   
    !Iterations end
   end do
   mCount=mCount-1
   do i=1,mCount
       avgM=avgM+m_bins(i)
   end do
   avgM=avgM/mCount
   
   do i=1,mCount
       diff=avgM-m_bins(i)
	   diffSq=diff**2
	   diffSqSum=diffSqSum+diffSq
   end do
       stdev=sqrt((diffSqSum/mCount))
	  
   write (out_unit2,*) T, avgM, stdev
   close(out_unit)
   
!Temp end
end do

do Tstep=3,4

   k=1

   T=Tstep
   beta=1/(kB*T)
   write(out_file_no,"(F3.1)") T
   out_file=out_file_base//out_file_no
   open (unit=out_unit,file=out_file,action="write",status="replace")
   do i=1,size
      do j=1,size
         spins(i,j)=1
      end do
   end do
   
   m_Avg=0
   everyHundred=1
   avgCount=1
   mCount=1 
   avgM=0   
   diffSqSum=0 
   do flipSpins=1,5000000
  
   rand = get_random_number()

   rand = get_random_number()
   i_swap_value = get_random_whole_number(rand,size)

   rand = get_random_number()
   j_swap_value = get_random_whole_number(rand,size)

   spins(i_swap_value, j_swap_value)=spins(i_swap_value, j_swap_value)*(-1)

   e_new=0
   e_old=0

   if (i_swap_value>1) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((i_swap_value-1),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((i_swap_value-1),j_swap_value)*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((size),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((size),j_swap_value)*(-1))
   end if
   if (j_swap_value>1) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value-1)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value-1))*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(size)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(size))*(-1))
   end if
   if (i_swap_value<size) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins((i_swap_value+1),j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins((i_swap_value+1),j_swap_value)*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(1,j_swap_value))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(1,j_swap_value)*(-1))
   end if
   if (j_swap_value<size) then
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value+1)))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,(j_swap_value+1))*(-1))
   else
      e_new=e_new+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,1))
      e_old=e_old+(spins(i_swap_value,j_swap_value)*spins(i_swap_value,1)*(-1))
   end if

   e_new=(-1)*e_new
   e_old=(-1)*e_old

   delta_e=e_new-e_old


   if (delta_e>0) then
      rand=get_random_number()
      g=exp((-beta*delta_e))
      if(g<=rand) then
         spins(i_swap_value,j_swap_value)=spins(i_swap_value,j_swap_value)*(-1)
      end if
   end if
   
   if (everyHundred == 100) then
       m=0
       do i=1,size
           do j=1,size
               m=m+spins(i,j) 
           end do
       end do

       m=m/(size*size)
	   if (m<0) then
	       m=(-1)*m
	   end if
       m_T(2,k)=m
       m_T(1,k)=flipSpins
       m_T(1,k)= m_T(1,k)/1000000
       write (out_unit,*) m_T(1,k),m_T(2,k)
	   
	   if (flipSpins>2999999) then
	       m_avg=m_avg+m
		   if (avgCount==1000) then
			   m_avg=m_avg/avgCount
			   avgCount=0
			   m_bins(mCount)=m_avg
			   mCount=mCount+1
		   end if
		   avgCount=avgCount+1
	   end if
	   
   end if
	   
   k=k+1
   if (everyHundred==100) then
       everyHundred=1
   else   
       everyHundred=everyHundred+1
   end if
   
    !Iterations end
   end do
   mCount=mCount-1
   do i=1,mCount
       avgM=avgM+m_bins(i)
   end do
   avgM=avgM/mCount
   
   do i=1,mCount
       diff=avgM-m_bins(i)
	   diffSq=diff**2
	   diffSqSum=diffSqSum+diffSq
   end do
       stdev=sqrt((diffSqSum/mCount))
	  
   write (out_unit2,*) T, avgM, stdev
   close(out_unit)
   
!Temp end
end do

close(out_unit2)


contains

   function get_random_number()
       integer :: i_seed
       integer, dimension(:), ALLOCATABLE :: Set_seed
       integer, dimension(1:8) :: dateSet_seed
       real :: get_random_number, r

       CALL RANDOM_SEED(size=i_seed)
       ALLOCATE(Set_seed(1:i_seed))
       CALL RANDOM_SEED(get=Set_seed)
       CALL DATE_AND_TIME(values=dateSet_seed)
       Set_seed(i_seed)=dateSet_seed(8)
       !Set_seed(2)=dateSet_seed(8)
       Set_seed(1)=dateSet_seed(8)*dateSet_seed(6)
       CALL RANDOM_SEED(put=Set_seed)
       DEALLOCATE(Set_seed)

       CALL RANDOM_NUMBER(r)
       get_random_number = r
    end function get_random_number

    function get_random_whole_number(y,x)
       integer :: get_random_whole_number, x
       real :: y
       get_random_whole_number = (x*y)+1
       if (get_random_whole_number==(x+1)) then
            get_random_whole_number = (x*y)+1
       end if
    end function


end program isingtest


!write (out_unit,*) spins(1,1),spins(1,2),spins(1,3),spins(1,4),spins(1,5),spins(1,6),spins(1,7),spins(1,8),spins(1,9),spins(1,10)
!write (out_unit,*) spins(2,1),spins(2,2),spins(2,3),spins(2,4),spins(2,5),spins(2,6),spins(2,7),spins(2,8),spins(2,9),spins(2,10)
!write (out_unit,*) spins(3,1),spins(3,2),spins(3,3),spins(3,4),spins(3,5),spins(3,6),spins(3,7),spins(3,8),spins(3,9),spins(3,10)
!write (out_unit,*) spins(4,1),spins(4,2),spins(4,3),spins(4,4),spins(4,5),spins(4,6),spins(4,7),spins(4,8),spins(4,9),spins(4,10)
!write (out_unit,*) spins(5,1),spins(5,2),spins(5,3),spins(5,4),spins(5,5),spins(5,6),spins(5,7),spins(5,8),spins(5,9),spins(5,10)
!write (out_unit,*) spins(6,1),spins(6,2),spins(6,3),spins(6,4),spins(6,5),spins(6,6),spins(6,7),spins(6,8),spins(6,9),spins(6,10)
!write (out_unit,*) spins(7,1),spins(7,2),spins(7,3),spins(7,4),spins(7,5),spins(7,6),spins(7,7),spins(7,8),spins(7,9),spins(7,10)
!write (out_unit,*) spins(8,1),spins(8,2),spins(8,3),spins(8,4),spins(8,5),spins(8,6),spins(8,7),spins(8,8),spins(8,9),spins(8,10)
!write (out_unit,*) spins(9,1),spins(9,2),spins(9,3),spins(9,4),spins(9,5),spins(9,6),spins(9,7),spins(9,8),spins(9,9),spins(9,10)
!write(out_unit,*)spins(10,1),spins(10,2),spins(10,3),spins(10,4),spins(10,5),spins(10,6),spins(10,7),spins(10,8),spins(10,9)
!write (out_unit,*) ' '
