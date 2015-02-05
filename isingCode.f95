program isingtest

IMPLICIT NONE

integer :: i, iteration, delta_energy, N, T
integer :: j, state, energy_new, energy_old
integer :: i_swap_value, j_swap_value
integer :: wholeRand, k
integer, dimension(10,10, 50) :: a
integer, parameter :: out_unit=20 
real, dimension(10,2) :: m

real :: rand, beta, g, kB

open (unit=out_unit,file="results.txt",action="write",status="replace")

! Define Boltzmann Constant

kB=1.3806488*(10**(-23))

! Total number of spots in square lattice

N=10*10


! The first random value is always the same, so it is
! created and not used
 
rand = get_random_number()
wholeRand = get_random_whole_number(rand)

! Loop through different temperatures
k=1
do T= 100,1000,100

    ! Define beta value for new temp

    beta = 1/(kB*T)

    ! Reset all spins in array a to positive 1

    do state = 1,50
        do i=1,10
           do j=1,10
              a(i,j, state)=1
           end do
        end do
    end do


     ! Loop to minimize 50 states
     do state = 1,50

          ! Iterate through random flips 100,000 times to find a minimized energy state
          do iteration=1,100000

             ! Picking random i and j coors to determine which spin to flip

             rand = get_random_number()
             i_swap_value = get_random_whole_number(rand)

             rand = get_random_number()
             j_swap_value = get_random_whole_number(rand)

             a(i_swap_value, j_swap_value, state) = a(i_swap_value, j_swap_value,state)*(-1)

             rand = get_random_number()
         
             if (i_swap_value>1) then
                  energy_new =energy_new+(a(i_swap_value, j_swap_value,state)* a((i_swap_value-1), j_swap_value, state))
                  energy_old = energy_old+ (a(i_swap_value, j_swap_value,state)* a((i_swap_value-1), j_swap_value,state)*(-1))
             end if
             if (i_swap_value<10) then
                  energy_new = energy_new+ (a(i_swap_value, j_swap_value,state)* a((i_swap_value+1), j_swap_value, state))
                  energy_old = energy_old+ (a(i_swap_value, j_swap_value,state)* a((i_swap_value+1), j_swap_value,state)*(-1))
             end if
             if (j_swap_value>1) then
                  energy_new =  energy_new+(a(i_swap_value, j_swap_value,state)* a(i_swap_value, (j_swap_value-1),state))
                  energy_old = energy_old+ (a(i_swap_value, j_swap_value,state)* a(i_swap_value, (j_swap_value-1),state)*(-1))
             end if
             if (j_swap_value<10) then
                  energy_new = energy_new+ (a(i_swap_value, j_swap_value,state)* a(i_swap_value, (j_swap_value+1),state))
                  energy_old = energy_old+ (a(i_swap_value, j_swap_value,state)* a(i_swap_value, (j_swap_value+1),state)*(-1))
             end if
         
             ! Multiplying the sum by -J
             energy_new = -1*energy_new
             energy_old = -1*energy_old

             ! Calculate the change in energy
             delta_energy = energy_new- energy_old
         
             if (delta_energy>0) then
                  rand = get_random_number()
                  g = exp((-beta*delta_energy))
                  if (g<=rand) then
                       a(i_swap_value, j_swap_value,state)= a(i_swap_value, j_swap_value,state)*(-1)
                  end if
             end if      
          
             ! Resetting energy values 
             energy_new=0
             energy_old=0      
        end do
         
     end do
     
     !Record m for given T value

     m(k,1) = T
     m(k,2) = 0

     do i=1,10
         do j=1,10
                m(k,2)=m(k,2)+a(i,j,state)
         end do
     end do

     m(k,2)=m(k,2)/N

     write (out_unit,*) m(k,1), m(k,2)
     
    k=k+1
end do

close(out_unit)



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

    function get_random_whole_number(y)
       integer :: get_random_whole_number
       real :: y
       get_random_whole_number = (10*y)+1
       if (get_random_whole_number==11) then
            get_random_whole_number = (10*y)+1
       end if
    end function

end program isingtest
