module global

  implicit none

  integer :: t_avg,temp_points
  integer :: latticeSize,iters,n1,n2
  integer, allocatable :: lattice(:,:)
  real(8) :: rand_loc(2)
  real(8) :: rand_temp
  real(8) :: temp_prob
  integer :: rand_loc_int(2)
  real(8) :: del_e
  real(8) :: h(2)
  real(8) :: j                  !magnetic coefficients
  real(8) :: min_temp,max_temp
  real(8) :: k                  !boltzman constant
  real(8) :: B                  !combination of temp and boltzman

  !Energy Variables
  real(8) :: energy_sum         !sum of total energy in lattice

end module global 
