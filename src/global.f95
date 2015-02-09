module global

  implicit none

  integer :: latticeSize,ii,iters,n1,n2,jj
  integer, allocatable :: lattice(:,:)
  real(8) :: rand_loc(2)
  real(8) :: rand_temp
  real(8) :: temp_prob
  integer :: rand_loc_int(2)
  real(8) :: del_e
  real(8) :: h(2)
  real(8) :: j                  !magnetic coefficients
  real(8) :: temp
  integer :: spin
  integer :: h1,h2,h3,h4 
  real(8) :: k                  !boltzman constant
  real(8) :: B                  !combination of temp and boltzman

  !Energy Variables
  real(8) :: energy_sum         !sum of total energy in lattice
end module global 
