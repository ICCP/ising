module global

  implicit none

  integer :: latticeSize,ii,iters,n
  integer, allocatable :: lattice(:,:)
  real(8) :: rand_loc(2)
  integer :: rand_loc_int(2)
  real(8) :: del_e
  real(8) :: h(2)
  real(8) :: j !magnetic coefficients
  real(8) :: temp
  integer :: spin
  integer :: h1,h2,h3,h4 
end module global 
