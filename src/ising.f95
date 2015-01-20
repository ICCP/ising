program ising
  implicit none

  integer :: nx, ny
  integer, allocatable :: sigma(:,:)
  integer :: i, j

  nx = 10
  ny = 10

  allocate (sigma(nx,ny))

  ! Should be lowest energy
  do j=1,ny
    do i=1,nx
      sigma(i,j) = 1
      write(*,"(i3)",advance='no') sigma(i,j)
    end do
    write(*,*)
  end do
  write (*,*) hamiltonian (nx, ny, sigma)

  ! Should be higher energy
  do j=1,ny
    do i=1,nx
      sigma(i,j) = -1
      write(*,"(i3)",advance='no') sigma(i,j)
    end do
    write(*,*)
  end do
  write (*,*) hamiltonian (nx, ny, sigma)

  ! Should be highest energy
  do j=1,ny
    do i=1,nx
      sigma(i,j) = (-1)**(i+j)
      write(*,"(i3)",advance='no') sigma(i,j)
    end do
    write(*,*)
  end do
  write (*,*) hamiltonian (nx, ny, sigma)

  deallocate (sigma)

contains
  function hamiltonian (nx, ny, sigma) result (energy)
    integer, intent(in) :: nx, ny
    integer, intent(in) :: sigma(nx,ny)
    real, parameter     :: inter = 1, field = 1
    integer             :: i, j
    real                :: energy

    energy = 0

    do j=1,ny
      do i=1,nx
        ! Interaction with nearest neighbors
        if (i /= 1) then
          energy = energy - inter * sigma(i,j) * sigma(i-1,j)
        end if
        if (j /= 1) then
          energy = energy - inter * sigma(i,j) * sigma(i,j-1)
        end if

        ! Interaction with field
        energy = energy - field * sigma(i,j)
      end do
    end do
  end function
end program
