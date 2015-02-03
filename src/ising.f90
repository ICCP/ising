program ising

use iso_fortran_env
use inputs

implicit none

!------------------------------------------------------------------------------
! Input parameters
!------------------------------------------------------------------------------

type(settings_type) :: settings

!------------------------------------------------------------------------------
! Lattice data
!------------------------------------------------------------------------------
integer, allocatable :: spin(:,:)

write(output_unit,'(A)') '==============================='
write(output_unit,'(A)') '* ICCP Project 1: Ising model *'
write(output_unit,'(A)') '==============================='

call read_settings(settings)
call print_settings_summary(settings)

allocate(spin(settings%Nx,settings%Ny))
spin = 1

call exit(0)

end program ising
