module inputs

use iso_fortran_env

type :: settings_type
  integer :: Nx, Ny, max_mc_iter, last_avg, n_temp
  logical :: log_flag
  real(REAL64) :: temp_min, temp_max
end type settings_type

contains
!---------------------------------------------------------------------
! Read settings from file 'settings.inp'
!---------------------------------------------------------------------
subroutine read_settings(settings)
  implicit none

  type(settings_type) :: settings
  integer :: stat

  open(unit=12,file='settings.inp',status='unknown',iostat=stat)
  if (stat /= 0) then
    write(error_unit,'(A)') 'Error opening settings file "settings.inp"'
    call exit(10)
  endif
  read(12,*) settings%Nx, settings%Ny
  read(12,*) settings%temp_min, settings%temp_max, settings%n_temp
  read(12,*) settings%max_mc_iter
  read(12,*) settings%last_avg
  read(12,*) settings%log_flag
  close(12)

  return
end subroutine read_settings
!---------------------------------------------------------------------
! Print a summary of the simulation settings for logging purposes
!---------------------------------------------------------------------
subroutine print_settings_summary(settings)
  implicit none
  
  type(settings_type) :: settings

  write(output_unit,'(A)') '============ Settings summary ============'
  write(output_unit,'(A,2I8)')  '-> Grid settings (Nx,Ny): ',&
                                settings%Nx,settings%Ny
  if (settings%log_flag) write(output_unit,'(A)') '-> Logging to log.out'
  write(output_unit,'(A)') '=========================================='

  return
end subroutine print_settings_summary
!---------------------------------------------------------------------
end module inputs
