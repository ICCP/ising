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
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
real(REAL64) :: ham, rn, flip_prb, temp, dE, mag, site_norm, dT, avg_mag, avg_E
integer :: i, j, cnt, nt

write(output_unit,'(A)') '*** ICCP Project 1: Ising model ***'

call read_settings(settings)
call print_settings_summary(settings)
call init_random_seed

if (settings%log_flag) open(unit=15,file='log.out',status='unknown')
open(unit=16,file='data.out',status='unknown')
open(unit=17,file='avgs.out',status='unknown')

site_norm = dble(settings%Nx * settings%Ny)
allocate(spin(settings%Nx,settings%Ny))

dT = 0.d0
if (settings%n_temp > 1) then
  dT = (settings%temp_max - settings%temp_min) / dble(settings%n_temp - 1)
endif

temp_loop: do nt = 1,settings%n_temp
  spin = 1
  mag = sum(spin)
  ham = 0.d0
  avg_mag = 0.d0
  avg_E = 0.d0

  temp = settings%temp_min + (nt-1) * dT

  write(output_unit,*) temp,'(',nt,'/',settings%n_temp,')'

  ! compute initial energy
  call hamiltonian(settings,spin,ham)
  if (settings%log_flag) write(15,*) 'Starting E =',ham

  ! run the monte carlo loop
  monte_loop: do cnt = 1,settings%max_mc_iter

    if (settings%log_flag) write(15,*) 'iter',cnt,' T =',temp,' E =',ham/site_norm,' M =',mag/site_norm

    ! write to a data file if there is only one temperature to simulate
    if (settings%n_temp == 1) write(16,*) cnt,temp,ham/site_norm,mag/site_norm

    if (cnt > settings%max_mc_iter - settings%last_avg + 1) then
      avg_mag = avg_mag + mag/site_norm/dble(settings%last_avg)
      avg_E = avg_E + ham/site_norm/dble(settings%last_avg)
    endif

    ! uniformly distributed random numbers give lattice site
    ! to potentially flip
    call random_number(rn)
    i = floor(real(settings%Nx)*rn) + 1
    call random_number(rn)
    j = floor(real(settings%Ny)*rn) + 1

    call nn_energy(settings,spin,i,j,dE)
    if (settings%log_flag) write(15,*) 'dE =',dE

    ! keep the new state if energy has decreased, otherwise
    ! retain the new configuration with probability flip_prb
    flip_prb = exp(-dE/temp)
    call random_number(rn)
    if (rn <= flip_prb) then
      spin(i,j) = -spin(i,j)
      ham = ham + dE
      mag = mag + 2.d0 * dble(spin(i,j))
    endif

  enddo monte_loop
  write(17,*) temp,avg_mag,avg_E
enddo temp_loop

if (settings%log_flag) close(15)
close(16)
close(17)

open(unit=8,file='spin.out',status='unknown')
do i = 1,settings%Nx
  do j = 1,settings%Ny
    write(8,*) i,j,spin(i,j)
  enddo
enddo
close(8)

call exit(0)

end program ising
!---------------------------------------------------------------------
! Compute the Hamiltonian of the spin lattice (no impressed field).
! Sum of nearest-neighbor products (no diagonals) can be computed
! using inner products between pairs of neighboring rows and columns.
!
! Inputs:
! - (settings) settings:  Settings data structure
! - (integer4) spin:      Spin lattice
! Outputs:
! - (double real) ham:    Hamiltonian (energy) of spin lattice
!---------------------------------------------------------------------
subroutine hamiltonian(settings,spin,ham)
  use iso_fortran_env
  use inputs

  implicit none

  type(settings_type) :: settings
  integer, dimension(settings%Nx,settings%Ny) :: spin
  integer, dimension(settings%Nx) :: rv1, rv2
  integer, dimension(settings%Ny) :: cv1, cv2
  real(REAL64) :: ham, Jcoup
  integer :: i

  Jcoup = 1.d0

  ham = 0.d0

  do i = 1,settings%Nx-1
    rv1 = spin(i,:)
    rv2 = spin(i+1,:)
    ham = ham - real(sum(rv1 * rv2))
  enddo
  
  do i = 1,settings%Ny-1
    cv1 = spin(:,i)
    cv2 = spin(:,i+1)
    ham = ham - real(sum(cv1 * cv2))
  enddo

  ham = ham * Jcoup

  return
end subroutine hamiltonian
!---------------------------------------------------------------------
! Stream data from /dev/urandom to use as seed fodder for RNG
! Code taken from:
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!---------------------------------------------------------------------
subroutine init_random_seed
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
  form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
    read(un) seed
    close(un)
  else
  ! Fallback to XOR:ing the current time and pid. The PID is
  ! useful in case one launches multiple instances of the same
  ! program in parallel.
    call system_clock(t)
    if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
      + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
      + dt(3) * 24_int64 * 60 * 60 * 1000 &
      + dt(5) * 60 * 60 * 1000 &
      + dt(6) * 60 * 1000 + dt(7) * 1000 &
      + dt(8)
    end if
    pid = getpid()
    t = ieor(t, int(pid, kind(t)))
    do i = 1, n
      seed(i) = lcg(t)
    end do
  end if
  call random_seed(put=seed)

  contains
  
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
    s = 104729
    else
    s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
!---------------------------------------------------------------------
! Compute nearest neighbors energy change
!---------------------------------------------------------------------
subroutine nn_energy(settings,spin,row,col,dE)
  use iso_fortran_env
  use inputs

  implicit none

  type(settings_type) :: settings
  integer, dimension(settings%Nx,settings%Ny) :: spin
  real(REAL64) :: dE
  integer :: row, col

  dE = 0.d0

  if (row-1 > 0) dE = dE + spin(row-1,col) * spin(row,col)
  if (row+1 <= settings%Nx) dE = dE + spin(row+1,col) * spin(row,col)
  if (col-1 > 0) dE = dE + spin(row,col-1) * spin(row,col)
  if (col+1 <= settings%Ny) dE = dE + spin(row,col+1) * spin(row,col)

  dE = 2.d0 * dE
  
  return
end subroutine nn_energy
