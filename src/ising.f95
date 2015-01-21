program ising
!  use rand_tools
  use plot
  implicit none

  integer              :: nx,ny
  integer, allocatable :: sigma(:,:)
  integer              :: i,j
  real, parameter      :: INTER = 1
  real, parameter      :: FIELD = 1
  real, parameter      :: BETA  = .1
  real, parameter      :: E     = exp(1.0)

  nx = 500
  ny = 500

  allocate (sigma(nx,ny))

  ! TODO: change this API
  call plot_init(nx,ny)
  call montecarlo(nx,ny,sigma)
  call plot_close()

  deallocate (sigma)

contains
  subroutine montecarlo(nx,ny,sigma)
    integer, intent(in)    :: nx,ny
    integer, intent(inout) :: sigma(nx,ny)
    integer                :: newsigma(nx,ny)
    integer                :: i,j,k
    real                   :: tmparr(nx,ny)
    real                   :: tmp
    real                   :: Ediff

    ! Initialize sigma randomly
    call init_random_seed()
    call random_number(tmparr)
    sigma = floor(tmparr*2)*2-1
    !call sigmaprint(nx,ny,sigma)
!    call plbop()
    call plot_lattice(sigma)
    call pleop()
    write(*,*)

    ! Repeat until a set number of states have been tried (or criteria met)
    do k=1,100000
      ! Determine new state by flipping a random spin
      call random_number(tmp)
      i = ceiling(tmp*nx)
      call random_number(tmp)
      j = ceiling(tmp*ny)

      ! Determine whether to accept or reject new state
      ! If the energy is less, stay
      ! Else, stay with probability of e^{-beta*(H_1-H_0)}
      newsigma = sigma
      newsigma(i,j) = -sigma(i,j)
      Ediff = neighb_contrib(i,j,nx,ny,newsigma) + field_contrib(i,j,newsigma) &
              - neighb_contrib(i,j,nx,ny,sigma) - field_contrib(i,j,sigma)
      call random_number(tmp)
      if (tmp < E**(-BETA*Ediff)) then
        sigma = newsigma
      end if
    end do
    !call sigmaprint(nx,ny,sigma)
    call plbop()
    call plot_lattice(sigma)
    call pleop()
  end subroutine

  subroutine sigmaprint(nx,ny,sigma)
    integer, intent(in) :: nx,ny
    integer, intent(in) :: sigma(nx,ny)
    integer             :: i,j

    do j=1,ny
      do i=1,nx
        if (sigma(i,j) == 1) then
          write(*,"(a)",advance='no') "+"
        end if
        if (sigma(i,j) == -1) then
          write(*,"(a)",advance='no') "-"
        end if
      end do
      write(*,*)
    end do
  end subroutine

  function neighb_contrib(i,j,nx,ny,sigma) result(contrib)
    integer, intent(in) :: i,j,nx,ny
    integer, intent(in) :: sigma(nx,ny)
    real                :: contrib

    contrib = 0

    ! Interaction between nearest neighbors
    if (i /= 1) then
      contrib = contrib - INTER*sigma(i,j)*sigma(i-1,j)
    end if
    if (j /= 1) then
      contrib = contrib - INTER*sigma(i,j)*sigma(i,j-1)
    end if
    if (i /= nx) then
      contrib = contrib - INTER*sigma(i,j)*sigma(i+1,j)
    end if
    if (j /= ny) then
      contrib = contrib - INTER*sigma(i,j)*sigma(i,j+1)
    end if

  end function

  function field_contrib(i,j,sigma) result(contrib)
    integer, intent(in) :: i,j
    integer, intent(in) :: sigma(nx,ny)
    real                :: contrib

    contrib = 0

    ! Interaction with external field
    contrib = contrib - FIELD*sigma(i,j)
  end function

  function hamiltonian(nx,ny,sigma) result(energy)
    integer, intent(in) :: nx,ny
    integer, intent(in) :: sigma(nx,ny)
    integer             :: i,j
    real                :: energy

    energy = 0

    do j=1,ny
      do i=1,nx
        ! neigb_contrib is halved to prevent double counting
        energy = energy + 0.5*neighb_contrib(i,j,nx,ny,sigma)
        energy = energy + field_contrib(i,j,sigma)
      end do
    end do
  end function

  subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = n)
    allocate(seed(n))
    open(newunit=un, file="/dev/urandom", access="stream",&
         form="unformatted", action="read", status="old", &
         iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       call system_clock(count)
       if (count /= 0) then
          t = transfer(count, t)
       else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = transfer(tms, t)
       end if
       s = ieor(t(1), t(2))
       pid = getpid() + 1099279 ! Add a prime
       s = ieor(s, pid)
       if (n >= 3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
       else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       end if
    end if
    call random_seed(put=seed)
  end subroutine init_random_seed
end program
