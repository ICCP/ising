program ising
  use rand_tools
  use plot
  implicit none

  integer, allocatable :: sigma(:,:)
  integer              :: nx,ny
  integer              :: nsteps
  real(8)              :: inter
  real(8)              :: field
  real(8)              :: beta
  integer              :: iostatus
  integer              :: u
  character(64)        :: a,b
  character(64)        :: fname

  ! Default options
  nx    = 100
  ny    = 100
  nsteps= 0
  inter = 1
  field = 0
  beta  = 1E-1

  if (iargc()>0) then
    call getarg(1,fname)
    open(unit=u,file=fname)
    ! Read in options
    ! Format is:
    ! <variable to change> <whitespace> <value to give variable>
    ! Input that doesn't make sense is ignored
    ! TODO: Work on handling weird input files
    ! TODO: Allow specifying output (e.g. Magnetization file, images, etc.)
    do
      ! Read in two strings
      read (u,*,IOSTAT=iostatus) a,b
      ! EOF ends input
      if (iostatus<0) then
        exit
      ! Parse input
      else
        if (a == 'nx') then
          read (b,*,IOSTAT=iostatus) nx
        end if
        if (a == 'ny') then
          read (b,*,IOSTAT=iostatus) ny
        end if
        if (a == 'nsteps') then
          read (b,*,IOSTAT=iostatus) nsteps
        end if
        if (a == 'inter') then
          read (b,*,IOSTAT=iostatus) inter
        end if
        if (a == 'field') then
          read (b,*,IOSTAT=iostatus) field
        end if
        if (a == 'beta') then
          read (b,*,IOSTAT=iostatus) beta
        end if
        if (iostatus/=0) then
          cycle
        end if
      end if
    end do
  end if

  ! Output input
  write (*,*) "# Simulation Parameters"
  write (*,*) "# nx=", nx
  write (*,*) "# ny=", ny
  write (*,*) "# nsteps=", nsteps
  write (*,*) "# inter=", inter
  write (*,*) "# field=", field
  write (*,*) "# beta=", beta

  allocate (sigma(nx,ny))
  ! TODO: initialize sigma separately

  call plot_init(nx,ny)
  call montecarlo()
  call plot_close()

  deallocate (sigma)

contains
  subroutine montecarlo()
    integer                :: i,j,k
    integer, allocatable   :: newsigma(:,:)
    real, allocatable      :: tmparr(:,:)
    real                   :: tmp
    real(8)                :: Ediff

    ! Initialize sigma randomly
    allocate(tmparr(nx,ny))
    call init_random_seed()
    call random_number(tmparr)
    sigma = floor(tmparr*2)*2-1
    deallocate(tmparr)
    ! Initialize sigma as all up
!    do j=1,ny
!      do i=1,nx
!        sigma(i,j)=1
!      end do
!    end do

    ! Print out initial conditions
    call plot_lattice(sigma)
    call pleop()

    ! Repeat until a set number of states have been tried (or criteria met)
    allocate(newsigma(nx,ny))
    do k=1,nsteps
      ! Determine new state by flipping a random spin
      call random_number(tmp)
      i = ceiling(tmp*nx)
      call random_number(tmp)
      j = ceiling(tmp*ny)

      ! Determine whether to accept or reject new state
      ! If the energy is less, stay
      ! Else, stay with probability of exp(-beta*(H_1-H_0))
      newsigma = sigma
      newsigma(i,j) = -sigma(i,j)
      Ediff =   neighb_contrib(i,j,newsigma) + field_contrib(i,j,newsigma) &
              - neighb_contrib(i,j,sigma)    - field_contrib(i,j,sigma)
      call random_number(tmp)
      if (tmp < exp(-beta*Ediff)) then
        sigma = newsigma
      end if

      write(*,*) magnetization()
    end do
    deallocate(newsigma)

    ! Create image
    call plbop()
    call plot_lattice(sigma)
    call pleop()
  end subroutine montecarlo

  function neighb_contrib(i,j,s) result(contrib)
    integer, intent(in) :: i,j
    integer, intent(in) :: s(:,:)
    real(8)             :: contrib

    contrib = 0

    ! Interaction between nearest neighbors
    if (i /= 1) then
      contrib = contrib - inter*s(i,j)*s(i-1,j)
    end if
    if (j /= 1) then
      contrib = contrib - inter*s(i,j)*s(i,j-1)
    end if
    if (i /= nx) then
      contrib = contrib - inter*s(i,j)*s(i+1,j)
    end if
    if (j /= ny) then
      contrib = contrib - inter*s(i,j)*s(i,j+1)
    end if

  end function neighb_contrib

  function field_contrib(i,j,s) result(contrib)
    integer, intent(in) :: i,j
    integer, intent(in) :: s(:,:)
    real(8)             :: contrib

    contrib = 0

    ! Interaction with external field
    contrib = contrib - field*s(i,j)
  end function field_contrib

  function hamiltonian() result(energy)
    integer :: i,j
    real(8) :: energy

    energy = 0

    do j=1,ny
      do i=1,nx
        ! neighb_contrib is halved to prevent double counting
        energy = energy + 0.5*neighb_contrib(i,j,sigma)
        energy = energy + field_contrib(i,j,sigma)
      end do
    end do
  end function hamiltonian

  function magnetization() result(mag)
    integer :: i,j
    integer :: m
    real(8) :: mag

    m = 0
    do j=1,ny
      do i=1,nx
        m=m+sigma(i,j)
      end do
    end do
    mag = dble(m)/(nx*ny)
  end function magnetization
end program ising
