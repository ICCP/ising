program ising
  use rand_tools
  implicit none

  integer,allocatable :: sigma(:,:)
  integer             :: nx,ny
  integer             :: nsteps
  real(8)             :: inter
  real(8)             :: field
  real(8)             :: beta
  integer             :: iostatus
  integer             :: uin,uout
  character(64)       :: a,b
  character(64)       :: finname,foutname
  integer             :: initmode
  integer             :: numout

  call init_random_seed()
  numout = 0

  uin  = 50
  uout = 51

  ! Default options
  nx       = 10
  ny       = 10
  nsteps   = 0
  inter    = 1
  field    = 0
  beta     = 1E0
  foutname = 'out'
  initmode = 1

  if(iargc()>0) then
    call getarg(1,finname)
    open(unit=uin,file=finname,action='read')
    ! Read in options
    ! Format is:
    ! <variable to change> <whitespace> <value to give variable>
    ! Input that doesn't make sense is ignored
    ! TODO: Work on handling weird input files
    do
      ! Read in two strings
      read(uin,*,IOSTAT=iostatus) a,b
      ! EOF ends input
      if(iostatus<0) then
        exit
      ! Parse input
      else
        if(a == 'nx')       read(b,*,IOSTAT=iostatus) nx
        if(a == 'ny')       read(b,*,IOSTAT=iostatus) ny
        if(a == 'nsteps')   read(b,*,IOSTAT=iostatus) nsteps
        if(a == 'inter')    read(b,*,IOSTAT=iostatus) inter
        if(a == 'field')    read(b,*,IOSTAT=iostatus) field
        if(a == 'beta')     read(b,*,IOSTAT=iostatus) beta
        if(a == 'initmode') read(b,*,IOSTAT=iostatus) initmode
        if(a == 'foutname') read(b,*,IOSTAT=iostatus) foutname
        if(iostatus/=0)     cycle
      end if
    end do
    close(unit=uin)
  end if

  ! Output input
  write(*,*) '# Simulation Parameters'
  write(*,*) '# nx=', nx
  write(*,*) '# ny=', ny
  write(*,*) '# nsteps=', nsteps
  write(*,*) '# inter=', inter
  write(*,*) '# field=', field
  write(*,*) '# beta=', beta
  write(*,*) '# initmode=', initmode
  write(*,*) '# foutname=', foutname

  allocate(sigma(nx,ny))
  open(unit=uout,file=foutname,status='replace',form='unformatted',access='stream',action='write')
  write(unit=uout) nx,ny,numout

  call initialize(initmode)
  write(unit=uout) sigma

  call montecarlo()

  rewind(unit=uout)
  write(unit=uout) nx,ny,numout

  write(*,*) nx,ny,numout

  close(unit=uout)
  deallocate (sigma)

contains
  ! Initialize sigma
  ! initmode - specifies the way in which we will initialize the array
  !   1: every spin is randomly assigned
  !   2: every spin is up
  !   3: every spin is down
  !   default: 1
  subroutine initialize(initmode)
    integer, intent(in) :: initmode
    real, allocatable   :: tmparr(:,:)

    if(initmode==2) then ! up
      sigma = 1
    else if(initmode==3) then ! down
      sigma = -1
    else ! random (default)
      allocate(tmparr(nx,ny))
      call random_number(tmparr)
      sigma = floor(tmparr*2)*2-1
      deallocate(tmparr)
    end if
  end subroutine initialize

  ! Flip a random spin
  ! If that was energetically favorable, keep it
  ! If not, keep it depending on temperature and energy change
  ! Repeat for nsteps
  subroutine montecarlo()
    integer             :: i,j,k
    integer,allocatable :: newsigma(:,:)
    real                :: tmp
    real(8)             :: Ediff


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
      if(tmp < exp(-beta*Ediff)) then
        sigma = newsigma
        numout=numout+1
        write(unit=uout) k,i,j
      end if
    end do
    deallocate(newsigma)

  end subroutine montecarlo

  function neighb_contrib(i,j,s) result(contrib)
    integer,intent(in) :: i,j
    integer,intent(in) :: s(:,:)
    real(8)            :: contrib

    contrib = 0

    ! Interaction between nearest neighbors
    if(i /= 1)  contrib = contrib - inter*s(i,j)*s(i-1,j)
    if(j /= 1)  contrib = contrib - inter*s(i,j)*s(i,j-1)
    if(i /= nx) contrib = contrib - inter*s(i,j)*s(i+1,j)
    if(j /= ny) contrib = contrib - inter*s(i,j)*s(i,j+1)

  end function neighb_contrib

  function field_contrib(i,j,s) result(contrib)
    integer,intent(in) :: i,j
    integer,intent(in) :: s(:,:)
    real(8)            :: contrib

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
end program ising
