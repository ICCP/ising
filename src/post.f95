program post
  use plot

  integer             :: i
  integer             :: nx,ny,nsteps
  integer,allocatable :: sigma(:,:)
  integer             :: uopts,uin,umag
  character(64)       :: a,b
  character(64)       :: foptsname,finname,magoutname
  integer             :: iostatus
  real(8),allocatable :: m(:)
  integer             :: step,x,y

  uopts = 50
  uin   = 51
  umag  = 52

  ! Default options
  finname = 'out'
  magoutname = 'mag.out'

  if(iargc()>0) then
    call getarg(1,foptsname)
    open(unit=uopts,file=foptsname,action='read')
    ! Read in options
    ! Format is:
    ! <variable to change> <whitespace> <value to give variable>
    ! Input that doesn't make sense is ignored
    ! TODO: Work on handling weird input files
    do
      ! Read in two strings
      read(uopts,*,IOSTAT=iostatus) a,b
      ! EOF ends input
      if(iostatus<0) then
        exit
      ! Parse input
      else
        if(a == 'finname') read(b,*,IOSTAT=iostatus) finname
        if(a == 'magoutname') read(b,*,IOSTAT=iostatus) magoutname
        if(iostatus/=0)    cycle
      end if
    end do
    close(unit=uopts)
  end if

  open(unit=uin,file=finname,form='unformatted',access='stream',action='read')
  open(unit=umag,file=magoutname,action='write')
  read(unit=uin) nx,ny,nsteps

  allocate(sigma(nx,ny))
  allocate(m(nsteps+1))
  read(unit=uin) sigma
  call plot_init(nx,ny)

  call plbop()
  call plot_lattice(sigma)
  call pleop()

  i=1
  m(i) = magnetization()
  write(umag,*) 0,m(i)
  do
    i = i+1
    read(unit=uin,iostat=iostatus) step,x,y

    if(iostatus>0) then
      write(0,*) 'There was a problem reading in the file'
      call exit(iostatus)
    else if(iostatus<0) then
      exit
    else
      sigma(x,y) = -sigma(x,y)
      m(i) = magnetization()
      write(umag,*) step,m(i)

!      call plbop()
!      call plot_lattice(sigma)
!      call pleop()
    end if
  end do

  call plbop()
  call plot_lattice(sigma)
  call pleop()

  call plot_close()
  deallocate(m)
  deallocate(sigma)
  close(unit=umag)
  close(unit=uin)

contains
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
end program post
