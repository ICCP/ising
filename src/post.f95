program post
  use plot

  integer             :: i
  integer             :: nx,ny,nsteps
  integer,allocatable :: sigma(:,:)
  integer             :: uin
  integer             :: iostatus
  real(8),allocatable :: m(:)
  integer             :: step,x,y

  uin = 50

  open(unit=uin,file='out',form='unformatted',access='stream',action='read')
  read(unit=uin) nx,ny,nsteps

  allocate(sigma(nx,ny))
  allocate(m(nsteps))
  call plot_init(nx,ny)

  read(unit=uin) sigma
  i=0
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
      write(*,*) step,m(i)

!      call plbop()
!      call plot_lattice(sigma)
!      call pleop()
    end if
  end do

  call plot_close()
  deallocate(m)
  deallocate(sigma)
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
