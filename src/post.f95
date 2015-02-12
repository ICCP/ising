program post
  use plot

  integer              :: nx,ny
  integer, allocatable :: sigma(:,:)
  integer              :: uin
  integer              :: iostatus

  uin = 50

  open(unit=uin,file='out',form='unformatted')
  read (unit=uin) nx,ny
  allocate(sigma(nx,ny))
  call plot_init(nx,ny)

  do
    read (unit=uin,iostat=iostatus) sigma
    if (iostatus>0) then
      write (*,*) "problem!"
      call exit(iostatus)
    else if (iostatus<0) then
      exit
    else
      call plbop()
      call plot_lattice(sigma)
      call pleop()
    end if
  end do

  call plot_close()
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
