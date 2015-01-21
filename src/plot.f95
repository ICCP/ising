module plot
  use plplot
  implicit none
  private

  public plot_init, plot_close, plot_lattice

contains
  subroutine plot_init(sizex,sizey)
    integer, intent(in) :: sizex,sizey

    !You can find default colors at
    !http://plplot.sourceforge.net/docbook-manual/plplot-html-5.9.9/plcol0.html

    call plscol0(0, 255, 255, 255)  ! white
    call plscol0(1, 255, 0, 0)      ! red
    call plscol0(2, 255, 77, 0)     ! orange
    call plscol0(3, 255, 255, 0)    ! yellow
    call plscol0(4, 0, 255, 0)      ! green
    call plscol0(5, 0, 0, 255)      ! blue
    call plscol0(6, 0, 255, 255)    ! cyan
    call plscol0(7, 255, 0, 255)    ! magenta
    call plscol0(8, 128, 128, 128)  ! grey
    call plscol0(9, 0, 0, 0)        ! black

    call plsetopt("geometry","1000x1000") ! Change?
    call plsdev("png")
    call plsfam(1,1,100000)
    call plsfnam("%n.png")
    call plinit()
    call pladv(0)
    call plvpas(0d0, 1d0, 0d0, 1d0, 1d0)
    call plwind(5d-1, sizex+5d-1, 5d-1, sizey+5d-1)
!    call plbox("",1d0,0,"",1d0,0)
!    call plenv(5d-1, sizex+5d-1, 5d-1, sizey+5d-1, 1, -2)
  end subroutine plot_init

  subroutine plot_close()
    call plspause(.false.)
    !call plend()
  end subroutine plot_close

  subroutine plot_lattice(lattice)
    ! Draw a rectangular grid of up- and down-arrows corresponding to Ising
    ! states. Because this redraws the entire screen, it is *very* slow for
    ! Metropolis models - instead, consider a routine that only redraws the
    ! sites that change.
    integer, intent(in) :: lattice(:,:)
    integer :: i, j
    real(8) :: x, y
    !Assumes 1 corresponds to up-spin, -1 to down-spin.

    call plclear()
    do i = 1, size(lattice, 1)
      do j = 1, size(lattice, 2)
        x = i; y = j
        if (lattice(i, j) .eq. 1) then
          call plcol0(8)        ! grey
          call plot_square(x,y)
        else ! Comment out the else clauses to speed drawing
          call plcol0(9)        ! black
          call plot_square(x,y)
        end if
      end do
    end do

    call plflush()
  end subroutine plot_lattice

  subroutine plot_square(i,j)
    real(8), intent(in) :: i,j
    real(8)             :: x(4),y(4)

    x = (/i-.5,i-.5,i+.5,i+.5/)
    y = (/j-.5,j+.5,j+.5,j-.5/)

    call plfill(x,y)
  end subroutine
end module plot
