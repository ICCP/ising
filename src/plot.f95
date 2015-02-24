module plot
  use plplot
  implicit none
  private

  public plot_init, plot_close, plot_lattice

contains

  subroutine plot_init(latticeSize)
    integer, intent(in) :: latticeSize
    call plsdev("xcairo")
    call plinit()

    !You can find default colors at
    !http://plplot.sourceforge.net/docbook-manual/plplot-html-5.9.9/plcol0.html

    !call plscol0(0, 255, 255, 255)  ! white
    !call plscol0(1, 255, 0, 0)      ! red
    !call plscol0(2, 255, 77, 0)     ! orange
    !call plscol0(3, 255, 255, 0)    ! yellow
    !call plscol0(4, 0, 255, 0)      ! green
    !call plscol0(5, 0, 0, 255)      ! blue
    !call plscol0(6, 0, 255, 255)    ! cyan
    !call plscol0(7, 255, 0, 255)    ! magenta
    !call plscol0(8, 128, 128, 128)  ! gray
    !call plscol0(9, 0, 0, 0)        ! black

    call plenv(0d0, latticeSize + 1d0, 0d0, latticeSize + 1d0, 0, 0)
  end subroutine plot_init

  subroutine plot_close()
    call plspause(.false.)
    call plend()
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
          call plcol0(1)            !default red
          call plpoin([x], [y], 30) !30 denotes up-arrow glyph
        else ! Comment out the else clauses to speed drawing
          call plcol0(11)           !default cyan
          call plpoin([x], [y], 31) !31 denotes down-arrow glyph
        end if
      end do
    end do

    call plflush()
  end subroutine plot_lattice

end module plot
