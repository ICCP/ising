module ising_routines
  
  contains
  pure function energy(lattice)
  
    use global
  
    implicit none

    integer, intent(inout) :: lattice(latticeSize,latticeSize)

    print*,'latice:',lattice
  end function
end module
