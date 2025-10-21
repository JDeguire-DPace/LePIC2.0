module mod_particles
    use iso_fortran_env
    use mpi
    use mod_constants
    use mod_physicsFunctions
    use mod_functionsText
    USE IFPORT

    implicit none
    private
    public :: Particle, find_particle_file

    type :: Particle
        real(real64), allocatable :: position1(:), position2(:), position3(:)
        real(real64), allocatable :: velocity1(:), velocity2(:), velocity3(:)
        real(real64) :: charge, mass, weight

    end type Particle

    ! interface Particle
    !     module procedure find_particle_file
    ! end interface Particle

    contains
        character function find_particle_file() result(path)
            character(len=256) :: filename, current_path
            character(len=512) :: line, key
            integer :: ios,status
            logical :: found = .false.
            
            status = GETCWD(current_path)
            write(*,*) 'Current working directory: ', trim(current_path)

            open(unit=10, file='./Utils/PathSimulation.txt', status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error: cannot open PathSimulation.txt'
                return
            end if

            do
                read(10,'(A)', iostat=ios) line
                if (ios == iostat_end) exit
                if (index(line, 'Particles =') > 0) then
                    ! extract text after '='
                    line = adjustl(trim(line))
                    call parse_path(line, filename)
                    found = .true.
                    exit
                end if
            end do
            close(10)
            if (.not. found) then
                print *, 'Warning: no "Particles =" entry found in PathSimulation.txt'
            else
                print *, 'Particle file found: ', trim(filename)
                path = trim(filename)
            end if
        end function find_particle_file
end module mod_particles