module mod_InputFiles
    use iso_fortran_env
    use mod_functionsText
    implicit none
    private

    public :: find_particles_file, find_conditions_file, find_geometry_file, find_boundary_file, find_magneticfield_file
    contains

    character(len=512) function find_particles_file() result(path)
        character(len=512) :: filename, line
        integer :: ios
        logical :: found
        path  = ''
        found = .false.

        ! (Portable; removed GETCWD to avoid IFPORT dependency)
        ! write(*,*) 'Working directory (not printed to avoid IFPORT).'

        open(unit=10, file='./Utils/PathSimulation.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open ./Utils/PathSimulation.txt'
            return
        end if
        do
            read(10,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Particles =') > 0) then
                call parse_value_after_equal(line, filename)
                path  = trim(filename)
                found = .true.
                exit
            end if
        end do
        close(10)
        if (.not. found) then
            write(*,*) 'Warning: no "Particles =" entry found in PathSimulation.txt'
        end if
    end function find_particles_file

    character(len=512) function find_conditions_file() result(path)
        character(len=512) :: filename, line
        integer :: ios
        logical :: found
        path  = ''
        found = .false.

        ! (Portable; removed GETCWD to avoid IFPORT dependency)
        ! write(*,*) 'Working directory (not printed to avoid IFPORT).'

        open(unit=10, file='./Utils/PathSimulation.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open ./Utils/PathSimulation.txt'
            return
        end if
        do
            read(10,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Conditions =') > 0) then
                call parse_value_after_equal(line, filename)
                path  = trim(filename)
                found = .true.
                exit
            end if
        end do
        close(10)
        if (.not. found) then
            write(*,*) 'Warning: no "Conditions =" entry found in PathSimulation.txt'
        end if
    end function find_conditions_file

    character(len=512) function find_geometry_file() result(path)
        character(len=512) :: filename, line
        integer :: ios
        logical :: found
        path  = ''
        found = .false.

        ! (Portable; removed GETCWD to avoid IFPORT dependency)
        ! write(*,*) 'Working directory (not printed to avoid IFPORT).'

        open(unit=10, file='./Utils/PathSimulation.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open ./Utils/PathSimulation.txt'
            return
        end if
        do
            read(10,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Geometry =') > 0) then
                call parse_value_after_equal(line, filename)
                path  = trim(filename)
                found = .true.
                exit
            end if
        end do
        close(10)
        if (.not. found) then
            write(*,*) 'Warning: no "Geometry =" entry found in PathSimulation.txt'
        end if
    end function find_geometry_file

    character(len=512) function find_boundary_file() result(path)
        character(len=512) :: filename, line
        integer :: ios
        logical :: found
        path  = ''
        found = .false.

        ! (Portable; removed GETCWD to avoid IFPORT dependency)
        ! write(*,*) 'Working directory (not printed to avoid IFPORT).'

        open(unit=10, file='./Utils/PathSimulation.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open ./Utils/PathSimulation.txt'
            return
        end if
        do
            read(10,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Boundary =') > 0) then
                call parse_value_after_equal(line, filename)
                path  = trim(filename)
                found = .true.
                exit
            end if
        end do
        close(10)
        if (.not. found) then
            write(*,*) 'Warning: no "Boundary =" entry found in PathSimulation.txt'
        end if
    end function find_boundary_file


    character(len=512) function find_magneticfield_file() result(path)
        character(len=512) :: filename, line
        integer :: ios
        logical :: found
        path  = ''
        found = .false.

        ! (Portable; removed GETCWD to avoid IFPORT dependency)
        ! write(*,*) 'Working directory (not printed to avoid IFPORT).'

        open(unit=10, file='./Utils/PathSimulation.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open ./Utils/PathSimulation.txt'
            return
        end if
        do
            read(10,'(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Field =') > 0) then
                call parse_value_after_equal(line, filename)
                path  = trim(filename)
                found = .true.
                exit
            end if
        end do
        close(10)
        if (.not. found) then
            write(*,*) 'Warning: no "Field =" entry found in PathSimulation.txt'
        end if
    end function find_magneticfield_file


end module mod_InputFiles