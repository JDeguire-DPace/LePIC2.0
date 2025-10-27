module mod_intro
        use iso_fortran_env
        use mpi
    implicit none
    private

    public :: print_introduction, print_license, draw_lepic, initialize_mpi

    contains

    subroutine draw_lepic()
        implicit none
        write(*,*) '                                      '
        write(*,'(a)') '                                                                       ⣶⣄⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀'
        write(*,'(a)') ' █████                ███████████  █████   █████████              ⠀⠀⠀⠀⢸⣿⣿⣷⣴⣿⡄⠀⠀⠀⠀⠀⢀⡀⠀⠀⠀'
        write(*,'(a)') '░░███                ░░███░░░░░███░░███   ███░░░░░███     ███    ⠀⠀⠀⠀⠰⣶⣾⣿⣿⣿⣿⣿⡇⠀⢠⣷⣤⣶⣿⡇⠀⠀⠀⠀⠀'
        write(*,'(a)') ' ░███         ██████  ░███    ░███ ░███  ███     ░░░     ░███     ⠀⠀  ⠙⣿⣿⣿⣿⣿⣿⣿⣀⣿⣿⣿⣿⣿⣧⣀⠀⠀⠀'
        write(*,'(a)') ' ░███        ███░░███ ░██████████  ░███ ░███          ███████████  ⣷⣦⣀⠘⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀'
        write(*,'(a)') ' ░███       ░███████  ░███░░░░░░   ░███ ░███         ░░░░░███░░░ ⢲⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠁⠀⠀'
        write(*,'(a)') ' ░███      █░███░░░   ░███         ░███ ░░███     ███    ░███   ⠀ ⠙⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⠁⠀'
        write(*,'(a)') ' ███████████░░██████  █████        █████ ░░█████████     ░░░    ⠀⠀⠚⠻⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠿⠿⠂⠀⠀'
        write(*,'(a)') '░░░░░░░░░░░  ░░░░░░  ░░░░░        ░░░░░   ░░░░░░░░░           ⠀⠀⠀⠀⠀⠀⠉⠙⢻⣿⣿⡿⠛⠉⡇⠀⠀⠀⠀⠀⠀⠀'
        write(*,'(a)') '                                                                        ⠋⠁⠀⠀⠸⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀'
        write(*,'(a)') '                                Laplace Explicit Particle-In-Cell code   '
        write(*,'(a)') '                                        (Modular Object-Oriented Edition)'
        write(*,'(a)') '                                     '

        return
    end subroutine draw_lepic

    subroutine print_license()
            implicit none
  
            write(*,'(A)') ' ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜     3D explicit parallel particle in cell code                             ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜   LePIC 3D algorithm version: 3.3.6.2                                      ⚜'
            write(*,'(A)') '⚜   Last modifications: Sep/2024                                             ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜   Copyright or © or Copr. Gwenael Fubiani (2022/11/22)                     ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜   gwenael.fubiani@cnrs.fr                                                  ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜   This software is a computer program whose purpose is to model low        ⚜'
            write(*,'(A)') '⚜   temperature plasmas with a Particle-In-Cell algorihm in 3D-3V            ⚜'
            write(*,'(A)') '⚜   dimensions.                                                              ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜   This software is governed by the CeCILL-C license under French law and   ⚜'
            write(*,'(A)') '⚜   abiding by the rules of distribution of free software.  You can  use,    ⚜'
            write(*,'(A)') '⚜   modify and/ or redistribute the software under the terms of the CeCILL-C ⚜'
            write(*,'(A)') '⚜   license as circulated by CEA, CNRS and INRIA at the following URL        ⚜'
            write(*,'(A)') '⚜   "http://www.cecill.info".                                                ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜   As a counterpart to the access to the source code and  rights to copy,   ⚜'
            write(*,'(A)') '⚜   modify and redistribute granted by the license, users are provided only  ⚜'
            write(*,'(A)') "⚜   with a limited warranty and the software's author, the holder of the     ⚜"
            write(*,'(A)') '⚜   economic rights,  and the successive licensors  have only  limited       ⚜'
            write(*,'(A)') '⚜   liability.                                                               ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') "⚜   In this respect, the user's attention is drawn to the risks associated   ⚜"
            write(*,'(A)') '⚜   with loading,  using,  modifying and/or developing or reproducing the    ⚜'
            write(*,'(A)') '⚜   software by the user in light of its specific status of free software,   ⚜'
            write(*,'(A)') '⚜   that may mean  that it is complicated to manipulate,  and  that  also    ⚜'
            write(*,'(A)') '⚜   therefore means  that it is reserved for developers  and  experienced    ⚜'
            write(*,'(A)') '⚜   professionals having in-depth computer knowledge. Users are therefore    ⚜'
            write(*,'(A)') "⚜   encouraged to load and test the software's suitability as regards their  ⚜"
            write(*,'(A)') '⚜   requirements in conditions enabling the security of their systems and/or ⚜'
            write(*,'(A)') '⚜   data to be ensured and,  more generally, to use and operate it in the    ⚜'
            write(*,'(A)') '⚜   same conditions as regards security.                                     ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') '⚜   The fact that you are presently reading this means that you have had     ⚜'
            write(*,'(A)') '⚜   knowledge of the CeCILL-C license and that you accept its terms.         ⚜'
            write(*,'(A)') '⚜                                                                            ⚜'
            write(*,'(A)') ' ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜ ⚜'
            write(*,'(A)') '                                                  '
            write(*,'(A)') '                                                  '
            
            return
    end subroutine print_license

    subroutine print_introduction()
        implicit none
        call draw_lepic()
        call print_license()
        return
    end subroutine print_introduction

    subroutine initialize_mpi()
        implicit none
        integer ierr,mpi_rank,nproc_mpi
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,mpi_rank,ierr)
        call mpi_comm_size(mpi_comm_world,nproc_mpi,ierr)
        write(*,*) ' MPI initialized with ', nproc_mpi, ' processes. This is process ', mpi_rank
        return
    end subroutine initialize_mpi

end module mod_intro