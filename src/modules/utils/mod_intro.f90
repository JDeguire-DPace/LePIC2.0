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
        write(*,*) '                                      '
        write(*,'(a)') ' █████                ███████████  █████   █████████             '
        write(*,'(a)') '░░███                ░░███░░░░░███░░███   ███░░░░░███     ███    '
        write(*,'(a)') ' ░███         ██████  ░███    ░███ ░███  ███     ░░░     ░███    '
        write(*,'(a)') ' ░███        ███░░███ ░██████████  ░███ ░███          ███████████'
        write(*,'(a)') ' ░███       ░███████  ░███░░░░░░   ░███ ░███         ░░░░░███░░░ '
        write(*,'(a)') ' ░███      █░███░░░   ░███         ░███ ░░███     ███    ░███    '
        write(*,'(a)') ' ███████████░░██████  █████        █████ ░░█████████     ░░░     '
        write(*,'(a)') '░░░░░░░░░░░  ░░░░░░  ░░░░░        ░░░░░   ░░░░░░░░░              '
        write(*,*) '                                      '
        write(*,*) '                                Laplace Explicit Particle-In-Cell code'
        write(*,*) '                                        (Modular Object-Oriented Edition)'
        write(*,*) '                                      '

        return
    end subroutine draw_lepic

    subroutine print_license()
            implicit none
  
            write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) '+                                                                          +'
            write(*,*) '+   3D explicit parallel particle in cell code                             +'
            write(*,*) '+                                                                          +'
            write(*,*) '+ LePIC 3D algorithm version: 3.3.6.2                                      +'
            write(*,*) '+ Last modifications: Sep/2024                                             +'
            write(*,*) '+                                                                          +'
            write(*,*) '+ Copyright or © or Copr. Gwenael Fubiani (2022/11/22)                     +'
            write(*,*) '+                                                                          +'
            write(*,*) '+ gwenael.fubiani@cnrs.fr                                                  +'
            write(*,*) '+                                                                          +'
            write(*,*) '+ This software is a computer program whose purpose is to model low        +'
            write(*,*) '+ temperature plasmas with a Particle-In-Cell algorihm in 3D-3V            +'
            write(*,*) '+ dimensions.                                                              +'
            write(*,*) '+                                                                          +'
            write(*,*) '+ This software is governed by the CeCILL-C license under French law and   +'
            write(*,*) '+ abiding by the rules of distribution of free software.  You can  use,    +'
            write(*,*) '+ modify and/ or redistribute the software under the terms of the CeCILL-C +'
            write(*,*) '+ license as circulated by CEA, CNRS and INRIA at the following URL        +'
            write(*,*) '+ "http://www.cecill.info".                                                +'
            write(*,*) '+                                                                          +'
            write(*,*) '+ As a counterpart to the access to the source code and  rights to copy,   +'
            write(*,*) '+ modify and redistribute granted by the license, users are provided only  +'
            write(*,*) "+ with a limited warranty and the software's author, the holder of the     +"
            write(*,*) '+ economic rights,  and the successive licensors  have only  limited       +'
            write(*,*) '+ liability.                                                               +'
            write(*,*) '+                                                                          +'
            write(*,*) "+ In this respect, the user's attention is drawn to the risks associated   +"
            write(*,*) '+ with loading,  using,  modifying and/or developing or reproducing the    +'
            write(*,*) '+ software by the user in light of its specific status of free software,   +'
            write(*,*) '+ that may mean  that it is complicated to manipulate,  and  that  also    +'
            write(*,*) '+ therefore means  that it is reserved for developers  and  experienced    +'
            write(*,*) '+ professionals having in-depth computer knowledge. Users are therefore    +'
            write(*,*) "+ encouraged to load and test the software's suitability as regards their  +"
            write(*,*) '+ requirements in conditions enabling the security of their systems and/or +'
            write(*,*) '+ data to be ensured and,  more generally, to use and operate it in the    +'
            write(*,*) '+ same conditions as regards security.                                     +'
            write(*,*) '+                                                                          +'
            write(*,*) '+ The fact that you are presently reading this means that you have had     +'
            write(*,*) '+ knowledge of the CeCILL-C license and that you accept its terms.         +'
            write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) '                                                  '
            write(*,*) '                                                  '
            
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