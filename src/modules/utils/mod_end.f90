module mod_end
    use mpi
    use mod_timer

    implicit none
    private

    public :: stop_calculation_mpi !, print_end_message

    contains
        subroutine stop_calculation_mpi
            !     ==============================================================
            !     VERSION:         0.1
            !     LAST MOD:      Oct/19
            !     MOD AUTHOR:    G. Fubiani
            !     COMMENTS: 
            !     -------------------------------------------------------------
            implicit none
            integer ierr

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            stop
            call MPI_Finalize(ierr)

            return
        end subroutine stop_calculation_mpi
end module mod_end