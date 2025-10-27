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

module mod_signals
  use, intrinsic :: iso_c_binding
  implicit none
  logical, volatile :: stop_requested = .false.
contains

  ! Handler must be bind(C) and take an int by VALUE
  subroutine sigint_handler(sig) bind(C)
    use, intrinsic :: iso_c_binding
    integer(c_int), value :: sig
    stop_requested = .true.
  end subroutine sigint_handler

  subroutine install_signal_handlers()
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_funptr) :: prev

    interface
      function c_signal(sig, handler) bind(C, name="signal") result(prev)
        import :: c_int, c_funptr
        integer(c_int), value :: sig
        type(c_funptr), value :: handler
        type(c_funptr) :: prev
      end function c_signal
    end interface

    ! 2 = SIGINT (Ctrl-C), 15 = SIGTERM (scancel/kill -15)
    prev = c_signal( 2_c_int, c_funloc(sigint_handler))
    prev = c_signal(15_c_int, c_funloc(sigint_handler))
    ! Optionally: 20 = SIGTSTP (Ctrl-Z) and 18 = SIGCONT, if you implement pause/resume
    ! prev = c_signal(20_c_int, c_funloc(sigint_handler))
    ! prev = c_signal(18_c_int, c_funloc(sigint_handler))
  end subroutine install_signal_handlers

end module mod_signals
