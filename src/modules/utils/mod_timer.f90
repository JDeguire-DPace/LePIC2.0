module mod_timer
    use iso_fortran_env, only: int64, real64
    implicit none
    private
    public :: Timer, wall_clock_resolution

    ! Clock characteristics (set lazily in submodule)
    integer(int64) :: g_rate = 0_int64
    integer(int64) :: g_cmax = huge(0_int64)

    type :: Timer
        private
        integer(int64) :: c_start = 0_int64
        integer(int64) :: c_end   = 0_int64
        logical        :: running = .false.
    contains
        ! Bind public method names to specific procedures implemented in submodule
        procedure :: start      => timer_start
        procedure :: stop       => timer_stop
        procedure :: reset      => timer_reset
        procedure :: elapsed_s  => timer_elapsed_s
        procedure :: elapsed_ms => timer_elapsed_ms
        procedure :: is_running => timer_is_running
    end type Timer

    ! ---- Declare interfaces for all procedures implemented in submodule ----
    interface
        module subroutine init_clock()
        end subroutine

        module function wall_clock_resolution() result(dt)
            real(real64) :: dt
        end function

        module pure function ticks_delta(c1, c2) result(dticks)
            integer(int64), intent(in) :: c1, c2
            integer(int64)             :: dticks
        end function

        module subroutine timer_start(self)
            class(Timer), intent(inout) :: self
        end subroutine

        module subroutine timer_stop(self)
            class(Timer), intent(inout) :: self
        end subroutine

        module subroutine timer_reset(self)
            class(Timer), intent(inout) :: self
        end subroutine

        module function timer_elapsed_s(self) result(sec)
            class(Timer), intent(in) :: self
            real(real64)             :: sec
        end function

        module function timer_elapsed_ms(self) result(ms)
            class(Timer), intent(in) :: self
            real(real64)             :: ms
        end function

        module pure function timer_is_running(self) result(flag)
            class(Timer), intent(in) :: self
            logical                   :: flag
        end function
    end interface
end module mod_timer
