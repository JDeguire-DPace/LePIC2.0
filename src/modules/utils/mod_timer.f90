module mod_timer
    use iso_fortran_env, only: int64, real64
    implicit none
    private
    public :: Timer, wall_clock_resolution

    ! Clock characteristics
    integer(int64) :: g_rate = 0_int64
    integer(int64) :: g_cmax = huge(0_int64)

    type :: Timer
        private
        integer(int64) :: c_start = 0_int64
        integer(int64) :: c_end   = 0_int64
        logical        :: running = .false.
    contains
        procedure :: start      => timer_start
        procedure :: stop       => timer_stop
        procedure :: reset      => timer_reset
        procedure :: elapsed_s  => timer_elapsed_s
        procedure :: elapsed_ms => timer_elapsed_ms
        procedure :: is_running => timer_is_running
    end type Timer

contains

    subroutine init_clock()
        integer :: cnt32, rate32, cmax32
        call system_clock(cnt32, rate32, cmax32)
        if (g_rate == 0_int64) g_rate = int(rate32, int64)
        if (cmax32 > 0)        g_cmax = int(cmax32, int64)
    end subroutine init_clock

    function wall_clock_resolution() result(dt)
        real(real64) :: dt
        if (g_rate == 0_int64) call init_clock()
        dt = 1.0_real64 / real(g_rate, real64)
    end function wall_clock_resolution

    pure function ticks_delta(c1, c2) result(dticks)
        integer(int64), intent(in) :: c1, c2
        integer(int64)             :: dticks
        if (g_cmax > 0_int64 .and. c2 < c1) then
            dticks = (g_cmax - c1) + c2 + 1_int64
        else
            dticks = c2 - c1
        end if
    end function ticks_delta

    subroutine timer_start(self)
        class(Timer), intent(inout) :: self
        integer :: cnow
        if (g_rate == 0_int64) call init_clock()
        call system_clock(cnow)
        self%c_start = int(cnow, int64)
        self%c_end   = self%c_start
        self%running = .true.
    end subroutine timer_start

    subroutine timer_stop(self)
        class(Timer), intent(inout) :: self
        integer :: cnow
        if (.not. self%running) return
        call system_clock(cnow)
        self%c_end   = int(cnow, int64)
        self%running = .false.
    end subroutine timer_stop

    subroutine timer_reset(self)
        class(Timer), intent(inout) :: self
        self%c_start = 0_int64
        self%c_end   = 0_int64
        self%running = .false.
    end subroutine timer_reset

    function timer_elapsed_s(self) result(sec)
        class(Timer), intent(in) :: self
        real(real64)             :: sec
        integer :: ccur32
        integer(int64) :: ccur

        if (g_rate == 0_int64) call init_clock()

        if (self%running) then
            call system_clock(ccur32)
            ccur = int(ccur32, int64)
        else
            ccur = self%c_end
        end if

        sec = real(ticks_delta(self%c_start, ccur), real64) / real(g_rate, real64)
    end function timer_elapsed_s

    function timer_elapsed_ms(self) result(ms)
        class(Timer), intent(in) :: self
        real(real64)             :: ms
        ms = 1000.0_real64 * self%elapsed_s()
    end function timer_elapsed_ms

    pure function timer_is_running(self) result(flag)
        class(Timer), intent(in) :: self
        logical :: flag
        flag = self%running
    end function timer_is_running

end module mod_timer
