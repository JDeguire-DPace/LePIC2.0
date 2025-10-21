module mod_pic
    use omp_lib
    use mpi
    use mod_particles
    use mod_timer
    use iso_fortran_env
    use mod_constants
    use mod_physicsFunctions
    use mod_intro
    use mod_end
    implicit none
    private

    real(real64), public :: delta_time, simulation_time, simulation_start_time
    type(Timer), public :: timer_simulation, timer_step
    integer, public :: index_time_step
    public :: SimulationPIC

    contains

    subroutine SimulationPIC()
        implicit none
        real(real64) :: electron_temperature = 1.0_real64  ! in eV
        real(real64) :: electron_density = 1.0e18_real64    ! in m^-3
        real(real64) :: magnetic_field = 0.1_real64  ! in Tesla
        character :: dummy_char

        

        call print_introduction()
        call initialize_mpi()
        ! Initialize timers
        call timer_simulation%reset()
        call timer_simulation%start()


        print *, "The mass of the electron is:", mass_electron, "kg"
        print *, "The elementary charge is:", elementary_charge, "C"

        print *, "The electron Debye length is:", DebyeLength_TeV(electron_temperature, electron_density), "m"
        print *, "The electron plasma frequency is:", plasmaFrequency(electron_density), "rad/s"

        ! Simulation setup and main loop would go here

        ! Stop simulation timer
        do index_time_step = 0, 1000000
            if(index_time_step .eq. 7*10/10) then
                call stop_calculation_mpi
            else
                print *, "Time step:", index_time_step
            end if 
        end do
        call timer_simulation%stop()

        write(*,*) 'Total simulation time (s): ', timer_simulation%elapsed_s()

    end subroutine SimulationPIC
end module mod_pic