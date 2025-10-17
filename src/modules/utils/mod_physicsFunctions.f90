module mod_physicsFunctions
    use iso_fortran_env
    use mod_constants
    implicit none
    private

    public :: DebyeLength_TKelvin, DebyeLength_TeV, plasmaFrequency, Gyrofrequency, Gyrroradius, &
              ThermalVelocity_TKelvin, ThermalVelocity_TeV
    contains



    function DebyeLength_TKelvin(temperature_K, density_m3) result(lambda_D)
        real(real64), intent(in) :: temperature_K  ! Temperature in Kelvin
        real(real64), intent(in) :: density_m3      ! Electron density in m^-3
        real(real64) :: lambda_D                    ! Debye length in meters

        lambda_D = sqrt((epsilon0 * k_Boltzmann * temperature_K) / &
                        (density_m3 * elementary_charge**2))
    end function DebyeLength_TKelvin

    function DebyeLength_TeV(temperature_eV, density_m3) result(lambda_D)
        real(real64), intent(in) :: temperature_eV  ! Temperature in electronvolts
        real(real64), intent(in) :: density_m3      ! Electron density in m^-3
        real(real64) :: lambda_D                    ! Debye length in meters

        lambda_D = sqrt((epsilon0 * temperature_eV * elementary_charge) / &
                        (density_m3 * elementary_charge**2))
    end function DebyeLength_TeV

    function plasmaFrequency(density_m3) result(omega_p)
        real(real64), intent(in) :: density_m3  ! Electron density in m^-3
        real(real64) :: omega_p                 ! Plasma frequency in radians per second

        omega_p = sqrt((density_m3 * elementary_charge**2) / &
                       (mass_electron * epsilon0))
    end function plasmaFrequency

    function Gyrofrequency(charge, magnetic_field) result(omega_c)
        real(real64), intent(in) :: magnetic_field, charge  ! Magnetic field in Tesla, charge in Coulombs
        real(real64) :: omega_c                       ! Gyrofrequency in radians per second

        omega_c = abs(charge) * magnetic_field / mass_electron
    end function Gyrofrequency

    function Gyrroradius(velocity_perp, mass, charge, magnetic_field) result(r_L)
        real(real64), intent(in) :: velocity_perp, magnetic_field, mass, charge  ! Perpendicular velocity in m/s, magnetic field in Tesla, mass in kg, charge in Coulombs
        real(real64) :: r_L                                   ! Gyrroradius in meters

        r_L = (mass * velocity_perp) / (abs(charge) * magnetic_field)
    end function Gyrroradius

    function ThermalVelocity_TKelvin(temperature_K, mass) result(v_th)
        real(real64), intent(in) :: temperature_K, mass  ! Temperature in Kelvin, mass in kg
        real(real64) :: v_th                             ! Thermal velocity in m/s

        v_th = sqrt((2.0_real64 * k_Boltzmann * temperature_K) / mass)
    end function ThermalVelocity_TKelvin

    function ThermalVelocity_TeV(temperature_eV, mass) result(v_th)
        real(real64), intent(in) :: temperature_eV, mass  ! Temperature in electronvolts, mass in kg
        real(real64) :: v_th                             ! Thermal velocity in m/s

        v_th = sqrt((2.0_real64 * temperature_eV * elementary_charge) / mass)
    end function ThermalVelocity_TeV



end module mod_physicsFunctions