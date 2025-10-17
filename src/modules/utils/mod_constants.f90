module mod_constants
    
    use iso_fortran_env
    
    implicit none
    private


    real(real64), parameter, public :: PI = 3.14159265358979323846_real64
    real(real64), parameter, public :: elementary_charge = 1.602176634e-19_real64 ! Coulombs
    real(real64), parameter, public :: k_Boltzmann = 1.380649e-23_real64 !  m2 kg s-2 K-1
    
    real(real64), parameter, public :: epsilon0 = 8.854187817e-12_real64 ! C2 kg-1 m-3 s2
    real(real64), parameter, public :: mu0 = 1.2566370614e-6_real64 ! kg m A-2 s-2
    real(real64), parameter, public :: speed_of_light = 2.99792458e8_real64 ! m s-1

    real(real64), parameter, public :: amu = 1.66053906660e-27_real64 ! kg
    real(real64), parameter, public :: mass_electron = 9.1093837015e-31_real64 ! kg  

end module mod_constants