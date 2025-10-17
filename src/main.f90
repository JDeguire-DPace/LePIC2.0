program main
  use mod_physicsFunctions
  use iso_fortran_env
  use mod_constants

  implicit none
  real(real64) :: electron_temperature = 1.0_real64  ! in eV
  real(real64) :: electron_density = 1.0e18_real64    ! in m^-3
  real(real64) :: magnetic_field = 0.1_real64

  

  ! Print numbers from 1 to 10
  print *, "The mass of the electron is:", mass_electron, "kg"
  print *, "The elementary charge is:", elementary_charge, "C"

  print *, "The electron Debye length is:", DebyeLength_TeV(electron_temperature, electron_density), "m"
  print *, "The electron plasma frequency is:", plasmaFrequency(electron_density), "rad/s"

end program main