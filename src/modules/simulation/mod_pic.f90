module mod_pic
  use iso_fortran_env, only: real64, int32, output_unit
  use mpi
  use mod_timer
  use mod_intro
  use mod_InputFiles            ! read_particle_conditions, find_conditions_file
  use mod_geometry             , only: Domain
  use mod_boundary            , only: Boundary, build_boundary
  use mod_io_hdf5             , only: save_array_h5, load_array_h5
  use mod_particles
  use mod_signals
  use mod_functionsText
  implicit none
  private
  public :: SimulationPIC

  contains

  subroutine SimulationPIC()
    implicit none

    ! ---- MPI / timing ----
    integer :: ierr_mpi, rank_local
    type(Timer) :: timer_simulation

    ! ---- Physics inputs (logging) ----
    real(real64) :: Te_eV, ne_m3, ngas_m3
    integer(int32) :: np_cell

    ! ---- Geometry & boundary ----
    type(Domain)   :: dom
    type(Boundary) :: bnd

    ! arrays for optional reload
    integer(int32), allocatable :: bcnd(:,:,:)
    real(real64)  , allocatable :: phi(:,:,:)
    real(real64) :: origin(3), spacing(3)
    integer(int32) :: ngh

    ! ---- Species / particles ----
    type(SpeciesTable) :: table_species
    type(Particle), allocatable :: parts(:)
    integer :: it

    ! ---- Init ----
    call MPI_Init(ierr_mpi)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_local, ierr_mpi)
    call install_signal_handlers()
    call timer_simulation%reset(); call timer_simulation%start()
    call print_introduction()

    ! ---- Geometry ----
    call dom%read_geometry()   ! reads nx,ny,nz & extents (meters)

    ! ---- Scalar conditions (for logging only) ----
    call read_particle_conditions(find_conditions_file(), Te_eV, ne_m3, ngas_m3, np_cell)
    if (rank_local == 0) then
      write(output_unit,*) 'Simulation parameters:'
      write(output_unit,*) '  Te (eV)   =', Te_eV
      write(output_unit,*) '  ne (m^-3) =', ne_m3
      write(output_unit,*) '  ngas (m^-3)=', ngas_m3
      write(output_unit,*) '  np_cell   =', abs(np_cell)
      write(output_unit,*) ' '
      call flush(output_unit)
    end if

    call build_boundary(bnd, dom)
 
 

    if (rank_local == 0) then
      write(output_unit,*) 'Boundary ready: size(bcnd)=', &
           size(bnd%bcnd,1),'x',size(bnd%bcnd,2),'x',size(bnd%bcnd,3)
      write(output_unit,*) 'Grid spacing (m): dx,dy,dz =', dom%dx1, dom%dx2, dom%dx3
      write(output_unit,*) ' '
      call flush(output_unit)
    end if

    ! ---- Species ----
    call load_species_from_file(table_species)

    ! ---- Particles (domain-aware) ----
    call build_particles_from_species_domain(table_species, dom, parts)

    if (rank_local==0) then
      write(output_unit,*) 'Allocations per species (N macro-particles):'
      do it = 1, size(parts)
        write(output_unit,'(A,1X,I0)') trim(parts(it)%name)//':', size(parts(it)%position1)
      end do
      write(output_unit,*) ' '
      call flush(output_unit)
    end if

    ! ---- Main loop (placeholder) ----
    do it = 0, 800
      if (rank_local==0 .and. mod(it,100)==0) then
        write(output_unit,*) 'Time step:', it
        call flush(output_unit)
      end if
      if (stop_requested) exit
      ! ... advance fields/particles/collisions here using dom, bnd, parts ...
    end do

    ! ---- Wrap up ----
    call timer_simulation%stop()
    if (rank_local==0) then
      write(output_unit,*) 'Elapsed time (s): ', timer_simulation%elapsed_s()
      call flush(output_unit)
    end if
    call MPI_Finalize(ierr_mpi)
  end subroutine SimulationPIC

end module mod_pic
