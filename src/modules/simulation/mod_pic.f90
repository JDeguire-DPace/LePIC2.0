module mod_pic
  use iso_fortran_env, only: real64, int32, int64, output_unit
  use mpi
  use mod_timer
  use mod_intro
  use mod_InputFiles
  use mod_particles
  use mod_signals
  implicit none
  private
  public :: SimulationPIC

contains

  subroutine SimulationPIC()
    implicit none
    ! ---- declarations FIRST ----
    integer :: ierr_mpi, rank_local
    type(Timer) :: timer_simulation
    real(real64) :: Te_eV, ne_m3, ngas_m3
    integer(int32) :: np_cell
    type(SpeciesTable) :: table_species
    type(Particle), allocatable :: parts(:)
    integer :: it

    ! ---- MPI init ----
    

    call MPI_Init(ierr_mpi)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_local, ierr_mpi)

    ! ---- enable Ctrl+C / SIGTERM graceful stop ----
    call install_signal_handlers()
    ! ---- timer ----
    call timer_simulation%reset()
    call timer_simulation%start()
    
    call print_introduction()

    ! ---- read inputs & species ----
    call read_particle_conditions(find_conditions_file(), Te_eV, ne_m3, ngas_m3, np_cell)
    call load_species_from_file(table_species)
    
    ! ---- (temporary) domain size on this rank ----
    ! Replace with your real geometry / local-domain once available.

    ! ---- allocate per-species particle arrays ----
    call build_particles_from_species(table_species, parts)

    if (rank_local==0) then
      write(output_unit,*) 'Simulation Parameters:'
      write(output_unit,*) '  Te (eV)=', Te_eV
      write(output_unit,*) '  ne (m^-3)=', ne_m3
      write(output_unit,*) '  ngas (m^-3)=', ngas_m3
      write(output_unit,*) '  np_cell =', np_cell
      write(output_unit,*) ' '
      write(output_unit,*) 'Allocations per species (N macro-particles):'
      do it = 1, size(parts)
        write(output_unit,'(A,1X,I0)') trim(parts(it)%name)//':', size(parts(it)%position1)
      end do
      write(output_unit,*) ' '
    end if

    ! ---- main loop (placeholder) ----
    do it = 0, 1000000
       if (rank_local==0 .and. mod(it,100000)==0) then
          write(output_unit,*) 'Time step:', it
       end if

       ! graceful interrupt check (Ctrl+C / SIGTERM)
       if (stop_requested) then
          call timer_simulation%stop()
          if (rank_local==0) then
             write(output_unit,*) 'Interrupted at step ', it
             write(output_unit,*) 'Elapsed time (s): ', timer_simulation%elapsed_s()
             call flush(output_unit)
          end if
          call MPI_Finalize(ierr_mpi)
          return   ! avoid double-finalize
       end if

       ! ... advance particles/fields here ...
    end do

    ! ---- normal end ----
    call timer_simulation%stop()
    if (rank_local==0) then
      write(output_unit,*) 'Total simulation time (s): ', timer_simulation%elapsed_s()
      call flush(output_unit)
    end if
    call MPI_Finalize(ierr_mpi)
    
  end subroutine SimulationPIC

end module mod_pic
