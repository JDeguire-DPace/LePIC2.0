module mod_pic
  use iso_fortran_env, only: real64, int32, output_unit
  use mpi
  use mod_timer
  use mod_intro
  use mod_InputFiles            ! read_particle_conditions, find_conditions_file
  use mod_geometry             , only: Domain
  use mod_boundary             , only: Boundary, build_boundary
  use mod_io_hdf5              , only: save_array_h5, load_array_h5
  use mod_particles
  use mod_signals
  use mod_functionsText
  ! NEW:
  use mod_MagneticField        , only: MagneticField
  implicit none
  private
  public :: SimulationPIC

contains

  subroutine SimulationPIC()
    implicit none
    ! ---- MPI / timing ----
    integer :: ierr_mpi, rank_local, nproc_mpi
    type(Timer) :: timer_simulation

    ! ---- Physics inputs (for logging only) ----
    real(real64) :: Te_eV, ne_m3, ngas_m3
    integer(int32) :: np_cell

    ! ---- Geometry & boundary ----
    type(Domain)   :: dom
    type(Boundary) :: bnd

    ! ---- Magnetic field ----
    type(MagneticField) :: bkgB
    real(real64), allocatable :: b1(:,:,:), b2(:,:,:), b3(:,:,:)
    real(real64) :: origin_dom(3), spacing_dom(3)
    real(real64) :: origin_B(3)  , spacing_B(3)
    integer(int32) :: ngh_dom

    ! ---- Species / particles ----
    type(SpeciesTable) :: table_species
    type(Particle), allocatable :: parts(:)
    integer :: it

    ! ------------------------------------------------------------------
    ! Init
    ! ------------------------------------------------------------------
    call MPI_Init(ierr_mpi)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_local, ierr_mpi)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc_mpi, ierr_mpi)
    call install_signal_handlers()
    call timer_simulation%reset(); call timer_simulation%start()
    if (rank_local == 0) call print_introduction()

    ! ------------------------------------------------------------------
    ! Geometry (nx,ny,nz, extents in meters)
    ! ------------------------------------------------------------------
    call dom%read_geometry()

    if (rank_local == 0) then
      write(output_unit,*) 'Geometry:'
      write(output_unit,'(A,3(I0,1X))') '  ncell = ', dom%ncell_x1, dom%ncell_x2, dom%ncell_x3
      write(output_unit,'(A,3(ES12.4,1X))') '  dx,dy,dz (m) = ', dom%dx1, dom%dx2, dom%dx3
      call flush(output_unit)
    end if

    origin_dom  = [0.0_real64, 0.0_real64, 0.0_real64]
    spacing_dom = [dom%dx1, dom%dx2, dom%dx3]
    ngh_dom     = dom%ncell_ghost

    ! ------------------------------------------------------------------
    ! Scalar particle conditions (for logging)
    ! ------------------------------------------------------------------
    call read_particle_conditions(find_conditions_file(), Te_eV, ne_m3, ngas_m3, np_cell)
    if (rank_local == 0) then
      write(output_unit,*) 'Plasma/gas (from conditions.inp):'
      write(output_unit,'(A,ES12.4)') '  Te (eV)   = ', Te_eV
      write(output_unit,'(A,ES12.4)') '  ne (m^-3) = ', ne_m3
      write(output_unit,'(A,ES12.4)') '  ngas (m^-3)= ', ngas_m3
      write(output_unit,'(A,I0)')     '  np_cell   = ', abs(np_cell)
      call flush(output_unit)
    end if

    ! ------------------------------------------------------------------
    ! Boundary (build bcnd/phi from geometry+boundary.inp) and save
    ! ------------------------------------------------------------------
    call build_boundary(bnd, dom)

    if (rank_local == 0) then
      write(output_unit,*) 'Boundary built:'
      write(output_unit,'(A,3(I0,1X))') '  size(bcnd) = ', size(bnd%bcnd,1), size(bnd%bcnd,2), size(bnd%bcnd,3)
      call flush(output_unit)
    end if

    call save_array_h5("./Outputs/state.h5","bcnd", bnd%bcnd, origin_dom, spacing_dom, ngh_dom, replace=.true.)
    call save_array_h5("./Outputs/state.h5","phi" , bnd%phi , origin_dom, spacing_dom, ngh_dom, replace=.true.)
    call save_array_h5("./Outputs/state.h5","wall_type", bnd%wall_type, replace=.true.)
    call save_array_h5("./Outputs/state.h5","potential_from_boundary", bnd%potential_from_boundary, replace=.true.)

    ! ------------------------------------------------------------------
    ! Magnetic field (background map or analytic)
    ! ------------------------------------------------------------------
    call bkgB%read_from_file( trim(find_magneticfield_file()) )
    call bkgB%allocate_arrays()
    call bkgB%build_field()

    ! spacing/origin for B-grid (can differ from domain)
    spacing_B = [ (bkgB%x1_max - bkgB%x1_min) / real(max(1,bkgB%nx1-1), real64), &
                  (bkgB%x2_max - bkgB%x2_min) / real(max(1,bkgB%nx2-1), real64), &
                  (bkgB%x3_max - bkgB%x3_min) / real(max(1,bkgB%nx3-1), real64) ]
    origin_B  = [ bkgB%x1_min, bkgB%x2_min, bkgB%x3_min ]

    ! save B-components separately (HDF5 commonly expects rank-3)
    allocate(b1(bkgB%nx1, bkgB%nx2, bkgB%nx3))
    allocate(b2(bkgB%nx1, bkgB%nx2, bkgB%nx3))
    allocate(b3(bkgB%nx1, bkgB%nx2, bkgB%nx3))
    b1 = bkgB%B(1,:,:,:)
    b2 = bkgB%B(2,:,:,:)
    b3 = bkgB%B(3,:,:,:)

    call save_array_h5("./Outputs/state.h5","B1", b1, origin_B, spacing_B, 0, replace=.true.)
    call save_array_h5("./Outputs/state.h5","B2", b2, origin_B, spacing_B, 0, replace=.false.)
    call save_array_h5("./Outputs/state.h5","B3", b3, origin_B, spacing_B, 0, replace=.false.)

    if (rank_local == 0) then
      write(output_unit,*) 'Magnetic field loaded:'
      write(output_unit,'(A,3(I0,1X))') '  grid = ', bkgB%nx1, bkgB%nx2, bkgB%nx3
      write(output_unit,'(A,3(ES12.4,1X))') '  spacing_B (m) = ', spacing_B
      write(output_unit,'(A,3(ES12.4,1X))') '  origin_B  (m) = ', origin_B
      write(output_unit,'(A,I0)') '  direction component populated = ', bkgB%direction
      call flush(output_unit)
    end if

    ! ------------------------------------------------------------------
    ! Species & particles (no domain decomposition yet)
    ! ------------------------------------------------------------------
    call load_species_from_file(table_species)
    call build_particles_from_species_domain(table_species, dom, parts)

    if (rank_local == 0) then
      write(output_unit,*) 'Particles allocated per species (macro count):'
      do it = 1, size(parts)
        write(output_unit,'(A,1X,I0)') trim(parts(it)%name)//':', size(parts(it)%position1)
      end do
      call flush(output_unit)
    end if

    ! ------------------------------------------------------------------
    ! Main loop (placeholder for now)
    ! ------------------------------------------------------------------
    do it = 0, 200
      if (rank_local==0 .and. mod(it,50)==0) then
        write(output_unit,'(A,I0)') 'timestep = ', it
        call flush(output_unit)
      end if
      if (stop_requested) exit
      ! TODO:
      !  - deposit charge from parts -> rho (grid)
      !  - solve Poisson -> phi (use dom/bnd)
      !  - E = -grad(phi)
      !  - push particles with (E, B from bkgB)
      !  - collisions / sources
      !  - diagnostics / saves
    end do

    ! ------------------------------------------------------------------
    ! Wrap up
    ! ------------------------------------------------------------------
    call timer_simulation%stop()
    if (rank_local==0) then
      write(output_unit,*) 'Elapsed time (s): ', timer_simulation%elapsed_s()
      call flush(output_unit)
    end if

    deallocate(b1, b2, b3)
    call MPI_Finalize(ierr_mpi)

  end subroutine SimulationPIC
end module mod_pic
